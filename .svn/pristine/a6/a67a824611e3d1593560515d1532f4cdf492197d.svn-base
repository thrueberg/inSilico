#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/CreateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/io/Format.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Element.hpp>
#include <base/dof/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/post/ErrorNorm.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <heat/Laplace.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/fe/Basis.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>


#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/NeumannForce.hpp>


//------------------------------------------------------------------------------
template<unsigned DIM>
class ReferenceSolution
{
public:
    
    ReferenceSolution( const double alpha,
                       const double beta  = 0.,
                       const double gamma = 0. )
        : alpha_( alpha ),
          beta_(  beta ),
          gamma_( gamma )
    { }

    typedef typename base::Vector<1  >::Type   VecDof;
    typedef typename base::Vector<DIM>::Type   VecDim;
    typedef typename base::Matrix<DIM,1>::Type GradType;


    VecDof evaluate( const VecDim& x ) const
    {
        VecDof result;
        result[0] = std::cos( alpha_ * x[0] );
        if ( DIM > 1 )
            result[0] *= std::cos( beta_ * x[1] );
        if ( DIM > 2 )
            result[0] *= std::cos( gamma_ * x[2] );

        return result;
    }

    GradType evaluateGradient( const VecDim& x ) const
    {
        GradType result;
        result( 0, 0 ) = -alpha_ * std::sin( alpha_ * x[0] );
        if ( DIM > 1 ) {
            result( 0, 0 ) *= std::cos( beta_ * x[1] );
            result( 1, 0 )  = -beta_ *
                std::cos( alpha_ * x[0] ) *
                std::sin( beta_ *  x[1] );
        }
        if ( DIM > 2 ) {
            result( 0, 0 ) *= std::cos( gamma_ * x[2] );
            result( 1, 0 ) *= std::cos( gamma_ * x[2] );
            result( 2, 0 ) = -gamma_ *
                std::cos( alpha_ * x[0] ) *
                std::cos( beta_  * x[1] ) *
                std::sin( gamma_ * x[2] );
        }

        return result;
    }

    VecDof laplacian( const VecDim& x ) const
    {
        const double factor = -1.0 *
            ( ( alpha_ * alpha_ ) +
              ( DIM > 1 ? beta_  * beta_  : 0.0 ) +
              ( DIM > 2 ? gamma_ * gamma_ : 0.0 ) );

        return factor * evaluate( x );
    }

private:
    const double alpha_;
    const double beta_;
    const double gamma_;
};


//------------------------------------------------------------------------------
template<unsigned DIM>
class PoissonProblem
{
public:
    PoissonProblem( const ReferenceSolution<DIM>& rs )
        : referenceSolution_( rs ) { }

    template<typename DOF>
    void dirichleBC( const typename ReferenceSolution<DIM>::VecDim& x,
                     DOF* doFPtr ) const
    {
        const double value = ( referenceSolution_.evaluate( x ) )[0];

        const double tol = 1.e-5;

        // make left and bottom boundaries as Dirichlet
        const bool onDirichletBdr = //true;
            ( std::abs( x[0] -  0. ) < tol ) or
            ( std::abs( x[1] -  0. ) < tol );

        if ( onDirichletBdr ) {
            if ( doFPtr -> isActive(0) )
                doFPtr -> constrainValue( 0, value );
        }
        return;
    }

    typename ReferenceSolution<DIM>::VecDof
    forceFun( const typename ReferenceSolution<DIM>::VecDim& x ) const
    {
        typename ReferenceSolution<DIM>::VecDof result;
        result = - referenceSolution_.laplacian( x );
        return result;
    }

    typename ReferenceSolution<DIM>::VecDof
    neumannBC( const typename ReferenceSolution<DIM>::VecDim& x,
               const typename ReferenceSolution<DIM>::VecDim& normal ) const
    {
        const typename ReferenceSolution<DIM>::GradType gradient
            = referenceSolution_.evaluateGradient( x );

        typename ReferenceSolution<DIM>::VecDof
            result = gradient.transpose() * normal;

        return result;
    }
    
private:
    const ReferenceSolution<DIM>& referenceSolution_;
};



//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf \n\n";
        return -1;
    }

    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    //--------------------------------------------------------------------------
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;
    const unsigned    doFSize  = 1;

    //--------------------------------------------------------------------------
    const unsigned    dim     = base::ShapeDim<shape>::value;
    typedef base::mesh::Node<dim>                 Node;
    typedef base::LagrangeShapeFun<geomDeg,shape> SFun;
    typedef base::mesh::Element<Node,SFun>        Element;
    typedef base::mesh::Unstructured<Element>     Mesh;

    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        base::io::smf::Reader<Mesh> smfReader;
        smfReader( mesh, smf ); 
        smf.close();
    }

    // Quadrature and surface quadrature
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;
 
    // DOF handling
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::dof::DegreeOfFreedom<doFSize>    DoF;
    typedef base::dof::Element<DoF,FEBasis::FEFun> FieldElement;
    typedef base::dof::Field<FieldElement>         Field;
    Field field;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, field );

    // Object of reference solution and the Poisson problem
    ReferenceSolution<dim> refSol( 3., 5., 4. );
    PoissonProblem<dim> pp( refSol );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.boundaryBegin(),
                                           meshBoundary.boundaryEnd(),
                                           mesh, field,
                                           boost::bind( &PoissonProblem<dim>::dirichleBC<DoF>,
                                                        &pp, _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << " Number of dofs " << numDofs << std::endl;

    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field, field );
    typedef FieldBinder::ElementPtrTuple FieldTuple;


    // Body force
    base::asmb::bodyForceComputation( quadrature, solver, fieldBinder,
                                      boost::bind( &PoissonProblem<dim>::forceFun,
                                                   &pp, _1 ) );

    // Create a connectivity out of this list
    typedef base::mesh::CreateBoundaryMesh<Element> CreateBoundaryMesh;
    typedef CreateBoundaryMesh::BoundaryMesh BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        CreateBoundaryMesh::apply( meshBoundary.boundaryBegin(),
                                   meshBoundary.boundaryEnd(),
                                   mesh,
                                   boundaryMesh );
    }


    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, field );

    // Neumann boundary condition
    base::asmb::neumannForceComputation( surfaceQuadrature, solver, surfaceFieldBinder,
                                         boost::bind( &PoissonProblem<dim>::neumannBC, &pp,
                                                      _1, _2 ) );


    // compute stiffness matrix
    typedef heat::Laplace<FieldTuple> Laplace;
    Laplace laplace( 1. );
    base::asmb::stiffnessMatrixComputation( quadrature, solver,
                                            fieldBinder, laplace );
    
    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.choleskySolve();

    // distribute results back to dofs
    base::dof::Distribute<DoF,Solver> distributeDoF( solver );
    std::for_each( field.doFsBegin(), field.doFsEnd(), distributeDoF );

    //--------------------------------------------------------------------------
    // output to a VTK file
    {
        // VTK Legacy
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        base::io::vtk::writePointData( vtkWriter, mesh, field, "temperature" );
        vtk.close();
    }

    //--------------------------------------------------------------------------
    // compute L2-error
    std::cout << "L2-error = "
              << base::post::errorComputation<0>(
                  quadrature, mesh, field,
                  boost::bind( &ReferenceSolution<dim>::evaluate,
                               &refSol, _1 ) )
              << '\n';

    //--------------------------------------------------------------------------
    // compute H1-error
    std::cout << "H1-error = "
              << base::post::errorComputation<1>(
                  quadrature, mesh, field,
                  boost::bind( &ReferenceSolution<dim>::evaluateGradient,
                               &refSol, _1 ) )
              << '\n';

    
    return 0;
}
