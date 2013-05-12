#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
#include <boost/lexical_cast.hpp> // lexical cast between objects

//------------------------------------------------------------------------------
#include <base/shape.hpp>
#include <base/BSplineShapeFun.hpp>
#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Structured.hpp>
#include <base/mesh/sampleStructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/CreateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>

#include <base/io/sgf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/fe/Basis.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Element.hpp>
#include <base/dof/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/post/evaluateAtNodes.hpp>
#include <base/post/ErrorNorm.hpp>

#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>

#include <base/kernel/Laplace.hpp>

#if 0
//------------------------------------------------------------------------------
//  Academic example:
//  The solution function
//  u(x) = \prod_{d=1}^DIM (x_d^4 - 2 x_d^3 + x_d^2)
//  has zero normal derivatives along the boundary of (0,1)^DIM
//  and the yields the forcing function
//  f(x) = -u''(x)
double fun(   const double x ) { return    x*x*x*x -  2.*x*x*x +    x*x;}
double funD(  const double x ) { return  4.*x*x*x  -  6.*x*x   + 2.*x;  }
double funDD( const double x ) { return 12.*x*x    - 12.*x     + 2.;    }
#else

const double alpha = 5.;
double fun(   const double x ) { return                std::sin( alpha * x ); }
double funD(  const double x ) { return        alpha * std::cos( alpha * x ); }
double funDD( const double x ) { return -alpha*alpha * std::sin( alpha * x ); }
#endif

//------------------------------------------------------------------------------
// Tensor-product realisations of 1D solution function
template<unsigned DIM>
base::Vector<1>::Type 
solution( const typename base::Vector<DIM>::Type& x )
{
    base::Vector<1>::Type u;
    u[0] = 1.;

    for ( unsigned d = 0; d < DIM; d++ ) u[0] *= fun( x[d] );

    return u;
}


template<unsigned DIM>
typename base::Matrix<DIM,1>::Type
gradient( const typename base::Vector<DIM>::Type& x )
{
    typename base::Matrix<DIM,1>::Type gradU;

    for ( unsigned d = 0; d < DIM; d++ ) {
        gradU[d] = funD( x[d] );

        for ( unsigned d2 = 0; d2 < DIM; d2++ ) {
            if ( d2 != d ) gradU[d] *= fun( x[d2] );
        }
    }

    return gradU;
}

template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM>::Type& x,
                  DOF* doFPtr ) 
{
    const typename base::Vector<1>::Type value = solution<DIM>( x );

    if ( doFPtr -> isActive(0) ) {
        doFPtr -> constrainValue( 0, value[0] );
    }

}

template<unsigned DIM>
base::Vector<1>::Type 
neumannBC( const typename base::Vector<DIM>::Type& x,
           const typename base::Vector<DIM>::Type& normal )
{
    const typename base::Vector<DIM>::Type gradU = gradient( x );
    return (gradU.transpose() * normal );
}


template<unsigned DIM>
base::Vector<1>::Type
forceFun( const typename base::Vector<DIM>::Type& x ) 
{
    base::Vector<1>::Type f;

    f[0] = 0.;

    for ( unsigned d = 0; d < DIM; d++ ) {
        double aux = funDD( x[d] );
        for ( unsigned d2 = 0; d2 < DIM; d2++ ) {
            if ( d2 != d ) aux *= fun( x[d2] );
        }
        f[0] += aux;
    }

    return -f;
}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // Check the number of input arguments
    if ( argc != 2 ) { 
        std::cout << "Usage:  " << argv[0] << " file.sgf \n";
        return 0;
    }
    
    // convert input argument to a string object
    const std::string sgfFileName = boost::lexical_cast<std::string>( argv[1] );
    // extract the basename of the file
    const std::string baseName = sgfFileName.substr( 0, sgfFileName.find( ".sgf" ) );

    //--------------------------------------------------------------------------
    const double penaltyFactor = 100.;
    const bool   nitsche = true;

    const unsigned dim      = 2;
    const unsigned geomDeg  = 1;
    const unsigned fieldDeg = 1;
    const unsigned doFSize  = 1;
    
    typedef base::mesh::Node<dim>                        Node;
    typedef base::BSplineShapeFun<dim,geomDeg>           GeomFun;
    typedef base::mesh::Element<Node,GeomFun>            Element;
    typedef base::mesh::Structured<Element>              Mesh;
    Mesh mesh;

    //--------------------------------------------------------------------------
    {
        std::ifstream sgf( sgfFileName.c_str() );
        base::io::sgf::Reader<Mesh> sgfReader;
        sgfReader( mesh, sgf );
        sgf.close();
    }


    //--------------------------------------------------------------------------
    // Quadrature 
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,Element::shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,Element::shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // Finite element basis
    typedef base::fe::Basis<Element::shape,fieldDeg,base::BSPLINE> FEBasis;
    
    // DOF handling
    typedef base::dof::DegreeOfFreedom<doFSize>    DoF;
    typedef base::dof::Element<DoF,FEBasis::FEFun> FieldElement;
    typedef base::dof::Field<FieldElement>         Field;
    Field field;
    base::dof::generate<FEBasis>( mesh, field );

    // fix a dof
    (*field.doFsBegin()) -> constrainValue( 0, 0. );

    // Number of DoFs 
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    //std::cout << " Number of dofs " << numDofs << std::endl;

    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field, field );
    typedef FieldBinder::ElementPtrTuple FieldTuple;

    // Body force
    base::asmb::bodyForceComputation( quadrature, solver, fieldBinder,
                                      boost::bind( &forceFun<dim>, _1 ) );

    //--------------------------------------------------------------------------
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create<dim>( mesh.gridSizes() );

    // Create a connectivity out of this list
    typedef base::mesh::CreateBoundaryMesh<Element> CreateBoundaryMesh;
    typedef CreateBoundaryMesh::BoundaryMesh BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        CreateBoundaryMesh::apply( meshBoundary.boundaryBegin(),
                                   meshBoundary.boundaryEnd(),
                                   mesh, boundaryMesh );
    }

    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, field, field );

#if 0
    // Neumann boundary condition
    base::asmb::neumannForceComputation( surfaceQuadrature, solver, surfaceFieldBinder,
                                         boost::bind( neumannBC<dim>, _1, _2 ) );
#endif

    //--------------------------------------------------------------------------
    // Stiffness matrix
    typedef base::kernel::Laplace<FieldTuple> Laplace;
    Laplace laplace( 1. );
    base::asmb::stiffnessMatrixComputation( quadrature, solver,
                                            fieldBinder, laplace );

    
    // Dirichlet BCs a la Nitsche
    base::nitsche::penaltyLHS( surfaceQuadrature, solver, surfaceFieldBinder,
                               penaltyFactor );
    
    base::nitsche::penaltyRHS( surfaceQuadrature, solver, surfaceFieldBinder,
                               boost::bind( solution<dim>, _1 ), penaltyFactor );

    if ( nitsche ) {
    
        base::nitsche::energyLHS( laplace, surfaceQuadrature, solver,
                                  surfaceFieldBinder );

        base::nitsche::energyRHS( laplace, surfaceQuadrature, solver, surfaceFieldBinder,
                                  boost::bind( solution<dim>, _1 ) );
    }


    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.choleskySolve();

    // distribute results back to dofs
    base::dof::Distribute<DoF,Solver> distributeDoF( solver );
    std::for_each( field.doFsBegin(), field.doFsEnd(), distributeDoF );

    //--------------------------------------------------------------------------
    // VTK file of the input mesh
    {
        const std::string vtkMeshFileName = baseName + ".vtk";
        std::ofstream vtk( vtkMeshFileName.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeStructuredGrid( mesh );
        {
            // Evaluate the solution field at every geometry node
            std::vector<base::Vector<doFSize>::Type> nodalValues;
            typedef base::GeomTraits<Element>::LocalVecDim LVD;

            base::post::sampleField( mesh, field, 
                                     std::back_inserter( nodalValues ) );
            vtkWriter.writePointData( nodalValues.begin(), nodalValues.end(), "heat" );
        }

        {
            // Evaluate the solution field at every geometry node
            std::vector<base::Matrix<dim,doFSize>::Type> cellValues;
            typedef base::GeomTraits<Element>::LocalVecDim LVD;

            base::post::sampleFieldGradient( mesh, field, 
                                             std::back_inserter( cellValues ), 0 );
            vtkWriter.writeCellData( cellValues.begin(), cellValues.end(), "flux" );
        }

        vtk.close();
    }
    
    // SMF file of the boundary mesh
    {
        // file name for boundary mesh output
        const std::string smfBoundaryFileName = baseName + "_boundary.smf";
        // output stream
        std::ofstream smf( smfBoundaryFileName.c_str() );
        // smf-writer object
        base::io::smf::Writer<CreateBoundaryMesh::BoundaryMesh> smfWriter;
        // write mesh to stream
        smfWriter( boundaryMesh, smf );
        // close stream
        smf.close();
    }

    // compute L2-error and tell it to the user
    std::cout //<< "L2-error = "
              << base::post::errorComputation<0>(
                  quadrature, mesh, field,
                  boost::bind( &solution<dim>, _1 ) )
              << "  " 
              << base::post::errorComputation<1>(
                  quadrature, mesh, field,
                  boost::bind( &gradient<dim>, _1 ) )
              << '\n';
    
    return 0;
}
