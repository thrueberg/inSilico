#include <iostream>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>

#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
#include <base/mesh/ExtractMeshFaces.hpp>
#include <base/mesh/GenerateMeshFromFaces.hpp>
#include <base/mesh/FaceIterator.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/CreateBoundaryMesh.hpp>

#include <base/Quadrature.hpp>
#include <base/LagrangeShapeFun.hpp>

#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/xdmf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/gp/Writer.hpp>

#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/Handler.hpp>
#include <base/dof/ConstrainBoundary.hpp>
#include <base/dof/EvaluateAtNodes.hpp>
#include <base/dof/generate.hpp>

#include <base/dof/Element.hpp>

#include <base/StiffnessMatrix.hpp>
#include <base/ForceIntegrator.hpp>
#include <base/BodyForce.hpp>
#include <base/NeumannForce.hpp>

#include <base/EvaluateField.hpp>
#include <base/ErrorNorm.hpp>

#include <heat/Laplace.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/fe/Basis.hpp>

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

    typedef typename base::VectorType<1  >::Type   VecDof;
    typedef typename base::VectorType<DIM>::Type   VecDim;
    typedef typename base::MatrixType<DIM,1>::Type GradType;


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

        //result[ 0 ] = 0.;
        return result;
    }
    
private:
    const ReferenceSolution<DIM>& referenceSolution_;
};



//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;
    
    const unsigned    dofSize = 1;

    //--------------------------------------------------------------------------
    const unsigned    dim     = base::ShapeDim<shape>::value;
    typedef base::mesh::Node<dim>                 Node;
    typedef base::LagrangeShapeFun<geomDeg,shape> SFun;
    typedef base::mesh::Element<Node,SFun>        Element;
    typedef base::mesh::Unstructured<Element>     Mesh;

    Mesh mesh;

    const unsigned kernelDegEstimate = 4; //(fieldDeg > 1 ? (fieldDeg-1)*(fieldDeg-1) : 2 );
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    
    // I
    {
        const std::string smfFile = boost::lexical_cast<std::string>( argv[1] );
        std::ifstream smf( smfFile.c_str() );
        base::io::smf::Reader<Mesh> smfReader;
        smfReader( mesh, smf ); 
        smf.close();
    }

    // Finite element basis
    typedef base::fe::Basis<shape,fieldDeg> FEBasis;
    FEBasis feBasis;


    // DOFs i
    typedef base::dof::DegreeOfFreedom<dofSize>         DoF;
    typedef base::dof::Handler<DoF,FEBasis>      DoFHandler;


    typedef base::dof::Container<DoF,FEBasis>     Temperature;
    Temperature temperature;
    base::dof::generate( mesh, temperature );
    
    DoFHandler doFHandler( mesh.elementsBegin(), mesh.elementsEnd() );


    // sparsity pattern
    {
        base::dof::IndexMap<FEBasis> indexMap;
        indexMap.generateDoFIndices( mesh.elementsBegin(), mesh.elementsEnd() );

        std::vector< std::vector<std::size_t> > sparsity;
        indexMap.generateSparsityPattern( sparsity );
    
        
        std::ofstream sp( "sparsity" );
        for ( unsigned d1 = 0; d1 < sparsity.size(); d1 ++ ) {
            for ( unsigned d2 = 0; d2 < sparsity[d1].size(); d2++ ) {
                sp << d1 << " " << sparsity[d1][d2] << "\n";
            }
        }
        sp.close();
    }
    

    ReferenceSolution<dim> refSol( 3., 5., 4. );
    PoissonProblem<dim> pp( refSol );
    

    {
        typedef heat::Laplace<Element,FEBasis> Laplace;
        typedef base::solver::Eigen3           Solver;

        //! Creates a list of <Element,faceNo> pairs
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        //! Object to constrain the boundary 
        base::dof::ConstrainBoundary<DoFHandler> constrainBoundary( doFHandler );

        constrainBoundary.apply( meshBoundary.boundaryBegin(),
                                 meshBoundary.boundaryEnd(),
                                 mesh,
                                 boost::bind( &PoissonProblem<dim>::dirichleBC<DoF>,
                                              &pp, _1, _2 ) );

        //! Number of DoFs after constraint application!
        const std::size_t numDofs =
            base::dof::numberDoFsConsecutively( doFHandler.container().doFsBegin(),
                                                doFHandler.container().doFsEnd() );
        std::cout << " Number of dofs " << numDofs << std::endl;

        //! Create a solver object
        Solver solver( numDofs );

        //! Stiffness matrix kernel for Laplace operator
        Laplace laplace( 1., feBasis, feBasis );

        Mesh::ElementPtrConstIter elem     = mesh.elementsBegin();
        Mesh::ElementPtrConstIter lastElem = mesh.elementsEnd();

        //! Body force stuff
        typedef base::VectorType<dofSize>::Type VecDof;
        typedef boost::function< VecDof( const Node::VecDim& ) > ForceFun;
        typedef base::BodyForce<ForceFun, Element, FEBasis> BodyForce;
        ForceFun forceFunct =
            boost::bind( &PoissonProblem<dim>::forceFun, &pp, _1 );
        BodyForce bodyForce( forceFunct, feBasis  );

        //! Force integrator
        typedef base::ForceIntegrator<Element,Quadrature,Solver,DoFHandler> ForceInt;
        ForceInt::ForceKernel forceKernel =
            boost::bind( bodyForce, _1, _2, _3, _4 );
        ForceInt forceInt( forceKernel, quadrature, solver, doFHandler );
        
        std::for_each( elem, lastElem, forceInt );

        //! Neumann boundary condition

        //! Create list of <Element,faceNo> pairs
        typedef base::mesh::CreateBoundaryMesh<Element> CreateBoundaryMesh;
        typedef CreateBoundaryMesh::BoundaryMesh BoundaryMesh;
        BoundaryMesh boundaryMesh;
        {
            base::mesh::MeshBoundary meshBoundary;
            meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

            //! Create a connectivity out of this list
            CreateBoundaryMesh::apply( meshBoundary.boundaryBegin(),
                                       meshBoundary.boundaryEnd(),
                                       mesh,
                                       boundaryMesh );
        }

        typedef BoundaryMesh::Element SurfaceElement;
        typedef boost::function< VecDof( const Node::VecDim&,
                                         const Node::VecDim& ) > TractionFun;
        typedef base::NeumannForce<TractionFun,SurfaceElement,
                                   FEBasis>                      NeumannForce;
        TractionFun tractionFunct =
            boost::bind( &PoissonProblem<dim>::neumannBC, &pp, _1, _2 );
        NeumannForce neumannForce( tractionFunct, feBasis );

        typedef base::ForceIntegrator<SurfaceElement,SurfaceQuadrature,Solver,DoFHandler>
            SurfaceForceInt;

        SurfaceForceInt::ForceKernel surfaceForceKernel =
            boost::bind( neumannForce, _1, _2, _3, _4 );

        SurfaceForceInt surfaceForceInt( surfaceForceKernel, surfaceQuadrature,
                                         solver, doFHandler );
        
        std::for_each( boundaryMesh.elementsBegin(),
                       boundaryMesh.elementsEnd(),
                       surfaceForceInt );
        
        //! Compute element stiffness matrices and assemble them
        typedef base::StiffnessMatrix<Element,Quadrature,Solver,DoFHandler> StiffMat;
        StiffMat::Kernel kernel =
            boost::bind( &Laplace::stiffness, &laplace, _1, _2, _3, _4);

        StiffMat stiffness( kernel, quadrature, solver, doFHandler );
        std::for_each( elem, lastElem, stiffness );
        
        //! Finalise assembly
        solver.finishAssembly();

        //! Solve
        solver.choleskySolve();

        // distribute results back to dofs
        DoFHandler::Container::DoFPtrIter dof    = doFHandler.container().doFsBegin();
        DoFHandler::Container::DoFPtrIter dofEnd = doFHandler.container().doFsEnd();

        base::dof::Distribute<DoFHandler::DoF,Solver> distributeDoF( solver );
        std::for_each( dof, dofEnd, distributeDoF );

    }

    // OUTPUT
    {
        //! Evaluate the solution field at every geometry node
        typedef base::dof::EvaluateAtNodes<DoFHandler> Evaluator;
        std::vector<Evaluator::ValueType>
            nodalValues( std::distance( mesh.nodesBegin(),
                                        mesh.nodesEnd() ), Evaluator::ValueType() );
        Evaluator evaluator( doFHandler, feBasis );
        evaluator.apply( mesh.elementsBegin(), mesh.elementsEnd(), nodalValues );

        // XDMF
        {
            base::io::xdmf::Writer xdmfWriter;
            xdmfWriter.writeGeomtry( mesh, "coord" );
            xdmfWriter.writeTopolgy( mesh, "connec" );
            xdmfWriter.writeAttribute( nodalValues.begin(), nodalValues.end(),
                                       "heat", "Heat" );
            
            std::ofstream xdmf( "test.xdmf" );
            xdmfWriter.finalise( xdmf );
            xdmf.close();
        }

        // VTK Legacy
        {
            std::ofstream vtk( "test.vtk" );
            base::io::vtk::LegacyWriter vtkWriter( vtk );
            vtkWriter.writeUnstructuredGrid( mesh );
            vtkWriter.writePointData( nodalValues.begin(), nodalValues.end(), "heat" );
            vtk.close();
        }

        // GnuPlot text file
        {
            std::ofstream gp( "test.dat" );
            base::io::gp::Writer::apply( mesh.elementsBegin(), mesh.elementsEnd(), gp );
            gp.close();
        }
    }

    // Boundary
    {
        //! Create list of <Element,faceNo> pairs
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        //! Create a connectivity out of this list
        typedef base::mesh::CreateBoundaryMesh<Element> CreateBoundaryMesh;
        CreateBoundaryMesh::BoundaryMesh boundaryMesh;
        CreateBoundaryMesh::apply( meshBoundary.boundaryBegin(),
                                   meshBoundary.boundaryEnd(),
                                   mesh,
                                   boundaryMesh );

        const std::string smfFile = "bdry.smf";
        std::ofstream smf( smfFile.c_str() );
        base::io::smf::Writer<CreateBoundaryMesh::BoundaryMesh> smfWriter;
        smfWriter( boundaryMesh, smf );
        smf.close();

#if 0
        typedef base::mesh::ExtractMeshFaces<shape> MeshFaces;
        std::vector<MeshFaces::Face> faces;

        MeshFaces::boundary( mesh.elementsBegin(), mesh.elementsEnd(), faces );


        typedef base::mesh::GenerateMeshFromFaces<Mesh,MeshFaces::nFace> GMFF;
        GMFF::FaceMesh faceMesh;

        GMFF::apply( mesh, faces, faceMesh );

#endif
    }

    // compute L2-error
    {
        base::EvaluateField<DoFHandler,Element> ef( doFHandler, feBasis );

        typedef
            boost::function<
                ReferenceSolution<dim>::VecDof( const ReferenceSolution<dim>::VecDim )>
            ExactField;

        ExactField exF = boost::bind( &ReferenceSolution<dim>::evaluate, &refSol,
                                      _1 );

        typedef base::ErrorNorm< base::EvaluateField<DoFHandler,Element>,
                                 ExactField, Quadrature> L2Error;

        L2Error l2Error( ef, exF, quadrature );

        std::vector<double> l2Errors( std::distance( mesh.elementsBegin(),
                                                     mesh.elementsEnd() ) );
        
        std::transform( mesh.elementsBegin(), mesh.elementsEnd(),
                        l2Errors.begin(), l2Error );

        const double l2ErrorSquared
            = std::accumulate( l2Errors.begin(), l2Errors.end(), 0. );

        std::cout << "L2-error = " << std::sqrt( l2ErrorSquared ) << "\n";
    }

    // compute H1-error
    {
        base::EvaluateFieldGradient<DoFHandler,Element> efg( doFHandler, feBasis );

        typedef
            boost::function<
                ReferenceSolution<dim>::GradType( const ReferenceSolution<dim>::VecDim )>
            ExactGradient;


        ExactGradient exG = boost::bind( &ReferenceSolution<dim>::evaluateGradient, &refSol,
                                         _1 );
        
        typedef base::ErrorNorm< base::EvaluateFieldGradient<DoFHandler,Element>,
                                 ExactGradient, Quadrature> H1Error;

        H1Error h1Error( efg, exG, quadrature );

        std::vector<double> h1Errors( std::distance( mesh.elementsBegin(),
                                                     mesh.elementsEnd() ) );
        
        std::transform( mesh.elementsBegin(), mesh.elementsEnd(),
                        h1Errors.begin(), h1Error );

        const double h1ErrorSquared
            = std::accumulate( h1Errors.begin(), h1Errors.end(), 0. );

        std::cout << "H1-error = " << std::sqrt( h1ErrorSquared ) << "\n";
    }
    
    return 0;
}
