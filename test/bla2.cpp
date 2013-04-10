#include <iostream>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>

#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Structured.hpp>

#include <base/Quadrature.hpp>

#include <base/sfun/BSpline.hpp>
#include <base/sfun/TensorProduct.hpp>

#include <base/io/sgf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/PropertiesParser.hpp>

#include <base/dof/access.hpp>
#include <base/dof/Container.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>

#include <base/ComputeAndAssembleMatrix.hpp>

#include <heat/Laplace.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/io/xdmf/Writer.hpp>

#include <base/io/vtk/LegacyWriter.hpp>

#include <base/mesh/SampleStructured.hpp>
#include <base/mesh/SurfaceElement.hpp>
#include <base/mesh/GenerateBoundaryTriangulationFromBox.hpp>

#include <base/io/OStreamIterator.hpp>

#include <base/EvaluateField.hpp>

//------------------------------------------------------------------------------
template<typename VEC>
struct WriteVector : boost::function< void( const VEC&,
                                            std::ostream & )>
{
    void operator()( const VEC & vec,
                     std::ostream & out ) const
    {
        for ( int d = 0; d < vec.size(); d ++ )
            out << vec[d] << " ";
        out << '\n';
    }
};

//------------------------------------------------------------------------------
template<typename VEC>
void writeVec( const VEC & vec,
               std::ostream & out )
{
    for ( int d = 0; d < vec.size(); d ++ )
        out << vec[d] << " ";
    out << '\n';
}

//------------------------------------------------------------------------------
template<typename NODE, typename DOFACC>
struct FixBoundaryValues
{
    FixBoundaryValues( DOFACC & da )
    : dofAccessor( da ) { }
    
    void operator()( const NODE* np ) const
    {
        typename DOFACC::DegreeOfFreedom * dof =
            dofAccessor.getDoFObject( np );

        typename NODE::VecDim x = np -> getX();

        const double tol = 1.e-5;
        
        const bool isOnBoundary =
            ( std::abs( x[0] - 0. ) < tol ) or
            ( std::abs( x[0] - 1. ) < tol ) or
            ( std::abs( x[1] - 0. ) < tol ) or
            ( std::abs( x[1] - 1. ) < tol );

        if ( isOnBoundary ) {
            typename NODE::VecDim xStar = NODE::VecDim::Constant( -1. );

            const double r = (x - xStar).norm();
            dof -> constrainValue( 0, std::log( r ) );

        }
    }

    DOFACC & dofAccessor;
};


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    const unsigned    dim     = 3;
    const unsigned    sfunDeg = 1;
    const base::Shape shape   = base::HyperCubeShape<dim>::value;
    const unsigned    dofSize = 1;
    
    typedef base::mesh::Node<dim>                 Node;
    //typedef base::LagrangeShapeFun<sfunDeg,shape> SFun;
    typedef base::sfun::BSpline<sfunDeg>          SFun1D;
    typedef base::sfun::TensorProduct<SFun1D,dim> SFun;
    typedef base::mesh::Element<Node,shape,SFun>  Element;
    typedef base::mesh::Structured<Element>       Mesh;

    Mesh mesh;

    typedef base::Quadrature<3,shape> Quadrature;
    Quadrature quadrature;


    // I
    {
        const std::string sgfFile = boost::lexical_cast<std::string>( argv[1] );
        std::ifstream sgf( sgfFile.c_str() );
        base::io::sgf::Reader<Mesh> sgfReader;
        sgfReader( mesh, sgf ); 
        sgf.close();
    }

    // O
    //{
    //    const std::string smfFile = "test.smf";
    //    std::ofstream smf( smfFile.c_str() );
    //    base::io::smf::Writer<Mesh> smfWriter;
    //    smfWriter( mesh, smf );
    //    smf.close();
    //}

    
    // DOFS
    typedef base::dof::DegreeOfFreedom<dofSize>         DoF;
    typedef base::dof::Container<DoF>                   DoFContainer;
    typedef base::dof::AccessViaNode<DoFContainer>      DoFAccessor;

    DoFContainer dofs;

    //! Generate a container full of dofs
    base::dof::GenerateFromNodes<Mesh,DoFContainer>()( mesh, dofs );

    DoFAccessor dofAccessor( dofs );
    {
        typedef heat::Laplace<Element> Laplace;
        typedef base::solver::Eigen3   Solver;



        FixBoundaryValues<Node,DoFAccessor> fbv( dofAccessor );

        std::for_each( mesh.nodesBegin(), mesh.nodesEnd(), fbv );

        //! Number them according to the nodes
        const std::size_t numDofs = 
            base::dof::NumberDoFsConsecutively<DoFContainer>()( dofs );
        //std::cout << "Number of Dofs = " << numDofs << std::endl;

        Solver solver( numDofs );

        Laplace laplace( 1. );

        //base::ComputeAndAssembleMatrix<Laplace,Quadrature,Solver,DoFAccessor,DoFAccessor>
        //    computeAndAssemble( laplace, quadrature, solver, dofAccessor, dofAccessor );

        base::ComputeAndAssembleMatrix<Laplace,Quadrature,Solver,DoFAccessor>
            computeAndAssemble( laplace, quadrature, solver, dofAccessor );

        Mesh::ElementPtrConstIter elem     = mesh.elementsBegin();
        Mesh::ElementPtrConstIter lastElem = mesh.elementsEnd();

        std::for_each( elem, lastElem,  computeAndAssemble );
        
        solver.finishAssembly();
        
        // solve
        solver.choleskySolve();

        // distribute results
        DoFContainer::DoFPtrIter dof    = dofs.dofsBegin();
        DoFContainer::DoFPtrIter dofEnd = dofs.dofsEnd();

        base::dof::Distribute<DoF,Solver> distributeDoF( solver );
        std::for_each( dof, dofEnd, distributeDoF );
    }

    // XDMF
    {
        base::io::xdmf::Writer xdmfWriter;
        xdmfWriter.writeGeomtry( mesh, "coords" );
        xdmfWriter.writeTopolgy( mesh, "connec" );
        xdmfWriter.writeAttribute( dofs, "heat", "Heat" );

        std::ofstream xdmf( "test.xdmf" );
        xdmfWriter.finalise( xdmf );
        xdmf.close();
    }

    // VTK Legacy
    {
        std::ofstream vtk( "test.vtk" );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        vtkWriter.writePointData( dofs, "heat" );
        vtk.close();
    }

    // try sampling
    {
        // Function object for evaluating the element geometry
        typedef base::GeomTraits<Element>::LocalVecDim LVD;
        typedef boost::function< Node::VecDim( const Element*,
                                               const LVD&) > GeomEval;
        GeomEval geometry = boost::bind( base::Geometry<Element>(), _1, _2 );

        // structured grid sampling function
        base::mesh::SampleStructured<Mesh,GeomEval> ss;

        // sample structured grid with sub-sampling
        const unsigned resolution = 2;
        ss( mesh, geometry,
            base::io::OStreamIterator< LVD >( std::cout, &writeVec<LVD> ),
            resolution );
    }

    // surface shit
    {
        typedef base::mesh::GenerateBoundaryTriangulationFromBox<Mesh> GBTFBB;
        GBTFBB::SurfaceMesh surfaceMesh;

        GBTFBB()( mesh, 0, true, surfaceMesh );

        //write smf
        const std::string smfFile = "testSurf.smf";
        std::ofstream smf( smfFile.c_str() );
        base::io::smf::Writer<GBTFBB::SurfaceMesh> smfWriter;
        smfWriter( surfaceMesh, smf );
        smf.close();


    }

    {
        SFun shapeFun;
        typedef base::EvaluateField<Element,DoFAccessor,SFun> EF;
        typedef EF::FieldValueType DoFVal;
        EF::VecDim xi = base::constantVector<dim>( 0.5 );

        EF ef( dofAccessor, shapeFun );

        const std::size_t numElements
            = std::distance( mesh.elementsBegin(), mesh.elementsEnd() );
        std::vector<DoFVal> fieldValues( numElements );

        boost::function<DoFVal( const Element* )> fieldOp
            = boost::bind( &EF::operator(), &ef, _1, xi );

        std::transform( mesh.elementsBegin(), mesh.elementsEnd(),
                        fieldValues.begin(), fieldOp );

        std::copy( fieldValues.begin(), fieldValues.end(),
                   base::io::OStreamIterator<DoFVal>( std::cout,
                                                      &writeVec<DoFVal> ) );
        
    }
    
    return 0;
}
