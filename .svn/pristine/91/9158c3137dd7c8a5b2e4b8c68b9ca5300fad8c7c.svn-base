#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>

#include <base/fe/Basis.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Element.hpp>
#include <base/dof/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/aux/algorithms.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/post/ErrorNorm.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/asmb/FieldIterator.hpp>
#include <heat/Laplace.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/FieldBinder.hpp>

//------------------------------------------------------------------------------
// Fundamental solution of the Laplace operator.
// This function is particularly useful, is it provides an exact solution for
// Laplace's equation for any geometry.
//[fSol]{
template<unsigned DIM>
class FundamentalSolution
{
public:
    typedef typename base::VectorType<DIM>::Type VecDim;
    typedef typename base::VectorType<1>::Type   VecDof;

    static VecDof evaluate( const VecDim& x )
    {
        const double dist = (sourcePoint_ - x).norm();

        const double aux =
            ( DIM == 2 ? -std::log( dist ) : 1./dist );

        VecDof result;
        result[0] = (1. / ((DIM-1) * M_PI) ) * aux;
        return result;
    }

private:
    static const VecDim sourcePoint_;
};

// initialise the source point (better outside of the domain!!)
template<unsigned DIM>
const typename FundamentalSolution<DIM>::VecDim
FundamentalSolution<DIM>::sourcePoint_ = base::constantVector<DIM>( -.5 );
//[fSol]}

//------------------------------------------------------------------------------
// Function for the point-wise constraint of the Boundary
//[dirichlet]{
template<unsigned DIM, typename DOF>
void dirichletBC( const typename FundamentalSolution<DIM>::VecDim& x,
                  DOF* doFPtr ) 
{
    const typename FundamentalSolution<DIM>::VecDof
        value = ( FundamentalSolution<DIM>::evaluate( x ) );

    if ( doFPtr -> isActive(0) ) {
        doFPtr -> constrainValue( 0, value[0] );
    }

}
//[dirichlet]}

struct PrintIDs
{
    template<typename FT>
    void operator()( const FT ft ) const
    {
        std::cout << (ft.template get<0>()) -> getID() <<" "
                  << (ft.template get<1>()) -> getID() <<" "
                  << (ft.template get<2>()) -> getID() <<"\n";
    }
};

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf \n\n";
        return -1;
    }
        
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = smfFile.substr( 0, smfFile.find( ".smf" ) );

    //--------------------------------------------------------------------------
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;
    
    const unsigned    doFSize = 1; // temperature

    //--------------------------------------------------------------------------
    const unsigned    dim     = base::ShapeDim<shape>::value;
    typedef base::mesh::Node<dim>                 Node;
    typedef base::LagrangeShapeFun<geomDeg,shape> SFun;
    typedef base::mesh::Element<Node,SFun>        Element;
    typedef base::mesh::Unstructured<Element>     Mesh;

    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::Reader<Mesh> smfReader;
        smfReader( mesh, smf ); 
        smf.close();
    }

    // Quadrature 
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // Finite element basis
    typedef base::fe::Basis<shape,fieldDeg> FEBasis;

    // DOF handling
    typedef base::dof::DegreeOfFreedom<doFSize>    DoF;
    typedef base::dof::Element<DoF,FEBasis::FEFun> FieldElement;
    typedef base::dof::Field<FieldElement>         Field;
    Field field;

    base::dof::generate<FEBasis>( mesh, field );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );
    
    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.boundaryBegin(),
                                           meshBoundary.boundaryEnd(),
                                           mesh, field, 
                                           boost::bind( &dirichletBC<dim,DoF>, _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << " Number of dofs " << numDofs << std::endl;
    
    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    typedef base::asmb::FieldBinder<Mesh,Field,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field, field );
    typedef FieldBinder::ElementPtrTuple FieldTuple;


    // Object providing the element stiffness matrix
    const double conductivity = 1.0;
    typedef heat::Laplace<FieldTuple> Laplace;
    Laplace laplace( conductivity );

    base::asmb::stiffnessMatrixComputation( quadrature, solver,
                                            fieldBinder, laplace );

    
    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.choleskySolve();

    // distribute results back to dofs
    {
        base::dof::Distribute<DoF,Solver> distributeDoF( solver );
        std::for_each( field.doFsBegin(), field.doFsEnd(), distributeDoF );
    }

    // output to a VTK file
    {
        // VTK Legacy
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        {
            // Evaluate the solution field at every geometry node
            std::vector<base::VectorType<doFSize>::Type> nodalValues;
            base::post::evaluateAtNodes( mesh, field, nodalValues );
            vtkWriter.writePointData( nodalValues.begin(), nodalValues.end(), "heat" );
        }
        vtk.close();
    }

    // compute L2-error and tell it to the user
    std::cout << "L2-error = "
              << base::post::errorComputation<0>(
                  quadrature, mesh, field,
                  boost::bind( &FundamentalSolution<dim>::evaluate, _1 ) )
              << '\n';

    return 0;
}
