#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/CreateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <base/dof/scaleConstraints.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/setField.hpp>
#include <base/fe/Basis.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include <base/asmb/BodyForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/Lame.hpp>

#include <solid/HyperElastic.hpp>
#include <fluid/Stokes.hpp>
#include <heat/Laplace.hpp>

#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

#include "Terzaghi.hpp"


const double coordTol = 1.e-5; // tolerance for coordinate comparison

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD1, typename FIELD2>
void writeVTK( const MESH& mesh,
               const FIELD1& displacement,
               const FIELD2& pressure,
               const std::string& meshFile,
               const unsigned numStep )
{
    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );
    // create file name with step number
    const std::string vtkFile = baseName + "." +
        base::io::leadingZeros( numStep ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );

    base::io::vtk::writePointData( vtkWriter, mesh, displacement, "disp" );
    base::io::vtk::writePointData( vtkWriter, mesh, pressure,     "pressure" );
    vtk.close();
}


//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void dirichletU( const typename base::Vector<DIM>::Type& x, DOF* doFPtr,
                 const double width )
{
    // fix bottom boundary
    if ( std::abs( x[1] - 0. ) < coordTol ) {
        for ( unsigned d = 0; d < DIM; d++ ) 
            doFPtr -> constrainValue( d, 0.0 );
    }

    // constraint in normal direction of the lateral boundaries
    if ( ( std::abs( x[0] - 0.    ) < coordTol ) or
         ( std::abs( x[0] - width )  < coordTol ) ) {
        doFPtr -> constrainValue( 0, 0.0 );
    }
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void dirichletP( const typename base::Vector<DIM>::Type& x, DOF* doFPtr,
                 const double height )
{
    // zero pressure on top boundary
    if ( std::abs( x[1] - height ) < coordTol ) {
        doFPtr -> constrainValue( 0, 0.0 );
    }
}


//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM>::Type
neumannBC( const typename base::Vector<DIM>::Type& x,
           const typename base::Vector<DIM>::Type& normal,
           const double height, const double force )
{
    typename base::Vector<DIM>::Type result = base::constantVector<DIM>( 0. );
    if ( std::abs( x[1] - height ) < coordTol ) {
        result[1] = -force;
    }
    
    return result;
}

//------------------------------------------------------------------------------
/** Overload the div U term in order to multiply it by the Biot coefficient
 */
template<typename FIELDTUPLE>
class PoroDivU
    : public fluid::VelocityDivergence<FIELDTUPLE>
{
    typedef fluid::VelocityDivergence<FIELDTUPLE> Basis;
    
public:
    PoroDivU( const double alpha)
    : Basis( true ),
      alpha_( alpha ) { }
    
    void tangentStiffness( const typename Basis::FieldTuple&     fieldTuple,
                           const typename Basis::LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        Basis::tangentStiffness( fieldTuple, xi, weight * alpha_, matrix );
    }

private:
    const double alpha_;
};


//------------------------------------------------------------------------------
template<typename MESH, typename DISP, typename PRESS>
std::pair<double,double> probe( const MESH& mesh, const DISP& displacement,
                                const PRESS& pressure )
{
    // pick an element
    const std::size_t numElements = std::distance( mesh.elementsBegin(), mesh.elementsEnd() );
    const std::size_t elemNum = (numElements +
                                 static_cast<std::size_t>(std::sqrt(numElements)) ) / 2;

    // geometry
    const typename MESH::Element* geomEp = mesh.elementPtr( elemNum );
    const typename base::Vector<MESH::Element::dim>::Type xi =
        base::constantVector<MESH::Element::dim>( 0. );

    const typename base::Geometry<typename MESH::Element>::result_type x =
        base::Geometry<typename MESH::Element>()( geomEp, xi );

    // displacement
    const typename DISP::Element* dispEp = displacement.elementPtr( elemNum );
    const typename base::Vector<DISP::DegreeOfFreedom::size>::Type uh =
        base::post::evaluateField( geomEp, dispEp, xi );
    
    // pressure
    const typename PRESS::Element* pressEp = pressure.elementPtr( elemNum );
    const typename base::Vector<PRESS::DegreeOfFreedom::size>::Type ph =
        base::post::evaluateField( geomEp, pressEp, xi );

    return std::make_pair( uh[1], ph[0] );
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // basic attributes of the computation
    const unsigned    geomDeg   = 1;
    const unsigned    fieldDegU = 2;
    const unsigned    fieldDegP = 1;
    const unsigned    tiOrder   = 2;   // order of time integrator
    const base::Shape shape     = base::QUAD;

    //typedef  base::time::BDF<tiOrder> MSM;
    typedef base::time::AdamsMoulton<tiOrder> MSM;

    // time stepping method determines the history size
    const unsigned nHist = MSM::numSteps;

    
    const std::string inputFile = "terz.inp";
    double E, nu, alpha, c0, k, H, F, tMax;
    unsigned numSteps, numIntervals;
    {
        base::io::PropertiesParser pp;
        pp.registerPropertiesVar( "E",     E     );
        pp.registerPropertiesVar( "nu",    nu    );
        pp.registerPropertiesVar( "alpha", alpha );
        pp.registerPropertiesVar( "c0",    c0    );
        pp.registerPropertiesVar( "k",     k     );
        pp.registerPropertiesVar( "H",     H     );
        pp.registerPropertiesVar( "F",     F     );

        pp.registerPropertiesVar( "tMax",         tMax  );
        pp.registerPropertiesVar( "numSteps",     numSteps );
        pp.registerPropertiesVar( "numIntervals", numIntervals );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str() );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        pp.readValues( inp );
        inp.close( );

    }

    const double deltaT = tMax / static_cast<double>( numSteps );

    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  mesh.smf \n";
        return 0;
    }

    const std::string meshFile = boost::lexical_cast<std::string>( argv[1] );

    

    // Analytic solution
    Terzaghi terzaghi( E, nu, alpha, c0, k, H, F );
    

    //--------------------------------------------------------------------------
    // define a mesh
    typedef base::Unstructured<shape,geomDeg>    Mesh;
    const unsigned dim = Mesh::Node::dim;

    // create a mesh and read from input
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // quadrature objects for volume and surface
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // Create a displacement field
    const unsigned    doFSizeU = dim;
    typedef base::fe::Basis<shape,fieldDegU>           FEBasisU;
    typedef base::Field<FEBasisU,doFSizeU,nHist>       Displacement;
    typedef Displacement::DegreeOfFreedom              DoFU;
    Displacement displacement;

    // Create a pressure field
    const unsigned    doFSizeP = 1;
    typedef base::fe::Basis<shape,fieldDegP>           FEBasisP;
    typedef base::Field<FEBasisP,doFSizeP,nHist>       Pressure;
    typedef Pressure::DegreeOfFreedom                  DoFP;
    Pressure pressure;

    // generate DoFs from mesh
    base::dof::generate<FEBasisU>( mesh, displacement );
    base::dof::generate<FEBasisP>( mesh, pressure );

    // Creates a list of <Element,faceNo> pairs along the boundary
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Create a boundary mesh from this list
    typedef base::mesh::CreateBoundaryMesh<Mesh::Element> CreateBoundaryMesh;
    typedef CreateBoundaryMesh::BoundaryMesh BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        CreateBoundaryMesh::apply( meshBoundary.begin(),
                                   meshBoundary.end(),
                                   mesh, boundaryMesh );
    }



    // material object
    typedef mat::hypel::StVenant Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    typedef base::asmb::FieldBinder<Mesh,Displacement,Pressure> Field;
    Field field( mesh, displacement, pressure );

    typedef Field::TupleBinder<1,1>::Type TopLeft;
    typedef Field::TupleBinder<1,2>::Type TopRight;
    typedef Field::TupleBinder<2,1>::Type BotLeft;
    typedef Field::TupleBinder<2,2>::Type BotRight;

    // surface displacement field
    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Displacement> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, displacement );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SFTB;
    
    // kernel objects
    typedef solid::HyperElastic<Material,TopLeft::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    typedef fluid::PressureGradient<TopRight::Tuple> GradP;
    GradP gradP;

    //typedef fluid::VelocityDivergence<FieldPUUTuple> DivU;
    //DivU divU( true ); // * alpha
    typedef PoroDivU<BotLeft::Tuple> DivU;
    DivU divU( alpha );

    typedef heat::Laplace<BotRight::Tuple> Laplace;
    Laplace laplace( k );

    // constrain the boundary
    base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                            meshBoundary.end(),
                                            mesh, displacement, 
                                            boost::bind( &dirichletU<dim,DoFU>, _1, _2, 1.0 ) );
    
    base::dof::constrainBoundary<FEBasisP>( meshBoundary.begin(),
                                            meshBoundary.end(),
                                            mesh, pressure, 
                                            boost::bind( &dirichletP<dim,DoFP>, _1, _2, H ) );


    // Number the degrees of freedom
    const std::size_t numDoFsU =
        base::dof::numberDoFsConsecutively( displacement.doFsBegin(), displacement.doFsEnd() );
    std::cout << "# Number of displacement dofs " << numDoFsU << std::endl;
    const std::size_t numDoFsP =
        base::dof::numberDoFsConsecutively( pressure.doFsBegin(), pressure.doFsEnd(), numDoFsU );
    std::cout << "# Number of pressure     dofs " << numDoFsP << std::endl;


    // write VTK file
    writeVTK( mesh, displacement, pressure, meshFile, 0 );

    
    for ( unsigned n = 0; n < numSteps; n++ ) {

        const double time = (n+1) * deltaT;
        std::cout << time << "  ";

        // initialise current displacement to zero (linear elasticity)
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoFU::clearValue, _1 ) );

        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDoFsU + numDoFsP );

        //------------------------------------------------------------------
        base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver, 
                                                         field, hyperElastic );

        base::asmb::stiffnessMatrixComputation<TopRight>( quadrature, solver,
                                                          field, gradP );

        base::asmb::stiffnessMatrixComputation<BotRight>( quadrature, solver,
                                                          field, laplace);

        // time integration terms
        base::time::computeReactionTerms<BotLeft,MSM>( divU, quadrature, solver,
                                                       field, deltaT, n );

        base::time::computeInertiaTerms<BotRight,MSM>( quadrature, solver,
                                                       field, deltaT, n, c0 );
        
        base::time::computeResidualForceHistory<BotRight,MSM>( laplace, 
                                                               quadrature, solver,
                                                               field, n );


        // Neumann boundary condition
        base::asmb::neumannForceComputation<SFTB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                                   boost::bind( &neumannBC<dim>, _1, _2,
                                                                H, F ) );


        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.superLUSolve();
            
        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, displacement );
        base::dof::setDoFsFromSolver( solver, pressure );
        
        // push history
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoFU::pushHistory, _1 ) );
        std::for_each( pressure.doFsBegin(), pressure.doFsEnd(),
                       boost::bind( &DoFP::pushHistory, _1 ) );


        // write VTK file
        writeVTK( mesh, displacement, pressure, meshFile, n+1 );
        
        // Finished time steps
        //--------------------------------------------------------------------------

        const std::pair<double,double> approx = probe( mesh, displacement, pressure );

        const double p = terzaghi.pressure( H/2., time );
        const double u = terzaghi.displacement( H/2., time );
        
        std::cout << u  << "  " << approx.first << "  "
                  << p  << "  " << approx.second << '\n';
    }
    
    return 0;
}
