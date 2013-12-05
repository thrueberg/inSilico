// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
// mesh related
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
// input/output
#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>
// quadrature
#include <base/Quadrature.hpp>
// FE basis
#include <base/fe/Basis.hpp>
// Field and degrees of freedom
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>
// assembly
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
// system solver
#include <base/solver/Eigen3.hpp>
// material
#include <mat/hypel/StVenant.hpp>
#include <mat/hypel/CompNeoHookean.hpp>
#include <mat/Lame.hpp>
// integral kernels
#include <solid/HyperElastic.hpp>
// time integration
#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>
// post processing
#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>

// tolerance for coordinate identification
static const double coordTol = 1.e-5;

// applied traction
double tractionValue( const double time )
{
    return ( time > 0.1 ? -10000.0 : 0.0 );
}

// Fix left side boundary (x_1=0)
template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM>::Type& x,
                  DOF* doFPtr ) 
{
    typedef typename base::Vector<DIM>::Type  VecDim;
    
    // location at x_1 = 0 or x_1 = 1
    const bool onLeftBdr  = ( std::abs( x[0] -  0. ) < coordTol );

    // Fix left boundary at x_0 = 0
    if ( onLeftBdr ) {
        for ( unsigned d = 0; d < DOF::size; d++ ) {
            if ( doFPtr -> isActive(d) )
                doFPtr -> constrainValue( d, 0.0 );
        }
    }

    return;
}

// Apply a normal traction to right side boundary (x_2=0)
template<unsigned DIM>
typename base::Vector<DIM>::Type
neumannBC( const typename base::Vector<DIM>::Type& x,
           const typename base::Vector<DIM>::Type& normal,
           const double value )
{
    typedef typename base::Vector<DIM>::Type VecDim;
    
    VecDim result = VecDim::Constant( 0. );

    const bool onRightBdr = ( std::abs( x[0] -  1. ) < coordTol );
    
    if ( onRightBdr ) result = value * normal;
        
    return result;
}

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
void writeVTKFile( const std::string& baseName,
                   const unsigned     step,
                   const MESH&        mesh,
                   const FIELD&       disp,
                   const FIELD&       veloc )
{
    // create file name with step number
    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );

    base::io::vtk::writePointData( vtkWriter, mesh, disp,  "disp" );
    base::io::vtk::writePointData( vtkWriter, mesh, veloc, "veloc" );
    
    vtk.close();
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  input.dat \n";
        return 0;
    }

    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 1;
    const base::Shape shape    = base::QUAD;
    const unsigned kernelDegEstimate = 3;
    const unsigned tiOrder = 1;

    // choice of material
    typedef mat::hypel::StVenant Material;

    // choose a time stepping method
    typedef base::time::AdamsMoulton<tiOrder> MSM;

    // time stepping method determines the history size
    const unsigned nHist = MSM::numSteps;

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double E, nu, density, stepSize;
    unsigned numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "density",          density );
        prop.registerPropertiesVar( "stepSize",         stepSize );
        prop.registerPropertiesVar( "numSteps",         numSteps );
                
        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        prop.readValuesAndCheck( inp );
        inp.close( );
    }

    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
    // Create a mesh
    typedef base::Unstructured<shape,geomDeg>    Mesh;
    const unsigned dim = Mesh::Node::dim;

    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // quadrature objects for volume and surface
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // Create a field
    const unsigned    doFSize = dim;
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>     Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field displacement, velocity;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, displacement );
    base::dof::generate<FEBasis>( mesh, velocity     );

    // Creates a list of <Element,faceNo> pairs along the boundary
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // constrain the boundary
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, displacement, 
                                           boost::bind( &dirichletBC<dim,DoF>,
                                                        _1, _2 ) );
    // same for velocity
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, velocity,
                                           boost::bind( &dirichletBC<dim,DoF>,
                                                        _1, _2 ) );

    // Create a boundary mesh from this list
    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
    }

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, displacement, velocity );
    typedef FieldBinder::TupleBinder<1,1>::Type DD;
    typedef FieldBinder::TupleBinder<1,2>::Type DV;
    typedef FieldBinder::TupleBinder<2,1>::Type VD;
    typedef FieldBinder::TupleBinder<2,2>::Type VV;
    
    // Bind the field to the boundary mesh
    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, displacement );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SFTB;

    // material object
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    // matrix kernel
    typedef solid::HyperElastic<Material,DD::Tuple> HyperElastic;
    HyperElastic linearElastic( material );

    typedef base::kernel::Mass<DV::Tuple> Mass;
    Mass mass( density );
    
    typedef base::kernel::Mass<VD::Tuple> Identity1;
    Identity1 identity1( -1. );
    
    typedef base::kernel::Mass<VV::Tuple> Identity2;
    Identity2 identity2( +1. );
            
    // Number the degrees of freedom
    const std::size_t numDoFsD =
        base::dof::numberDoFsConsecutively( displacement.doFsBegin(), displacement.doFsEnd() );
    const std::size_t numDoFsV =
        base::dof::numberDoFsConsecutively( velocity.doFsBegin(), velocity.doFsEnd(), numDoFsD );
    
    std::cout << "# Number of dofs (" << numDoFsD << ", "
              << numDoFsV << ") " << std::endl;

    // write a vtk file
    writeVTKFile( baseName, 0, mesh, displacement, velocity );

    // point to check
    typedef Mesh::Node::VecDim VecDim;
    VecDim x;
    for ( unsigned d = 0; d < dim; d++ ) x[d] = 0.5;

    // find point in mesh
    std::pair<Mesh::Element*,VecDim> probe =
        base::post::findLocationInMesh( mesh, x, coordTol, 10 );

    // prepare a monitor
    base::post::Monitor<Mesh::Element,Field::Element>
        monitorD( probe.first, 
                  displacement.elementPtr( (probe.first) -> getID() ),
                  probe.second );

    base::post::Monitor<Mesh::Element,Field::Element>
        monitorV( probe.first, 
                  velocity.elementPtr( (probe.first) -> getID() ),
                  probe.second );
    
    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < numSteps; step++ ) {

        // linear elasticity !!
        base::dof::clearDoFs( displacement );

        // time
        const double time = step * stepSize;

        // applied traction
        const double tracValue = tractionValue( time );

        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDoFsD + numDoFsV );

        // value of applied traction
        base::asmb::neumannForceComputation<SFTB>( surfaceQuadrature, solver,
                                                   surfaceFieldBinder,
                                                   boost::bind( &neumannBC<dim>,
                                                                _1, _2, tracValue ) );

        // Compute element stiffness matrices and assemble them
        base::asmb::stiffnessMatrixComputation<DD>( quadrature, solver,
                                                    fieldBinder, linearElastic );

        // Compute element stiffness matrices and assemble them
        base::asmb::stiffnessMatrixComputation<VV>( quadrature,  solver,
                                                    fieldBinder, identity2 );


        // Reaction terms of time integrator
        base::time::computeReactionTerms<DV,MSM>( mass, quadrature, solver,
                                                  fieldBinder, stepSize, step );

        base::time::computeReactionTerms<VD,MSM>( identity1, quadrature, solver,
                                                  fieldBinder, stepSize, step );

        // History terms of time integrator
        base::time::computeResidualForceHistory<DD,MSM>( linearElastic,
                                                         quadrature, solver,
                                                         fieldBinder, step );

        base::time::computeResidualForceHistory<VV,MSM>( identity2, 
                                                         quadrature, solver,
                                                         fieldBinder, step );

        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.superLUSolve();
            
        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, displacement );
        base::dof::setDoFsFromSolver( solver, velocity );

        // tell user something
        std::cout << time << "  " << tracValue << "  ";
        monitorD.solution( std::cout );
        monitorV.solution( std::cout );
        std::cout << std::endl;

        // write a vtk file
        writeVTKFile( baseName, step+1, mesh, displacement, velocity );

        // push history
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

        // push history
        std::for_each( velocity.doFsBegin(), velocity.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

    }
    // Finished load steps
    //--------------------------------------------------------------------------
    
    return 0;
}
