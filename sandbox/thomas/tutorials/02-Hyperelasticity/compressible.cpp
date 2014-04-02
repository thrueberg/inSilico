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
// FE Basis
#include <base/fe/Basis.hpp>
// Field and degrees of freedom
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>
// assembly
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/FieldBinder.hpp>
// linear system solver
#include <base/solver/Eigen3.hpp>
// material
#include <mat/hypel/StVenant.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>
// integral kernel and stress
#include <solid/HyperElastic.hpp>
#include <solid/Stress.hpp>

// tolerance for coordinate identification
static const double coordTol = 1.e-5;

//------------------------------------------------------------------------------
// Fix x_0=0 and pull at x_1=1
template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM>::Type& x,
                  DOF* doFPtr,
                  const double value ) 
{
    // location at x_1 = 0 or x_1 = 1
    const bool onLeftBdr  = ( std::abs( x[0] -  0. ) < coordTol );
    const bool onRightBdr = ( std::abs( x[0] -  1. ) < coordTol );

    // Fix left boundary at x_0 = 0
    if ( onLeftBdr ) {
        for ( unsigned d = 0; d < DOF::size; d++ ) {
            if ( doFPtr -> isActive(d) )
                doFPtr -> constrainValue( d, 0.0 );
        }
    }

    // If assked for, apply normal displacement at x_1=1
    if (  onRightBdr ) { 
        if ( doFPtr -> isActive(0) ) {
            doFPtr -> constrainValue( 0, 0.0 );
            doFPtr -> constrainValue( 1, value );
        }
    }
        
    return;
}

//------------------------------------------------------------------------------
template<typename MESH, typename DISP, typename MATERIAL>
void writeVTKFile( const std::string& baseName,
                   const unsigned     step,
                   const MESH&        mesh,
                   const DISP&        disp,
                   const MATERIAL&    material )
{
    // create file name with step number
    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );

    base::io::vtk::writePointData( vtkWriter, mesh, disp, "disp" );
    base::io::vtk::writeCellData( vtkWriter, mesh, disp, 
                                  boost::bind( solid::cauchy<typename MESH::Element,
                                                             typename DISP::Element,
                                                             MATERIAL>,
                                               _1, _2, material ), "sigma" );
            
    vtk.close();
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // spatial dimension
    const unsigned    dim = 2;
    
    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 1;
    const base::Shape shape    = base::HyperCubeShape<dim>::value;
    const unsigned kernelDegEstimate = 5;

    // Type of material
    //typedef mat::hypel::StVenant Material;
    typedef mat::hypel::NeoHookeanCompressible Material;

    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  input.dat \n";
        return 0;
    }

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double E, nu, pull, tolerance;
    unsigned maxIter, loadSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "pull",             pull );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "loadSteps",        loadSteps );
        prop.registerPropertiesVar( "tolerance",        tolerance );

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

    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // quadrature objects for volume and surface
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // Create a field
    const unsigned    doFSize = dim;
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize>           Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field displacement;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, displacement );

    // Creates a list of <Element,faceNo> pairs along the boundary
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // constrain the boundary
    const double firstPull = pull / static_cast<double>( loadSteps );
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, displacement, 
                                           boost::bind( &dirichletBC<dim,DoF>,
                                                        _1, _2, firstPull ) );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, displacement );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // material object
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    // matrix kernel
    typedef solid::HyperElastic<Material,FTB::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );
            
    // Number the degrees of freedom
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( displacement.doFsBegin(), displacement.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs << std::endl;

    // create table for writing the convergence behaviour of the nonlinear solves
    base::io::Table<4>::WidthArray widths = {{ 2, 10, 10, 10 }};
    base::io::Table<4> table( widths );
    table % "Load step" % "iteration" % "|F|"  % "|x|";
    std::cout << "#" << table;

    // write a vtk file
    writeVTKFile( baseName, 0, mesh, displacement, material );

    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < loadSteps; step++ ) {

        const double pullFactor =
            (step == 0 ? 1. :
             static_cast<double>(step+1.)/static_cast<double>(step) );

        base::dof::scaleConstraints( displacement, pullFactor );

        //----------------------------------------------------------------------
        // Nonlinear iterations
        //----------------------------------------------------------------------
        unsigned iter = 0;
        while ( iter < maxIter ) {

            table % step % iter;
    
            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( numDofs );

            // Residual forces
            base::asmb::computeResidualForces<FTB>( quadrature, solver,
                                                    fieldBinder,
                                                    hyperElastic );
            
            // Compute element stiffness matrices and assemble them
            base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                         fieldBinder,
                                                         hyperElastic );

            // Finalise assembly
            solver.finishAssembly();

            // norm of residual 
            const double conv1 = solver.norm();
            table % conv1;

            // convergence via residual norm
            if ( conv1 < tolerance * E ) { // note the tolerance multiplier
                std::cout << table;
                break;
            }

            // Solve
            solver.choleskySolve();
            //solver.cgSolve();
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, displacement );

            // norm of displacement increment
            const double conv2 = solver.norm();
            table % conv2;
            std::cout << table;
            iter++;
            
            // convergence via increment
            if ( conv2 < tolerance ) break;
        }
        // Finished non-linear iterations
        //----------------------------------------------------------------------

        // warning
        if ( iter == maxIter ) {
            std::cout << "# (WW) Step " << step << " has not converged within "
                      << maxIter << " iterations \n";
        }

        // write a vtk file
        writeVTKFile( baseName, step+1, mesh, displacement, material );
        
    }
    // Finished load steps
    //--------------------------------------------------------------------------
    
    return 0;
}
