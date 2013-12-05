#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/CreateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>
#include <base/post/evaluateAtNodes.hpp>

#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/scaleConstraints.hpp>
#include <base/dof/generate.hpp>
#include <base/fe/Basis.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include <base/solver/Eigen3.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/hypel/CompNeoHookean.hpp>
#include <mat/Lame.hpp>

#include <solid/HyperElastic.hpp>
#include <solid/Stress.hpp>


//------------------------------------------------------------------------------
//  Bock of material, occupying (0,1)^DIM, fix x_1 = 0 and pull at x_1 = 1.
//  Optionally, at x_1 = 1, a surface traction is applied or a normal
//  displacement.
template<unsigned DIM>
class PulledSheetProblem
{
public:
    typedef typename base::Vector<DIM>::Type VecDim;

    // Fix x_0=0 and optionally pull at x_1=1
    template<typename DOF>
    static void dirichletBC( const VecDim& x, DOF* doFPtr,
                             const bool pullRightSide,
                             const double value ) 
    {
        // tolerance for coordinate identification
        const double tol = 1.e-5;

        // location at x_1 = 0 or x_1 = 1
        const bool onLeftBdr = ( std::abs( x[0] -  0. ) < tol );
        const bool onRightBdr = ( std::abs( x[0] -  1. ) < tol );

        // Fix left boundary at x_0 = 0
        if ( onLeftBdr ) {
            for ( unsigned d = 0; d < DOF::size; d++ ) {
                if ( doFPtr -> isActive(d) )
                    doFPtr -> constrainValue( d, 0.0 );
            }
        }

        // If assked for, apply normal displacement at x_1=1
        if (  onRightBdr and pullRightSide ) { 
            if ( doFPtr -> isActive(0) )
                doFPtr -> constrainValue( 0, value );
            
        }
        
        return;
    }

    // apply surface traction at x_1 = 1 in normal direction only
    static VecDim neumannBC( const VecDim& x,
                             const VecDim& normal,
                             const double value )
    {
        VecDim result = VecDim::Constant( 0. );

        const double tol = 1.e-5;
        
        const bool onTractionBdr =
            ( std::abs( x[0] -  1. ) < tol );

        if ( onTractionBdr ) result[1] = value;
        
        return result;
    }
    
};

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
    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;

    // typedef mat::hypel::StVenant Material;
    typedef mat::hypel::CompNeoHookean Material;

    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  input.dat \n";
        return 0;
    }

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double E, nu, pull, traction, tolerance;
    unsigned maxIter, loadSteps;
    bool dispControlled;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "pull",             pull );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "loadSteps",        loadSteps );
        prop.registerPropertiesVar( "traction",         traction );
        prop.registerPropertiesVar( "dispControlled",   dispControlled );
        prop.registerPropertiesVar( "tolerance",        tolerance );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        prop.readValues( inp );
        inp.close( );

        // Make sure all variables have been found
        if ( not prop.isEverythingRead() ) {
            prop.writeUnread( std::cerr );
            VERIFY_MSG( false, "Could not find above variables" );
        }
    }

    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

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

    // Create a field
    const unsigned    doFSize = dim;
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize>           Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field field;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, field );

    // Creates a list of <Element,faceNo> pairs along the boundary
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Create a boundary mesh from this list
    typedef base::mesh::CreateBoundaryMesh<Mesh::Element> CreateBoundaryMesh;
    typedef CreateBoundaryMesh::BoundaryMesh BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        CreateBoundaryMesh::apply( meshBoundary.boundaryBegin(),
                                   meshBoundary.boundaryEnd(),
                                   mesh, boundaryMesh );
    }

    // constrain the boundary
    const double firstPull = pull / static_cast<double>( loadSteps );
    base::dof::constrainBoundary<FEBasis>( meshBoundary.boundaryBegin(),
                                           meshBoundary.boundaryEnd(),
                                           mesh, field, 
                                           boost::bind( &PulledSheetProblem<dim>::dirichletBC<DoF>,
                                                        _1, _2, dispControlled, firstPull ) );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, field );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SFTB;

    // material object
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    // matrix kernel
    typedef solid::HyperElastic<Material,FTB::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );
            
    // Number the degrees of freedom
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs << std::endl;

    // create table for writing the convergence behaviour of the nonlinear solves
    base::io::Table<4>::WidthArray widths = {{ 2, 10, 10, 10 }};
    base::io::Table<4> table( widths );
    table % "Load step" % "iteration" % "|F|"  % "|x|";
    std::cout << "#" << table;

    // write a vtk file
    writeVTKFile( baseName, 0, mesh, field, material );


    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < loadSteps; step++ ) {

        //----------------------------------------------------------------------
        // Nonlinear iterations
        //----------------------------------------------------------------------
        unsigned iter = 0;
        while ( iter < maxIter ) {

            table % step % iter;
    
            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( numDofs );

            // apply traction boundary condition, if problem is not disp controlled
            if ( not dispControlled ) {
                // value of applied traction
                const double tracValue =
                    static_cast<double>(step+1) / static_cast<double>( loadSteps );
                base::asmb::neumannForceComputation<SFTB>( surfaceQuadrature, solver,
                                                           surfaceFieldBinder,
                                                           boost::bind( &PulledSheetProblem<dim>::
                                                                        neumannBC,
                                                                        _1, _2, tracValue ) );
            }
            
            base::asmb::computeResidualForces<FTB>( quadrature, solver,
                                                    fieldBinder,
                                                    hyperElastic );
            
            // Compute element stiffness matrices and assemble them
            base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                         fieldBinder,
                                                         hyperElastic,
                                                         iter > 0 );

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
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, field, iter > 0 );

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
        writeVTKFile( baseName, step+1, mesh, field, material );
        
    }
    // Finished load steps
    //--------------------------------------------------------------------------
    
    return 0;
}
