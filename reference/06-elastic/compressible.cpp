#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>
#include <base/fe/Basis.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include <base/solver/Eigen3.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include <solid/HyperElastic.hpp>
#include <solid/Stress.hpp>

const double coordTol = 1.e-6;

#include "PulledSheet.hpp"

//------------------------------------------------------------------------------
namespace ref06{

    template<typename MESH, typename DISP, typename MATERIAL>
    void writeVTKFile( const std::string& baseName,
                       const unsigned     step,
                       const MESH&        mesh,
                       const DISP&        disp,
                       const MATERIAL&    material );

    int compressible( int argc, char * argv[] );

}


//------------------------------------------------------------------------------
// Write displacement and stresses
template<typename MESH, typename DISP, typename MATERIAL>
void ref06::writeVTKFile( const std::string& baseName,
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

    const typename base::Vector<MESH::Node::dim>::Type xi =
        base::ShapeCentroid<MESH::Element::shape>::apply();

    // Bind the fields together
    typedef base::asmb::FieldBinder<const MESH,const DISP> FieldBinder;
    FieldBinder fieldBinder( mesh, disp );
    typedef typename FieldBinder::template TupleBinder<1,1>::Type FTB;

    base::io::vtk::writeCellData<FTB>( vtkWriter, fieldBinder, 
                                        boost::bind( solid::cauchy<
                                                     typename FTB::Tuple,
                                                     MATERIAL>,
                                                     _1, material, xi ), "sigma" );
    vtk.close();
}


//------------------------------------------------------------------------------
/** Compute the elastic (large!) deformation of a block of material.
 *  See problem description in ref06::PulledSheetProblem for a the boundary
 *  conditions. The problem is non-linear and therefore a Newton method is
 *  used. Moreover, the applied load (displacement or traction) is divided into
 *  load steps.
 *
 *  New features:
 *  - vectorial problem: DoFs are not scalar anymore
 *  - hyper-elasticity
 *
 */
int ref06::compressible( int argc, char * argv[] )
{
    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::HyperCubeShape<SPACEDIM>::value;

    // typedef mat::hypel::StVenant Material;
    typedef mat::hypel::NeoHookeanCompressible Material;

    // usage message
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0] << "  mesh.smf input.dat \n";
        return 0;
    }

    // read name of input file
    const std::string meshFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double E, nu, pull, traction, tolerance;
    unsigned maxIter, loadSteps;
    bool dispControlled;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
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
    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
    }

    // initial pull = (total amount) / (number of steps)
    const double firstPull = pull / static_cast<double>( loadSteps );

    // constrain the boundary
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, field, 
                                           boost::bind(
                                               &ref06::PulledSheet<dim>::dirichletBC<DoF>,
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
    base::io::Table<4>::WidthArray widths = {{ 2, 5, 5, 15 }};
    base::io::Table<4> table( widths );
    table % "Step" % "Iter" % "|F|"  % "|x|";
    std::cout << "#" << table;

    // write a vtk file
    ref06::writeVTKFile( baseName, 0, mesh, field, material );


    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < loadSteps; step++ ) {

        // rescale constraints in every load step: (newValue / oldValue)
        const double  pullFactor =
            (step == 0 ?
             static_cast<double>( step+1 ) :
             static_cast<double>( step+1 )/ static_cast<double>(step) );

        // scale constraints
        base::dof::scaleConstraints( field, pullFactor );

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
                const double tractionFactor =
                    traction * 
                    static_cast<double>(step+1) / static_cast<double>( loadSteps );

                // apply traction load
                base::asmb::neumannForceComputation<SFTB>(
                    surfaceQuadrature, solver,
                    surfaceFieldBinder,
                    boost::bind( &ref06::PulledSheet<dim>::neumannBC,
                                 _1, _2, tractionFactor ) );
            }

            // residual forces
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
            //solver.choleskySolve();
            solver.cgSolve();
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, field );

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
        ref06::writeVTKFile( baseName, step+1, mesh, field, material );
        
    }
    // Finished load steps
    //--------------------------------------------------------------------------
    
    return 0;
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    return ref06::compressible( argc, argv );
}
