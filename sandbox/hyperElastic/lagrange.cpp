// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
// mesh related
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
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
#include <base/dof/generate.hpp>
#include <base/dof/constrainBoundary.hpp>
// assembly
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>
// linear system solver
#include <base/solver/Eigen3.hpp>
// material
#include <mat/hypel/StVenant.hpp>
#include "NeoHookeanCompressible.hpp"
#include <mat/Lame.hpp>
// integral kernel and stress
#include <solid/Stress.hpp>

#include <base/post/ErrorNorm.hpp>

#include "Lagrange.hpp"
#include <solid/HyperElastic.hpp>

static const double coordTol = 1.e-8;

//------------------------------------------------------------------------------
// Function for the point-wise constraint of the Boundary
//[dirichlet]{
template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM>::Type& x, DOF* doFPtr )
{
    for ( unsigned d = 0; d < DIM; d++ ) {
        const bool clamp ( (std::abs( x[d] - 0. ) < coordTol ) or
                           (std::abs( x[d] - 1. ) < coordTol ) );

        if ( clamp ) {
            if ( doFPtr -> isActive(d) ) {
                doFPtr -> constrainValue( d, 0.0 );
            }
        }
    }
}

//------------------------------------------------------------------------------
template<typename MESH, typename DISP>
void writeVTKFile( const std::string& baseName,
                   const unsigned     step,
                   const MESH&        mesh,
                   const DISP&        disp )
{
    // create file name with step number
    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );

    base::io::vtk::writePointData( vtkWriter, mesh, disp, "disp" );

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
    const unsigned kernelDegEstimate = 3;

    // usage message
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0] << "  mesh.smf input.dat \n"
                  << "Compiled for dim=" << dim << "\n\n";
        return 0;
    }

    // read name of input file
    const std::string meshFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string inputFile = boost::lexical_cast<std::string>( argv[2] );


    // read from input file
    double E, nu, ubar1, ubar2, ubar3, tolerance;
    unsigned maxIter;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "ubar1",            ubar1 );
        prop.registerPropertiesVar( "ubar2",            ubar2 );
        prop.registerPropertiesVar( "ubar3",            ubar3 );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "tolerance",        tolerance );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );
    }

    // find base name from mesh file
    const std::string baseName = "lagrange";

    const double lambda = mat::Lame::lambda( E, nu );
    const double mu     = mat::Lame::mu(     E, nu );

    boost::array<double,dim> ubar;
    {
        const boost::array<double,3> aux = {{ ubar1, ubar2, ubar3 }};
        for ( unsigned d = 0; d < dim; d++ ) ubar[d] = aux[d];
    }

    const AnalyticLagrangeTensor<dim> analyticLagrange( lambda, mu, ubar );

    //--------------------------------------------------------------------------
    // Create a mesh
    typedef base::Unstructured<shape,geomDeg>    Mesh;

    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
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

    // constrain the boundary
    {
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                               meshBoundary.end(),
                                               mesh, displacement,
                                               boost::bind( &dirichletBC<dim,
                                                            Field::DegreeOfFreedom>,
                                                            _1, _2 ) );
        
    }

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, displacement );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // matrix kernel
    typedef mat::hypel::NeoHookeanCompressible Material;
    Material material( lambda, mu );
    typedef solid::HyperElastic<Material,FTB::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );
            
    // Number the degrees of freedom
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( displacement.doFsBegin(), displacement.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs << std::endl;

    // write a vtk file
    writeVTKFile( baseName, 0, mesh, displacement );

    //----------------------------------------------------------------------
    // Nonlinear iterations
    //----------------------------------------------------------------------
    unsigned iter = 0;
    while ( iter < maxIter ) {

        std::cout << iter << "  ";
        
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
                                                     hyperElastic,
                                                     iter > 0 );

        // Body force
        base::asmb::bodyForceComputation<FTB>( quadrature, solver, fieldBinder,
                                               boost::bind( &AnalyticLagrangeTensor<dim>::force,
                                                            &analyticLagrange, _1 ) );

        // Finalise assembly
        solver.finishAssembly();

        // norm of residual 
        const double conv1 = solver.norm();
        std::cout << conv1 << "  ";

        // convergence via residual norm
        if ( conv1 < tolerance * E ) { // note the tolerance multiplier
            std::cout << std::endl;
            break;
        }

        // Solve
        solver.choleskySolve();
            
        // distribute results back to dofs
        base::dof::addToDoFsFromSolver( solver, displacement, iter > 0 );

        // write a vtk file
        writeVTKFile( baseName, iter+1, mesh, displacement );
        
        // norm of displacement increment
        const double conv2 = solver.norm();
        std::cout << conv2 << std::endl;
        iter++;
            
        // convergence via increment
        if ( conv2 < tolerance ) break;
    }
    // Finished non-linear iterations
    //----------------------------------------------------------------------
    
    // warning
    if ( iter == maxIter ) {
        std::cout << "# (WW) Solution has not converged within "
                  << maxIter << " iterations \n";
    }
    
    // Finished load steps
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // compute L2-error
    std::cout << "L2-error = "
              << base::post::errorComputation<0>(
                  quadrature, mesh, displacement,
                  boost::bind( &AnalyticLagrangeTensor<dim>::solution, 
                               &analyticLagrange, _1 ) )
              << '\n';

    return 0;
}
