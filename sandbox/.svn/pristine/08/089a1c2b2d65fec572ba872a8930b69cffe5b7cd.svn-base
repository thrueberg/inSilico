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
#include <base/io/Format.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/Format.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
// quadrature
#include <base/Quadrature.hpp>
// FE basis, field and degrees of freedom
#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
// assembly
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
// fluid terms
#include <fluid/Stokes.hpp>
#include <fluid/Convection.hpp>
// food terms
#include <mat/thermal/IsotropicConstant.hpp>
#include <heat/Static.hpp>
#include <heat/Convection.hpp>
// system solver
#include <base/solver/Eigen3.hpp>
// time integration
#include <base/time/BDF.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

static const double coordTol = 1.e-6;

//------------------------------------------------------------------------------
// Driven cavity boundary conditions
template<unsigned DIM, typename DOF>
void dirichletBCVelocity( const typename base::Vector<DIM>::Type& x,
                          DOF* doFPtr, const double value ) 
{
    // if d-th coordinate has the value 1.0
    bool onLid = ( std::abs( x[DIM-1] - 1.0 ) < coordTol );
    
    // remove the corner/edge locations
    for ( unsigned d = 0; d < DIM-1; d ++ ) {
        if ( std::abs( x[d] - 0.0 ) < coordTol ) onLid = false;
        if ( std::abs( x[d] - 1.0 ) < coordTol ) onLid = false;
    }

    // boundary condition is either 0 or the e_1 vector 
    for ( unsigned d = 0; d < DIM; d ++ ) {
        const double tmp = ( (d==0) and onLid ) ? value : 0.0;
        
        doFPtr -> constrainValue( d, tmp );
        
    }
}

//------------------------------------------------------------------------------
// Provide some food concentration at the top
template<unsigned DIM, typename DOF>
void dirichletBCFood( const typename base::Vector<DIM>::Type& x,
                      DOF* doFPtr, const double value ) 
{
    // if d-th coordinate has the value 1.0
    bool onLid = ( std::abs( x[DIM-1] - 1.0 ) < coordTol );
    
    // remove the corner/edge locations
    for ( unsigned d = 0; d < DIM-1; d ++ ) {
        if ( std::abs( x[d] - 0.0 ) < coordTol ) onLid = false;
        if ( std::abs( x[d] - 1.0 ) < coordTol ) onLid = false;
    }

    if ( onLid ) doFPtr -> constrainValue( 0, value );
    else         doFPtr -> constrainValue( 0, 0.0 );

}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input.dat \n\n";
        return -1;
    }

    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
    // degrees of lowest-order TH element
    const unsigned    fieldDegU = 2; 
    const unsigned    fieldDegP = 1;
    const unsigned    fieldDegF = 1; 
    const base::Shape shape     = base::QUAD;
    const unsigned    tiOrder   = 1;

    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    //--------------------------------------------------------------------------
    std::string meshFile;
    double viscosity, density, tolerance, stepSize, kappa,
        foodDensity, foodBoundary;
    unsigned maxIter, numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile     );
        prop.registerPropertiesVar( "viscosity",        viscosity    );
        prop.registerPropertiesVar( "density",          density      );
        prop.registerPropertiesVar( "maxIter",          maxIter      );
        prop.registerPropertiesVar( "tolerance",        tolerance    );
        prop.registerPropertiesVar( "stepSize",         stepSize     );
        prop.registerPropertiesVar( "numSteps",         numSteps     );
        prop.registerPropertiesVar( "kappa",            kappa        );
        prop.registerPropertiesVar( "foodDensity",      foodDensity  );
        prop.registerPropertiesVar( "foodBoundary",     foodBoundary );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input erro" );
        inp.close( );
    }

    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
    // create a mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned dim = Mesh::Node::dim;

    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature 
    const unsigned kernelDegEstimate = 4;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // time integration
    typedef base::time::BDF<tiOrder> MSM;
    const unsigned nHist = MSM::numSteps;

    //--------------------------------------------------------------------------
    // Finite element bases

    // Velocity
    const unsigned    doFSizeU = dim;
    typedef base::fe::Basis<shape,fieldDegU>         FEBasisU;
    typedef base::Field<FEBasisU,doFSizeU,nHist>     Velocity;
    typedef Velocity::DegreeOfFreedom                DoFU;
    Velocity velocity;
    base::dof::generate<FEBasisU>( mesh, velocity );

    // Pressure 
    const unsigned    doFSizeP = 1;
    typedef base::fe::Basis<shape,fieldDegP>         FEBasisP;
    typedef base::Field<FEBasisP,doFSizeP,nHist>     Pressure;
    typedef Pressure::DegreeOfFreedom                DoFP;
    Pressure pressure;
    base::dof::generate<FEBasisP>( mesh, pressure );

    // Food
    const unsigned    doFSizeF = 1;
    typedef base::fe::Basis<shape,fieldDegF>         FEBasisF;
    typedef base::Field<FEBasisF,doFSizeF,nHist>     Food;
    typedef Food::DegreeOfFreedom                    DoFF;
    Food food;
    base::dof::generate<FEBasisF>( mesh, food );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Constrain the boundary 
    base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                            meshBoundary.end(),
                                            mesh, velocity,
                                            boost::bind( &dirichletBCVelocity<dim,DoFU>,
                                                         _1, _2, 1.0 ) );
    // Fix one pressure dof
    Pressure::DoFPtrIter pIter = pressure.doFsBegin();
    std::advance( pIter, std::distance( pressure.doFsBegin(), pressure.doFsEnd() )/2 );
    (*pIter) -> constrainValue( 0, 0.0 );

    // BC for the food
    base::dof::constrainBoundary<FEBasisF>( meshBoundary.begin(),
                                            meshBoundary.end(),
                                            mesh, food,
                                            boost::bind( &dirichletBCFood<dim,DoFF>,
                                                         _1, _2, foodBoundary ) );

    // Number of DoFs after constraint application!
    const std::size_t numDoFsU =
        base::dof::numberDoFsConsecutively( velocity.doFsBegin(), velocity.doFsEnd() );
    std::cout << " Number of velocity dofs " << numDoFsU << std::endl;

    const std::size_t numDoFsP =
        base::dof::numberDoFsConsecutively( pressure.doFsBegin(), pressure.doFsEnd(),
            numDoFsU );
    std::cout << " Number of pressure dofs " << numDoFsP << std::endl;

    const std::size_t numDoFsF =
        base::dof::numberDoFsConsecutively( food.doFsBegin(), food.doFsEnd() );
    std::cout << " Number of food dofs " << numDoFsF << std::endl;

    // the composite field with geometry, velocity and pressure
    typedef base::asmb::FieldBinder<Mesh,Velocity,Pressure,Food> Field;
    Field field( mesh, velocity, pressure, food );

    // define the system blocks (U,U), (U,P), and (P,U)
    typedef Field::TupleBinder<1,1,1>::Type UUU;
    typedef Field::TupleBinder<1,2>::Type   UP;
    typedef Field::TupleBinder<2,1>::Type   PU;
    typedef Field::TupleBinder<3,3,1>::Type FFU;

    // types of integral kernels
    typedef fluid::VectorLaplace<     UUU::Tuple> VecLaplace;
    typedef fluid::Convection<        UUU::Tuple> Convection;
    typedef fluid::PressureGradient<  UP::Tuple>  GradP;
    typedef fluid::VelocityDivergence<PU::Tuple>  DivU;

    typedef mat::thermal::IsotropicConstant       Material;
    typedef heat::Static<Material,FFU::Tuple>     Diffusion;
    typedef heat::Convection<FFU::Tuple>          FoodConvection;

    VecLaplace vecLaplace( viscosity );
    Convection convection( density );
    GradP      gradP;
    DivU       divU;

    Material   material( kappa );
    Diffusion  diffusion(          material );
    FoodConvection foodConvection( foodDensity );

    // system solver
    typedef base::solver::Eigen3           Solver;

    const bool incremental = false;
    
    //--------------------------------------------------------------------------
    // Time loop
    for ( unsigned step = 0; step < numSteps; step++ ) {

        // current time
        const double time = step * stepSize;

        //
        std::cout << "Time step " << step << " at time " << time
                  << std::endl;

        //base::dof::clearDoFs( velocity );
        //base::dof::clearDoFs( pressure );

#if 1
        // ramp function
        const double appliedVelocity = 
            ( time < 0.5 ? 0.0 : ( time < 1.0 ? time - 0.5 : 1.0 ) );

        // apply new velocity on the boundary
        base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                                meshBoundary.end(),
                                                mesh, velocity,
                                                boost::bind( &dirichletBCVelocity<dim,DoFU>,
                                                             _1, _2, appliedVelocity ) );
#endif        
        double prevResNorm;
        double prevSolNorm;
        
        //----------------------------------------------------------------------
        // Nonlinear fixpoint iterations
        unsigned iter = 0;
        while( iter < maxIter ) {

            // Create a solver object
            Solver solverFluid( numDoFsU + numDoFsP );

            std::cout << " - Iteration " << iter << std::flush;
    
            // Compute system matrix
            base::asmb::stiffnessMatrixComputation<UUU>( quadrature, solverFluid,
                                                         field, vecLaplace,
                                                         incremental );

            base::asmb::stiffnessMatrixComputation<UUU>( quadrature, solverFluid,
                                                         field, convection,
                                                         incremental );

            base::asmb::stiffnessMatrixComputation<UP>( quadrature, solverFluid,
                                                        field, gradP,
                                                        incremental );

            base::asmb::stiffnessMatrixComputation<PU>( quadrature, solverFluid,
                                                        field, divU,
                                                        incremental );

            if ( incremental ) {
                base::asmb::computeResidualForces<UUU>( quadrature, solverFluid,
                                                        field, vecLaplace );
                base::asmb::computeResidualForces<UUU>( quadrature, solverFluid,
                                                        field, convection );
                base::asmb::computeResidualForces<UP >( quadrature, solverFluid,
                                                        field, gradP );
                base::asmb::computeResidualForces<PU >( quadrature, solverFluid,
                                                        field, divU );
            }
            
            // compute inertia terms, d/dt, due to time integration
            base::time::computeInertiaTerms<UUU,MSM>( quadrature, solverFluid,
                                                      field, stepSize, step,
                                                      density, incremental );

            // Finalise assembly
            solverFluid.finishAssembly();

            // check convergence via solver norms
            const double residualNorm = solverFluid.norm( );
            std::cout << " |R| = " << residualNorm << std::flush;

            bool isConverged = false;
            
            if ( ( incremental ) or ( iter== 0 ) ) {
                if ( residualNorm < tolerance * viscosity ) isConverged = true;
            }
            else {
                if ( std::abs( (residualNorm - prevResNorm)/prevResNorm ) < tolerance )
                    isConverged = true;
            }
            
            if ( isConverged ) { std::cout << std::endl; break; }
            prevResNorm = residualNorm;

            // Solve
            solverFluid.superLUSolve();

            // distribute results back to dofs
            if ( incremental ) {
                base::dof::addToDoFsFromSolver( solverFluid, velocity );
                base::dof::addToDoFsFromSolver( solverFluid, pressure );
            }
            else {
                base::dof::setDoFsFromSolver( solverFluid, velocity );
                base::dof::setDoFsFromSolver( solverFluid, pressure );
            }

            // check convergence via solver norms
            const double solNorm = solverFluid.norm( );
            std::cout << " |dU| = " << solNorm << std::endl;

            if ( ( incremental ) or ( iter== 0 ) ) {
                if ( solNorm < tolerance * viscosity ) isConverged = true;
            }
            else {
                if ( std::abs( (solNorm - prevSolNorm)/prevSolNorm ) < tolerance )
                    isConverged = true;
            }
            
            if ( isConverged ) break; 
            prevSolNorm = solNorm;

            iter++;

        }

        // solve the food equation
        {
            std::cout << " * Solve food equation " << std::endl;

            base::dof::clearDoFs( food );
            
            // Create a solver object
            Solver solverFood( numDoFsF );

            // Compute system matrix from Laplacian
            base::asmb::stiffnessMatrixComputation<FFU>( quadrature, solverFood,
                                                         field, diffusion );
            
            // Compute system matrix from Convection
            base::asmb::stiffnessMatrixComputation<FFU>( quadrature, solverFood, 
                                                         field, foodConvection );

            // compute inertia terms, d/dt, due to time integration
            base::time::computeInertiaTerms<FFU,MSM>( quadrature, solverFood,
                                                      field, stepSize, step,
                                                      foodDensity );

            // Finalise assembly
            solverFood.finishAssembly();

            // Solve a possibly non-symmetric system
            solverFood.superLUSolve();

            // distribute results back to dofs
            base::dof::setDoFsFromSolver( solverFood, food );
            
            // push history
        }

        // Move in history storage (n+1 -> n)
        {
            base::dof::pushHistory( velocity );
            base::dof::pushHistory( pressure );
            base::dof::pushHistory( food     );
        }

        // output to a VTK file
        {
            // VTK Legacy
            const std::string vtkFile = baseName
                + "." + base::io::leadingZeros( step+1 ) + ".vtk";
            std::ofstream vtk( vtkFile.c_str() );
            base::io::vtk::LegacyWriter vtkWriter( vtk );
            vtkWriter.writeUnstructuredGrid( mesh );
            base::io::vtk::writePointData( vtkWriter, mesh, velocity, "U" );
            base::io::vtk::writePointData( vtkWriter, mesh, pressure, "P" );
            base::io::vtk::writePointData( vtkWriter, mesh, food,     "F" );
            vtk.close();
        }

    } // Finished time loop

    return 0;
}
