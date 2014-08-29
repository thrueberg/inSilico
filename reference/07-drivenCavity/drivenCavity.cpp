#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>

#include <base/mesh/MeshBoundary.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/Format.hpp>
#include <base/io/PropertiesParser.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>

#include <base/io/Format.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <fluid/Stokes.hpp>
#include <fluid/Convection.hpp>

#include <base/solver/Eigen3.hpp>

#define INCREMENTAL

//------------------------------------------------------------------------------
namespace ref07{

    //--------------------------------------------------------------------------
    // Function for the point-wise constraint of the Boundary
    template<unsigned DIM, typename DOF>
    void dirichletBC( const typename base::Vector<DIM>::Type& x,
                      DOF* doFPtr ) 
    {
        const double tol = 1.e-5;

        // if d-th coordinate has the value 1.0
        bool onLid = ( std::abs( x[DIM-1] - 1.0 ) < tol );
        // remove the corner/edge locations
        for ( unsigned d = 0; d < DIM-1; d ++ ) {
            if ( std::abs( x[d] - 0.0 ) < tol ) onLid = false;
            if ( std::abs( x[d] - 1.0 ) < tol ) onLid = false;
        }

        // boundary condition is either 0 or the e_1 vector 
        for ( unsigned d = 0; d < DIM; d ++ ) {
            const double value = ( (d==0) and onLid ) ? 1.0 : 0.0;
            
            if ( doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, value );
        }
    }


    int drivenCavity( int argc, char * argv[] );
}


//------------------------------------------------------------------------------
/** Solve the Driven Cavity problem of the incompressible Navier-Stokes system.
 *
 *  New features are
 *   - a mixed finite element formulation, here \f$( \vec{u}, p )\f$
 *   - use of different polynomial degrees per field
 *
 *  The Navier-Stokes system reads
 *  \f{eqnarray*}{
 *     \rho u_t + u \cdot \nabla u - \nabla \cdot \sigma(u,p) &=& f \cr
 *                                   \nabla \cdot u           &=& 0
 *  \f}
 *  where the first equation is the momentum balance of the fluid in which the
 *  total time derivative occurs. The stress is commonly due to a Newtonian
 *  material law and becomes
 *  \f[
 *        \sigma (u,p) = -p I + \mu ( \nabla u + \nabla^T u )
 *  \f]
 *  The second equation is the mass balance for an incompressible fluid.
 *  Material parameters are thus the mass density \f$ \rho \f$ and the
 *  viscosity \f$ \mu \f$. In this application, the velocity field \f$ u \f$
 *  and the pressure field \f$ p \f$ are approximated independently by
 *  Lagrangian shape functions. Moreover, the famous inf-sup stability
 *  condition does not allow as to use an equal-order approximation, that
 *  is the polynomial degress of the basis functions for \f$ u \f$ cannot
 *  be equal to the ones of \f$ p \f$. We make use of the lowest-order
 *  Taylor-Hood element (proven to be stable) and therefore use a quadratic
 *  approximation of the velocity and a linear approximation of the pressure.
 *
 *  In this application, the stationary Navier-Stokes system (\f$ u_t = 0 \f$)
 *  is solved for the case of a 'driven cavity'. This is a closed fluid-filled
 *  box which is subject to a horizontal boundary velocity on the top face.
 *  The other boundaries are zero-velocity walls. Note that this case comprises
 *  a pure Dirichlet problem and therefore the pressure is only determined up
 *  to a constant. In order to remove this null-space an arbitrarily chose
 *  pressure degree of freedom is set to zero. 
 *
 *  The figure shows the numerical solution for the Reynolds number
 *  \f$ Re = 1000 \f$. For this case, the flow field is stationary and
 *  therefore, no time integration is necessary. In the figure, some
 *  streamlines are plotted on top of the pressure contour coloring. 
 *
 *  \image html drivenCavity.png "Flow visualisation"
 *
 */
int ref07::drivenCavity( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
    // degrees of lowest-order TH element
    const unsigned    fieldDegU = 2; 
    const unsigned    fieldDegP = 1;
    const base::Shape shape     = base::QUAD;

    //--------------------------------------------------------------------------
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input.dat \n\n";
        return -1;
    }

    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    //--------------------------------------------------------------------------
    std::string meshFile;
    double viscosity, density, tolerance;
    unsigned maxIter;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "viscosity",        viscosity );
        prop.registerPropertiesVar( "density",          density );
        prop.registerPropertiesVar( "maxIter",          maxIter );
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

    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
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
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // Finite element basis
    const unsigned    doFSizeU = dim;
    typedef base::fe::Basis<shape,fieldDegU>         FEBasisU;
    typedef base::Field<FEBasisU,doFSizeU>           Velocity;
    typedef Velocity::DegreeOfFreedom                DoFU;
    Velocity velocity;
    base::dof::generate<FEBasisU>( mesh, velocity );
    
    const unsigned    doFSizeP = 1;
    typedef base::fe::Basis<shape,fieldDegP>         FEBasisP;
    typedef base::Field<FEBasisP,doFSizeP>           Pressure;
    Pressure pressure;
    base::dof::generate<FEBasisP>( mesh, pressure );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                            meshBoundary.end(),
                                            mesh, velocity,
                                            boost::bind( &dirichletBC<dim,DoFU>, _1, _2 ) );

    // Fix one pressure dof
    Pressure::DoFPtrIter pIter = pressure.doFsBegin();
    std::advance( pIter, std::distance( pressure.doFsBegin(), pressure.doFsEnd() )/2 );
    (*pIter) -> constrainValue( 0, 0.0 );

    // Number of DoFs after constraint application!
    const std::size_t numDoFsU =
        base::dof::numberDoFsConsecutively( velocity.doFsBegin(), velocity.doFsEnd() );
    std::cout << " Number of velocity dofs " << numDoFsU << std::endl;

    const std::size_t numDoFsP =
        base::dof::numberDoFsConsecutively( pressure.doFsBegin(), pressure.doFsEnd(),
            numDoFsU );
    std::cout << " Number of pressure dofs " << numDoFsP << std::endl;

    // the composite field with geometry, velocity and pressure
    typedef base::asmb::FieldBinder<Mesh,Velocity,Pressure> Field;
    Field field( mesh, velocity, pressure );

    // define the system blocks (U,U), (U,P), and (P,U)
    typedef Field::TupleBinder<1,1,1>::Type UU;
    typedef Field::TupleBinder<1,2>::Type   UP;
    typedef Field::TupleBinder<2,1>::Type   PU;
    
    typedef fluid::VectorLaplace<     UU::Tuple>  VecLaplace;
    typedef fluid::Convection<        UU::Tuple>  Convection;
    typedef fluid::PressureGradient<  UP::Tuple>  GradP;
    typedef fluid::VelocityDivergence<PU::Tuple>  DivU;

    VecLaplace vecLaplace( viscosity );
    Convection convection( density );
    GradP      gradP;
    DivU       divU;

    // for fixed-point iterations
    double prevSolNorm;
    double prevResNorm;

#ifdef INCREMENTAL
    const bool incremental = true;
#else
    const bool incremental = false;
#endif

    //--------------------------------------------------------------------------
    // Nonlinear iterations
    unsigned iter = 0;
    while( iter < maxIter ) {

        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDoFsU + numDoFsP );
        solver.registerFields<UU>(  field );
        solver.registerFields<UP>(  field );
        solver.registerFields<PU>(  field );

        std::cout << "Iteration " << iter << ": " << std::flush;
    
        // Compute system matrix
        base::asmb::stiffnessMatrixComputation<UU>( quadrature, solver,
                                                    field, vecLaplace, incremental );

        base::asmb::stiffnessMatrixComputation<UU>( quadrature, solver,
                                                    field, convection, incremental );

        base::asmb::stiffnessMatrixComputation<UP>( quadrature, solver,
                                                    field, gradP, incremental );

        base::asmb::stiffnessMatrixComputation<PU>( quadrature, solver,
                                                    field, divU, incremental );

        if ( incremental ) {
            base::asmb::computeResidualForces<UU >( quadrature, solver, field, vecLaplace );
            base::asmb::computeResidualForces<UU >( quadrature, solver, field, convection );
            base::asmb::computeResidualForces<UP>(  quadrature, solver, field, gradP );
            base::asmb::computeResidualForces<PU >( quadrature, solver, field, divU );
        }
        
        // Finalise assembly
        solver.finishAssembly();

        // check convergence via solver norms
        const double residualNorm = solver.norm( );
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
        solver.superLUSolve();
        //solver.umfPackLUSolve();

        // distribute results back to dofs
        if ( incremental ) {
            base::dof::addToDoFsFromSolver( solver, velocity );
            base::dof::addToDoFsFromSolver( solver, pressure );
        }
        else {
            base::dof::setDoFsFromSolver( solver, velocity );
            base::dof::setDoFsFromSolver( solver, pressure );
        }
    
        // check convergence via solver norms
        const double solNorm = solver.norm( );
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

    // output to a VTK file
    {
        // VTK Legacy
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        base::io::vtk::writePointData( vtkWriter, mesh, velocity, "U" );
        base::io::vtk::writePointData( vtkWriter, mesh, pressure, "P" );
        vtk.close();
    }

    return 0;
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    return ref07::drivenCavity( argc, argv );
}
