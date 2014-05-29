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
#include <fluid/GalerkinLeastSquares.hpp>

#include <base/solver/Eigen3.hpp>



#define SDIV
#define STABIL

const double tol = 1.e-5;

//------------------------------------------------------------------------------
// Driven cavity
template<unsigned DIM, typename DOF>
void drivenCavity( const typename base::Vector<DIM>::Type& x,
                   DOF* doFPtr ) 
{
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

//------------------------------------------------------------------------------
// Simple shear
template<unsigned DIM, typename DOF>
void shear( const typename base::Vector<DIM>::Type& x,
            DOF* doFPtr ) 
{
    const bool left   = (std::abs( x[0] - 0. ) < tol );
    const bool topBot =
        ( std::abs( x[DIM-1] - 1. ) < tol ) or
        ( std::abs( x[DIM-1] - 0. ) < tol );
        

    //
    if ( left or topBot ) {
        for ( unsigned d = 0; d < DIM; d ++ ) {
            const double value = (d==0 ? x[1] : 0. );
            if ( doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, value );
        }
    }
}

//------------------------------------------------------------------------------
// Poiseuille
template<unsigned DIM, typename DOF>
void poiseuille( const typename base::Vector<DIM>::Type& x,
                 DOF* doFPtr ) 
{
    const bool left   = (std::abs( x[0] - 0. ) < tol );
    const bool topBot =
        ( std::abs( x[DIM-1] - 1. ) < tol ) or
        ( std::abs( x[DIM-1] - 0. ) < tol );
        

    //
    if ( left or topBot ) {
        for ( unsigned d = 0; d < DIM; d ++ ) {
            const double value = (d==0 ? x[1]*(1.-x[1]) : 0. );
            if ( doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, value );
        }
    }
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
#ifdef STABIL
    const unsigned    fieldDegU = 1; 
    const unsigned    fieldDegP = 1;
#else
    const unsigned    fieldDegU = 2; 
    const unsigned    fieldDegP = 1;
#endif
    const base::Shape shape     = base::QUAD;

    //--------------------------------------------------------------------------
    if ( argc != 5 ) {
        std::cout << "Usage:  " << argv[0] << " alpha viscosity dw meshFile.smf \n";
        return -1;
    }

    const double alpha         = boost::lexical_cast<double>(      argv[1] );
    const double viscosity     = boost::lexical_cast<double>(      argv[2] );
    const bool   douglasWang   = boost::lexical_cast<double>(      argv[3] );
    const std::string meshFile = boost::lexical_cast<std::string>( argv[4] );
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
                                            boost::bind( &poiseuille<dim,DoFU>, _1, _2 ) );
                                            //boost::bind( &shear<dim,DoFU>, _1, _2 ) );

    // Fix one pressure dof
    Pressure::DoFPtrIter pIter = pressure.doFsBegin();
    //std::advance( pIter, std::distance( pressure.doFsBegin(), pressure.doFsEnd() )/2 );
    //(*pIter) -> constrainValue( 0, 0.0 );

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
    typedef Field::TupleBinder<1,1>::Type   TopLeft;
    typedef Field::TupleBinder<1,2>::Type   TopRight;
    typedef Field::TupleBinder<2,1>::Type   BotLeft;
    typedef Field::TupleBinder<2,2>::Type   BotRight;
    
    typedef fluid::VectorLaplace<     TopLeft::Tuple>  VecLaplace;
    typedef fluid::StressDivergence<  TopLeft::Tuple>  StressDiv;
    typedef fluid::PressureGradient< TopRight::Tuple>  GradP;
    typedef fluid::VelocityDivergence<BotLeft::Tuple>  DivU;

    VecLaplace vecLaplace( viscosity );
    StressDiv  stressDiv( viscosity );
    GradP      gradP;
    DivU       divU;

    //--------------------------------------------------------------------------
    typedef fluid::gls::StressDivergence<   TopLeft::Tuple>  StressDivStabil;
    typedef fluid::gls::VectorLaplace<      TopLeft::Tuple>  VecLaplaceStabil;
    typedef fluid::gls::PressureGradient<  TopRight::Tuple>  GradPStabil;
    typedef fluid::gls::VelocityDivergence< BotLeft::Tuple>  DivUStabil;
    typedef fluid::gls::PressureGradient2< TopRight::Tuple>  GradPStabil2;
    typedef fluid::gls::VelocityDivergence2<BotLeft::Tuple>  DivUStabil2;
    typedef fluid::gls::PressureLaplace<   BotRight::Tuple>  PressureLaplace;

    StressDivStabil  stressDivStabil(  alpha, viscosity, douglasWang );
    VecLaplaceStabil vecLaplaceStabil( alpha, viscosity, douglasWang );
    GradPStabil      gradPStabil(      alpha, viscosity, douglasWang );
    DivUStabil       divUStabil(       alpha, viscosity );
    GradPStabil2     gradPStabil2(     alpha, viscosity, douglasWang );
    DivUStabil2      divUStabil2(      alpha, viscosity );
    PressureLaplace  pressureLaplace(  alpha );

    //--------------------------------------------------------------------------
    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDoFsU + numDoFsP );
    solver.registerFields<TopLeft>(  field );
    solver.registerFields<TopRight>( field );
    solver.registerFields<BotLeft>(  field );
    solver.registerFields<BotRight>( field );

    // Compute system matrix
#ifdef SDIV
    base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver,
                                                     field, stressDiv );
#else
    base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver,
                                                     field, vecLaplace );
#endif

    base::asmb::stiffnessMatrixComputation<TopRight>( quadrature, solver,
                                                      field, gradP);

    base::asmb::stiffnessMatrixComputation<BotLeft>( quadrature, solver,
                                                     field, divU );

    // Stabilisation terms
#ifdef STABIL
    
#ifdef SDIV
    base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver,
                                                     field, stressDivStabil );
    base::asmb::stiffnessMatrixComputation<TopRight>( quadrature, solver,
                                                      field, gradPStabil2 );
    base::asmb::stiffnessMatrixComputation<BotLeft>( quadrature, solver,
                                                     field, divUStabil2 );
#else
    base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver,
                                                     field, vecLaplaceStabil );
    base::asmb::stiffnessMatrixComputation<TopRight>( quadrature, solver,
                                                      field, gradPStabil );
    base::asmb::stiffnessMatrixComputation<BotLeft>( quadrature, solver,
                                                     field, divUStabil );
#endif

    base::asmb::stiffnessMatrixComputation<BotRight>( quadrature, solver,
                                                      field, pressureLaplace );

#endif
    
    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.superLUSolve();

    // distribute results back to dofs
    base::dof::setDoFsFromSolver( solver, velocity );
    base::dof::setDoFsFromSolver( solver, pressure );

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
