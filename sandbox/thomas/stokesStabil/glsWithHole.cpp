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
#include <base/cut/ScaledField.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/location.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>

#include <base/io/Format.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/bruteForce.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/stabiliseBasis.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>


#include <fluid/Stokes.hpp>
#include <fluid/Convection.hpp>
#include <fluid/GalerkinLeastSquares.hpp>

#include <base/solver/Eigen3.hpp>


#define SDIV
#define HOELLIG
#define STABIL

const double coordTol = 1.e-5;


//------------------------------------------------------------------------------
//! Analytic level set function for a spherical domain
template<unsigned DIM>
bool spherical( const typename base::Vector<DIM>::Type& x,
                typename base::Vector<DIM>::Type& xClosest,
                const typename base::Vector<DIM>::Type& c,
                const double R )
{
    const typename base::Vector<DIM>::Type y = x - c;
    
    if ( y.norm() < coordTol ) {
        xClosest =    c;
        xClosest[0] = R;
    }
    else{
        xClosest = (R / y.norm()) * y + c;
    }

    if ( y.norm() > R ) return true;
    return false;
}

//------------------------------------------------------------------------------
//! Analytic level set function for a spherical domain
template<unsigned DIM>
bool vertical( const typename base::Vector<DIM>::Type& x,
               typename base::Vector<DIM>::Type& xClosest,
               const double x1 )
{
    xClosest    = x;
    xClosest[0] = x1;

    return x[0] < x1;
}


//------------------------------------------------------------------------------
// Driven cavity
template<unsigned DIM, typename DOF>
void drivenCavity( const typename base::Vector<DIM>::Type& x,
                   DOF* doFPtr ) 
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
    const bool left   = (std::abs( x[0] - 0. ) < coordTol );
    const bool topBot =
        ( std::abs( x[DIM-1] - 1. ) < coordTol ) or
        ( std::abs( x[DIM-1] - 0. ) < coordTol );
        

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
    const bool left   = (std::abs( x[0] - 0. ) < coordTol );
    const bool topBot =
        ( std::abs( x[DIM-1] - 1. ) < coordTol ) or
        ( std::abs( x[DIM-1] - 0. ) < coordTol );
        

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
        std::cout << "Usage:  " << argv[0]
                  << " alpha viscosity dw meshFile.smf \n";
        return -1;
    }

    const double alpha         = boost::lexical_cast<double>( argv[1] );
    const double viscosity     = boost::lexical_cast<double>( argv[2] );
    const bool   douglasWang   = boost::lexical_cast<bool>(   argv[3] );
    const std::string meshFile = boost::lexical_cast<std::string>( argv[4] );
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned dim = Mesh::Node::dim;
    typedef Mesh::Node::VecDim VecDim;

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
    //typedef base::Field<FEBasisU,doFSizeU>           Velocity;
    typedef base::cut::ScaledField<FEBasisU,doFSizeU> Velocity;
    typedef Velocity::DegreeOfFreedom                DoFU;
    Velocity velocity;
    base::dof::generate<FEBasisU>( mesh, velocity );
    
    const unsigned    doFSizeP = 1;
    typedef base::fe::Basis<shape,fieldDegP>         FEBasisP;
    //typedef base::Field<FEBasisP,doFSizeP>           Pressure;
    typedef base::cut::ScaledField<FEBasisP,doFSizeP> Pressure;
    Pressure pressure;
    base::dof::generate<FEBasisP>( mesh, pressure );

    // find geometry association for the dofs
    std::vector<std::pair<std::size_t,VecDim> > doFLocationP, doFLocationU;
    base::dof::associateLocation( pressure, doFLocationP );
    base::dof::associateLocation( velocity, doFLocationU );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    const VecDim c = base::constantVector<dim>( 0.5 );
    const double R = 0.3;
    base::cut::analyticLevelSet( mesh,
                                 boost::bind( &spherical<dim>, _1, _2, c, R ),
                                 //boost::bind( &vertical<dim>, _1, _2, 0.77777777 ),
                                 true, levelSet );

    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells;
    base::cut::generateCutCells( mesh, levelSet, cells );

    typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
    CutQuadrature cutQuadrature( cells );


// Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                            meshBoundary.end(),
                                            mesh, velocity,
                                            boost::bind( &poiseuille<dim,DoFU>, _1, _2 ) );
                                            //boost::bind( &shear<dim,DoFU>, _1, _2 ) );

    // Fix one pressure dof
    //Pressure::DoFPtrIter pIter = pressure.doFsBegin();
    //std::advance( pIter, std::distance( pressure.doFsBegin(), pressure.doFsEnd() )/2 );
    //(*pIter) -> constrainValue( 0, 0.0 );

    // compute supports
    std::vector<double> supportsU, supportsP;
    base::cut::supportComputation( mesh, velocity, cutQuadrature, supportsU );
    base::cut::supportComputation( mesh, pressure, cutQuadrature, supportsP );

#ifdef HOELLIG
    base::cut::stabiliseBasis( mesh, velocity, supportsU, doFLocationU );
    base::cut::stabiliseBasis( mesh, pressure, supportsP, doFLocationP );
#else
    velocity.tagBasis( supportsU, 1.e-10 );
    pressure.tagBasis( supportsP, 1.e-10 );
#endif

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
    base::asmb::stiffnessMatrixComputation<TopLeft>( cutQuadrature, solver,
                                                     field, stressDiv );
#else
    base::asmb::stiffnessMatrixComputation<TopLeft>( cutQuadrature, solver,
                                                     field, vecLaplace );
#endif

    base::asmb::stiffnessMatrixComputation<TopRight>( cutQuadrature, solver,
                                                      field, gradP);

    base::asmb::stiffnessMatrixComputation<BotLeft>( cutQuadrature, solver,
                                                     field, divU );

    // Stabilisation terms
#ifdef STABIL
    
#ifdef SDIV
    base::asmb::stiffnessMatrixComputation<TopLeft>(cutQuadrature, solver,
                                                    field, stressDivStabil );
    base::asmb::stiffnessMatrixComputation<TopRight>( cutQuadrature, solver,
                                                      field, gradPStabil2 );
    base::asmb::stiffnessMatrixComputation<BotLeft>( cutQuadrature, solver,
                                                     field, divUStabil2 );
#else
    base::asmb::stiffnessMatrixComputation<TopLeft>( cutQuadrature, solver,
                                                     field, vecLaplaceStabil );
    base::asmb::stiffnessMatrixComputation<TopRight>( cutQuadrature, solver,
                                                      field, gradPStabil  );
    base::asmb::stiffnessMatrixComputation<BotLeft>( cutQuadrature, solver,
                                                     field, divUStabil  );
#endif

    base::asmb::stiffnessMatrixComputation<BotRight>( cutQuadrature, solver,
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
        std::vector<double> distances;
        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( distances ),
                        boost::bind( &LevelSet::getSignedDistance, _1 ) );
        vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
        vtk.close();
    }

    return 0;
}
