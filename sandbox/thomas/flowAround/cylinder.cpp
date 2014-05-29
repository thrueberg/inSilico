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
#include <base/io/Format.hpp>
#include <base/io/PropertiesParser.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/setField.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/SimpleIntegrator.hpp>

#include <base/io/Format.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <fluid/Stokes.hpp>
#include <fluid/Convection.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>
#include <base/nitsche/Parameters.hpp>

#include <base/cut/LevelSet.hpp>
#include <base/cut/bruteForce.hpp>
#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/ScaledField.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/TransferSurfaceDatum.hpp>
#include <base/cut/ComputeSurfaceForces.hpp>

#include <base/time/BDF.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

// tolerance for coordiantes
const double coordTol = 1.e-5;

//------------------------------------------------------------------------------
// Filter of Dirichlet BC
template<unsigned DIM>
bool boundaryFilter( const typename base::Vector<DIM>::Type& x )
{
    if ( std::abs( x[0] + 1.  ) < coordTol ) return true;
    if ( std::abs( x[1] + 0.5 ) < coordTol ) return true;
    if ( std::abs( x[1] - 1.5 ) < coordTol ) return true;

    return false;
}

//------------------------------------------------------------------------------
// Function for the point-wise constraint of the Boundary
//[dirichlet]{
template<unsigned DIM>
typename base::Vector<DIM>::Type
dirichlet( const typename base::Vector<DIM>::Type& x, const double factor )
{
    
    typename base::Vector<DIM>::Type result = base::constantVector<DIM>( 0. );

    // if 1st coordinate has the value -1
    const bool isInflow = ( std::abs( x[0] + 1.0 ) < coordTol );

    // parabolic inflow
    if ( isInflow ) {
        const double xi = x[1] - 0.5;
        const double inflowU = 1.5 * (1. - xi*xi);
        result[0] = inflowU * factor;
    }

    return result;
}
//[dirichlet]}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void surfaceVelocity( const typename base::Vector<DIM,double>::Type& x,
                      DOF* doFPtr ) 
{
    typename base::Vector<DIM,double>::Type u;
    for ( unsigned d = 0; d < DIM; d++) u[d] = 0.;
    //u[0] =  (x[1]-0.5);
    //u[1] = -(x[0]-0.5);
    
    for ( unsigned d = 0; d < DIM; d++ )
        doFPtr -> setValue( d, u[d] );
}

//------------------------------------------------------------------------------
// output to a VTK file
template<typename MESH, typename VELOCITY, typename PRESSURE, typename LSET>
void writeVTKFile( const std::string& baseName,
                   const unsigned iter,
                   MESH& mesh,
                   VELOCITY& velocity, 
                   PRESSURE& pressure,
                   const std::vector<LSET>& levelSet,
                   const double viscosity )
{
    // VTK Legacy
    const std::string vtkFile = baseName
        + "." + base::io::leadingZeros( iter ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );
    base::io::vtk::writePointData( vtkWriter, mesh, velocity, "U" );
    base::io::vtk::writePointData( vtkWriter, mesh, pressure, "P" );

    {
        std::vector<double> distances;
        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( distances ),
                        boost::bind( &LSET::getSignedDistance, _1 ) );

        vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
    }

    {
        typename base::Vector<MESH::Node::dim>::Type xi =
            base::ShapeCentroid<MESH::Element::shape>::apply();
        
        std::vector<typename base::Matrix<MESH::Node::dim,MESH::Node::dim>::Type> sigma;

        typedef base::asmb::FieldBinder<MESH,VELOCITY,PRESSURE> Field;
        Field field( mesh, velocity, pressure );

        typedef typename Field::template TupleBinder<1,2>::Type UP;
        fluid::Stress<typename UP::Tuple> stress( viscosity );

        typename Field::FieldIterator first = field.elementsBegin();
        typename Field::FieldIterator  last = field.elementsEnd();
        for ( ; first != last; ++first ) {

            sigma.push_back( stress( UP::makeTuple( *first ), xi ) );
        }
        

        vtkWriter.writeCellData( sigma.begin(), sigma.end(), "sigma" );

    }

        
    vtk.close();
}

//------------------------------------------------------------------------------
// output to a VTK file
template<typename SMESH, typename SFIELD>
void writeSurfaceVTKFile( const std::string& baseName,
                          const unsigned iter,
                          const SMESH&  mesh,
                          const SFIELD& surfVelocity,
                          const SFIELD& surfForces )
{
    // VTK Legacy
    const std::string vtkFile = baseName + ".surf."
        + base::io::leadingZeros( iter ) + ".vtk";
    
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );
    base::io::vtk::writePointData( vtkWriter, mesh, surfVelocity, "U" );
    base::io::vtk::writePointData( vtkWriter, mesh, surfForces,   "F" );
    
    vtk.close();
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
    const unsigned    dim       = 2;
    // degrees of lowest-order TH element
    const unsigned    fieldDegU = 2; 
    const unsigned    fieldDegP = 1;

    const unsigned    tiOrder   = 1;
    typedef  base::time::BDF<tiOrder> MSM;
    
    const base::Shape shape     = base::SimplexShape<dim>::value;
    const base::Shape surfShape = base::SimplexShape<dim-1>::value;

    //--------------------------------------------------------------------------
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input.dat \n\n";
        return -1;
    }

    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    //--------------------------------------------------------------------------
    std::string meshFile, surfFile;
    double viscosity, density, tolerance, penaltyFactor, stepSize;
    unsigned maxIter, numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "surfFile",         surfFile );
        prop.registerPropertiesVar( "viscosity",        viscosity );
        prop.registerPropertiesVar( "density",          density );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "tolerance",        tolerance );
        prop.registerPropertiesVar( "penaltyFactor",    penaltyFactor );
        prop.registerPropertiesVar( "stepSize",         stepSize );
        prop.registerPropertiesVar( "numSteps",         numSteps );

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

    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    //--------------------------------------------------------------------------
    // Surface mesh
    typedef base::Unstructured<surfShape,1,dim>    SurfMesh;

    SurfMesh surfMesh;
    {
        std::ifstream smf( surfFile.c_str() );
        base::io::smf::readMesh( smf, surfMesh );
        smf.close();
    }

    //--------------------------------------------------------------------------
    // Compute the level set data
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    const bool isSigned = true;
    base::cut::bruteForce( mesh, surfMesh, isSigned, levelSet );

    const unsigned kernelDegEstimate = 5;

    //--------------------------------------------------------------------------
    // Make cut cell structure
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells;
    base::cut::generateCutCells( mesh, levelSet, cells );

    // Quadrature 
    typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
    CutQuadrature quadrature( cells, true );

    // for surface
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    //------------------------------------------------------------------------------
    // Finite element bases
    const unsigned    nHist = MSM::numSteps;
    
    const unsigned    doFSizeU = dim;
    typedef base::fe::Basis<shape,fieldDegU>                FEBasisU;
    typedef base::cut::ScaledField<FEBasisU,doFSizeU,nHist> Velocity;
    Velocity velocity;
    base::dof::generate<FEBasisU>( mesh, velocity );
    
    const unsigned    doFSizeP = 1;
    typedef base::fe::Basis<shape,fieldDegP>                FEBasisP;
    typedef base::cut::ScaledField<FEBasisP,doFSizeP,nHist> Pressure;
    Pressure pressure;
    base::dof::generate<FEBasisP>( mesh, pressure );

    const unsigned    doFSizeS = dim;
    typedef base::fe::Basis<surfShape,1>              FEBasisS;
    typedef base::Field<FEBasisS,doFSizeS>            SurfField;
    SurfField surfVelocity, surfForces;
    base::dof::generate<FEBasisS>( surfMesh, surfVelocity );
    base::dof::generate<FEBasisS>( surfMesh, surfForces   );

    // set initial condition to the identity 
    base::dof::setField( surfMesh, surfVelocity,
                         boost::bind( &surfaceVelocity<dim,
                                      SurfField::DegreeOfFreedom>, _1, _2 ) );


    // boundary datum    
    base::cut::TransferSurfaceDatum<SurfMesh,SurfField,Mesh::Element>
        s2d( surfMesh, surfVelocity, levelSet );
    
    //--------------------------------------------------------------------------
    //  surface mesh
    typedef base::mesh::BoundaryMeshBinder<Mesh>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh, immersedMesh;

    // from boundary
    {
        // identify list of element boundary faces
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        // generate a mesh from that list (with a filter)
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh,
                                          boost::bind( &boundaryFilter<dim>, _1 ) );

        // make a surface mesh from the implicit surface
        base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, immersedMesh );
    }

    // the composite field with geometry, velocity and pressure
    typedef base::asmb::FieldBinder<Mesh,Velocity,Pressure> Field;
    Field field( mesh, velocity, pressure );

    // define the system blocks (U,U), (U,P), and (P,U)
    typedef Field::TupleBinder<1,1,1>::Type TopLeft;
    typedef Field::TupleBinder<1,2>::Type   TopRight;
    typedef Field::TupleBinder<2,1>::Type   BotLeft;


    std::vector<double> supportsU, supportsP;
    std::size_t numDoFsU = std::distance( velocity.doFsBegin(), velocity.doFsEnd() );
    supportsU.resize( numDoFsU );
    std::size_t numDoFsP = std::distance( pressure.doFsBegin(), pressure.doFsEnd() );
    supportsP.resize( numDoFsP );

    base::cut::supportComputation( mesh, velocity, quadrature, supportsU );
    base::cut::supportComputation( mesh, pressure, quadrature, supportsP );

    velocity.scaleAndTagBasis( supportsU, 1.e-8 );
    pressure.scaleAndTagBasis( supportsP, 1.e-8 ); 
    //velocity.tagBasis( supportsU, 1.e-8 );
    //pressure.tagBasis( supportsP, 1.e-8 ); 

    // Fix one pressure dof
    // Pressure::DoFPtrIter pIter = pressure.doFsBegin();
    // std::advance( pIter, std::distance( pressure.doFsBegin(), pressure.doFsEnd() )/5 );
    // (*pIter) -> constrainValue( 0, 0.0 );

    // Number of DoFs after constraint application!
    numDoFsU =
        base::dof::numberDoFsConsecutively( velocity.doFsBegin(), velocity.doFsEnd() );
    std::cout << "# Number of velocity dofs " << numDoFsU << std::endl;

    numDoFsP =
        base::dof::numberDoFsConsecutively( pressure.doFsBegin(), pressure.doFsEnd(),
            numDoFsU );
    std::cout << "# Number of pressure dofs " << numDoFsP << std::endl;

    // kernels
    typedef fluid::StressDivergence<  TopLeft::Tuple>  StressDivergence;
    typedef fluid::Convection<        TopLeft::Tuple>  Convection;
    typedef fluid::PressureGradient<  TopRight::Tuple> GradP;
    typedef fluid::VelocityDivergence<BotLeft::Tuple>  DivU;

    StressDivergence stressDivergence( viscosity );
    Convection       convection(       density );
    GradP            gradP;
    DivU             divU( true );

    // for surface fields
    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Velocity,Pressure> SurfaceFieldBinder;
    typedef SurfaceFieldBinder::TupleBinder<1,1,1>::Type STBUU;
    typedef SurfaceFieldBinder::TupleBinder<1,2  >::Type STBUP;
    SurfaceFieldBinder   boundaryFieldBinder(  boundaryMesh, velocity, pressure );
    SurfaceFieldBinder   immersedFieldBinder(  immersedMesh, velocity, pressure );

    std::ofstream forces( "forces.dat" );

    for ( unsigned step = 0; step < numSteps; step++ ) {

        const double time = step * stepSize;
        const double factor = ( time < 1.0 ? time : 1.0 );
    
        std::cout << step << ":  time=" << time << ", factor=" << factor
                  << "\n";

        //base::dof::clearDoFs( velocity );
        //base::dof::clearDoFs( pressure );
        
        //--------------------------------------------------------------------------
        // Nonlinear iterations
        unsigned iter = 0;
        while( iter < maxIter ) {

            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( numDoFsU + numDoFsP );

            std::cout << "* Iteration " << iter << std::flush;

            // compute inertia terms, d/dt, due to time integration
            base::time::computeInertiaTerms<TopLeft,MSM>( quadrature, solver,
                                                          field, stepSize, step,
                                                          density );
    
            // Compute system matrix
            base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver,
                                                             field, stressDivergence );

            base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver,
                                                             field, convection );

            base::asmb::stiffnessMatrixComputation<TopRight>( quadrature, solver,
                                                              field, gradP );

            base::asmb::stiffnessMatrixComputation<BotLeft>( quadrature, solver,
                                                             field, divU );
            // compute residual forces
            base::asmb::computeResidualForces<TopLeft >( quadrature, solver, field,
                                                         stressDivergence );
            base::asmb::computeResidualForces<TopLeft >( quadrature, solver, field, convection );
            base::asmb::computeResidualForces<TopRight>( quadrature, solver, field, gradP );
            base::asmb::computeResidualForces<BotLeft >( quadrature, solver, field, divU );
        
            // Parameter classes
            base::nitsche::OuterBoundary ob( viscosity );
            base::nitsche::ImmersedBoundary<Cell> ib( viscosity, cells );

            // Penalty method
            base::nitsche::penaltyLHS<STBUU>( surfaceQuadrature, solver,
                                              boundaryFieldBinder, ob, penaltyFactor );
        
            base::nitsche::penaltyRHS<STBUU>( surfaceQuadrature, solver, boundaryFieldBinder, 
                                              boost::bind( &dirichlet<dim>, _1, factor),
                                              ob, penaltyFactor );

            base::nitsche::penaltyLHS<STBUU>( surfaceQuadrature, solver,
                                              immersedFieldBinder, ib, penaltyFactor );
        
            base::nitsche::penaltyRHS2<STBUU>( surfaceQuadrature, solver, immersedFieldBinder,
                                               s2d, ib, penaltyFactor );

            // Nitsche terms
            base::nitsche::primalEnergyLHS<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                                   boundaryFieldBinder, ob );
            base::nitsche::dualEnergyLHS<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                                 boundaryFieldBinder, ob );
            base::nitsche::energyRHS<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                             boundaryFieldBinder,
                                             boost::bind( &dirichlet<dim>, _1, factor),
                                             ob );
        
            base::nitsche::primalEnergyLHS<STBUP>( gradP, surfaceQuadrature, solver,
                                                   boundaryFieldBinder, ob );
            base::nitsche::dualEnergyLHS<STBUP>( gradP, surfaceQuadrature, solver,
                                                 boundaryFieldBinder, ob );
         
            base::nitsche::energyRHS<STBUP>( gradP, surfaceQuadrature, solver,
                                             boundaryFieldBinder,
                                             boost::bind( &dirichlet<dim>, _1, factor), ob );

            base::nitsche::energyResidual<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                                  boundaryFieldBinder, ob );
            
            base::nitsche::energyResidual<STBUP>( gradP, surfaceQuadrature, solver,
                                                  boundaryFieldBinder, ob );

            base::nitsche::primalEnergyLHS<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                                   immersedFieldBinder, ib );
            base::nitsche::dualEnergyLHS<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                                 immersedFieldBinder, ib );

            base::nitsche::energyRHS2<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                              immersedFieldBinder, s2d, ib );
            
            base::nitsche::primalEnergyLHS<STBUP>( gradP, surfaceQuadrature, solver,
                                                   immersedFieldBinder, ib );
            base::nitsche::dualEnergyLHS<STBUP>( gradP, surfaceQuadrature, solver,
                                                 immersedFieldBinder, ib );

            base::nitsche::energyRHS2<STBUP>( gradP, surfaceQuadrature, solver,
                                              immersedFieldBinder, s2d, ib );

            base::nitsche::energyResidual<STBUU>( stressDivergence, surfaceQuadrature, solver,
                                                  immersedFieldBinder, ib );

            base::nitsche::energyResidual<STBUP>( gradP, surfaceQuadrature, solver,
                                                  immersedFieldBinder, ib );
        
            // Finalise assembly
            solver.finishAssembly();

            // check convergence via solver norms
            const double residualNorm = solver.norm();
            std::cout << " |R| = " << residualNorm << std::flush;

            if ( residualNorm < tolerance * viscosity) {
                std::cout << std::endl;
                break;
            }

            // Solve
            solver.superLUSolve();

            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, velocity );
            base::dof::addToDoFsFromSolver( solver, pressure );
            //base::dof::setDoFsFromSolver( solver, pressure );
        
            // check convergence via solver norms
            const double incrementNorm = solver.norm(0, numDoFsU );
            std::cout << " |dU| = " << incrementNorm << std::endl;

            // push history
            base::dof::pushHistory( velocity );
            base::dof::pushHistory( pressure );
        
            if ( incrementNorm < tolerance ) break;

            iter++;

        }

        writeVTKFile( baseName, step, mesh, velocity, pressure, levelSet, viscosity );
        {
            base::Vector<dim>::Type sumOfForces = base::constantVector<dim>( 0. );
            
            typedef Field::TupleBinder<1,2>::Type UP;
            //typedef fluid::Stress<UP::Tuple> Stress;
            //Stress stress( viscosity );
            typedef fluid::Traction<UP::Tuple> Traction;
            Traction traction( viscosity );
                
            base::cut::ComputeSurfaceForces<SurfMesh,SurfField,
                                            SurfaceQuadrature,STBUP::Tuple,Traction>
                computeSurfaceForces( surfMesh, surfForces, surfaceQuadrature, levelSet, traction );

            SurfaceFieldBinder::FieldIterator first = immersedFieldBinder.elementsBegin();
            SurfaceFieldBinder::FieldIterator  last = immersedFieldBinder.elementsEnd();
            for ( ; first != last; ++first ) {
                
                sumOfForces +=
                    computeSurfaceForces( STBUP::makeTuple( *first ) );
                
            }

            writeSurfaceVTKFile( baseName, step, surfMesh, surfVelocity, surfForces );


            std::cout << "  F= " << sumOfForces.transpose() << " \n";

            forces << time << " " << sumOfForces.transpose() << std::endl;;
        }

    }

    forces.close();

    return 0;
}
