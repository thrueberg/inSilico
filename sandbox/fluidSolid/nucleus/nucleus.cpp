#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/verify.hpp>
#include <base/mesh/Size.hpp>
#include <base/auxi/BoundingBox.hpp>
#include <base/cut/ImplicitFunctions.hpp>

#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/bruteForce.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/stabiliseBasis.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/location.hpp>
#include <base/dof/constrainBoundary.hpp>

#include <base/asmb/SimpleIntegrator.hpp>
#include <base/kernel/FieldIntegral.hpp> 

#include <base/solver/Eigen3.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>
#include <base/post/ParticleTracer.hpp>

#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include <surf/Moments.hpp> 

//------------------------------------------------------------------------------
// Local includes
#include "../findLocation.hpp"

#include "../extrapolateToOutside.hpp"
#include "../advectField.hpp"
#include "../moveSurface.hpp"

#include "../writeVTK.hpp"

#include "../BoundaryConditions.hpp"
#include "../TypesAndAttributes.hpp"

#include "../Fluid.hpp"
#include "../Solid.hpp"
#include "../SolidFluidInterface.hpp"
#include "../FluidFluidInterface.hpp"
#include "../generateMesh.hpp"

#include "InputNucleus.hpp"


//==============================================================================
//! Write data to VTK
template<typename MESH, typename DISP, typename VELOC, typename PRESS>
void writeVTKFile( const std::string& baseName,
                   const unsigned step,
                   const MESH&    mesh,
                   const DISP&    displacement,
                   const VELOC&   velocityCS,
                   const PRESS&   pressureCS,
                   const VELOC&   velocityGel,
                   const PRESS&   pressureGel, 
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSetN,
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSetM,
                   const double dt )
{
    VTKWriter<MESH> vtkWriter( mesh, baseName, step );

    vtkWriter.writeDistances( levelSetN,         "distNucl" );
    vtkWriter.writeDistances( levelSetM,         "distMemb" );
    vtkWriter.writeField( displacement,          "disp" );
    vtkWriter.writeField( velocityCS,            "velocityCS" );
    vtkWriter.writeField( pressureCS,            "pressureCS" );
    vtkWriter.writeField( velocityGel,           "velocityGel" );
    vtkWriter.writeField( pressureGel,           "pressureGel" );
    vtkWriter.writePreviousField( displacement,  "prevDisp" );

    // generate a global velocity field
    static const unsigned doFSize = DISP::DegreeOfFreedom::size;
    typedef typename base::Vector<doFSize>::Type VecDoF;
    std::vector<VecDoF> globalVelocity;
    {
        std::vector<VecDoF> disp, velCS, velGel;
        base::post::evaluateFieldAtNodes( mesh, displacement, disp );
        base::post::evaluateFieldAtNodes( mesh, velocityCS,   velCS );
        base::post::evaluateFieldAtNodes( mesh, velocityGel,  velGel );

        std::vector<double> dN, dM;
        std::transform( levelSetN.begin(), levelSetN.end(),
                        std::back_inserter( dN ),
                        boost::bind(
                            &base::cut::LevelSet<MESH::Node::dim>::getSignedDistance, _1 ) );
        std::transform( levelSetM.begin(), levelSetM.end(),
                        std::back_inserter( dM ),
                        boost::bind(
                            &base::cut::LevelSet<MESH::Node::dim>::getSignedDistance, _1 ) );

        for ( std::size_t s = 0; s < disp.size(); s++ ) {
            const VecDoF nuclU = disp[s]/dt;
            const VecDoF globalU = (dN[s] >= 0. ? nuclU :
                                    (dM[s] >= 0. ? velCS[s] : velGel[s] ) );
            globalVelocity.push_back( globalU );
        }
    }
    vtkWriter.writePointData( globalVelocity.begin(),
                              globalVelocity.end(), "globalU" );
}

//------------------------------------------------------------------------------
//! Write the VTK file of surface mesh
template<typename SMESH>
void writeSurfaceVTKFile( const std::string& name,
                          const unsigned num, 
                          const SMESH& immersedSurface )
{
    VTKWriter<SMESH> vtkWriter( immersedSurface, name, num );
}

//==============================================================================
//------------------------------------------------------------------------------
// Compute the enclosed volume, the centroid and principal directions
template<typename BOUNDMESH, typename SURFMESH, typename SQUAD, typename SCUTQUAD>
void surfaceFeatures( const BOUNDMESH& boundaryMesh, 
                      const SURFMESH&  membrane,
                      const SURFMESH&  envelope, 
                      const SQUAD&     surfaceQuadrature,
                      const SCUTQUAD&  surfaceCutQuadrature, 
                      std::ostream& out )
{
    static const unsigned dim = SURFMESH::Node::dim;
    typedef typename base::Vector<dim>::Type     VecDim;
    typedef typename base::Matrix<dim,dim>::Type MatDimDim;
 
    // get shape features of the membrane
    const double volumeM =
        surf::enclosedVolume( membrane, surfaceQuadrature ) +
        surf::enclosedVolume( boundaryMesh, surfaceCutQuadrature );

    const VecDim momentM =
        surf::enclosedVolumeMoment( membrane, surfaceQuadrature ) +      
        surf::enclosedVolumeMoment( boundaryMesh, surfaceCutQuadrature );

    const MatDimDim secondMomentM =
        surf::enclosedVolumeSecondMoment( membrane, surfaceQuadrature ) +      
        surf::enclosedVolumeSecondMoment( boundaryMesh, surfaceCutQuadrature );

    const VecDim centroidM = momentM / volumeM;
     
    MatDimDim inertiaM;
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            inertiaM(d1, d2) = secondMomentM(d1,d2) - volumeM * centroidM[d1] * centroidM[d2];
        }
    }
         
    // compute principal values of inertia and the angle-axis rotation
    MatDimDim evecM;
    const VecDim princValM = base::eigenPairs( inertiaM, evecM );
    const std::pair<double,base::Vector<3>::Type> aaM = base::angleAxis( evecM );
     
    out << volumeM << "  " << centroidM.transpose() << "  "
        << princValM.transpose() << "  "
        << aaM.first << "  "  << (aaM.second).transpose() << "  ";

    //--------------------------------------------------------------------------
    // get shape features of the envelope
    const double volumeE =
        surf::enclosedVolume( envelope, surfaceQuadrature );

    const VecDim momentE =
        surf::enclosedVolumeMoment( envelope, surfaceQuadrature );

    const MatDimDim secondMomentE =
        surf::enclosedVolumeSecondMoment( envelope, surfaceQuadrature );

    const VecDim centroidE = momentE / volumeE;
     
    MatDimDim inertiaE;
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            inertiaE(d1, d2) = secondMomentE(d1,d2) - volumeE * centroidE[d1] * centroidE[d2];
        }
    }
         
    // compute principal values of inertia and the angle-axis rotation
    MatDimDim evecE;
    const VecDim princValE = base::eigenPairs( inertiaE, evecE );
    const std::pair<double,base::Vector<3>::Type> aaE = base::angleAxis( evecE );
     
    out << volumeE << "  " << centroidE.transpose() << "  "
        << princValE.transpose() << "  "
        << aaE.first << "  "  << (aaE.second).transpose();

}

//==============================================================================
// MAIN
int main( int argc, char* argv[] )
{
    // spatial dimension
    const unsigned    dim = SPACEDIM; 

    typedef base::solver::Eigen3           Solver;
    
    // Check input arguments
    if ( argc != 2 ) {
        std::cerr << "Usage: " << argv[0] << " input.dat\n"
                  << "(Compiled for dim=" << dim << ")\n\n";
        return -1;
    }
    
    // read name of input file
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[1] );
    const InputNucleus<dim> ui( inputFile );

    // Simulation attributes
    const bool simplex = false;
    const bool stokesStabil = false;
    typedef mat::hypel::NeoHookeanCompressible Material;
    typedef TypesAndAttributes<dim,simplex>    TA;

    // basic attributes of the computation
    base::auxi::BoundingBox<dim> bbox( ui.bbmin, ui.bbmax );
    
    TA::Mesh mesh;
    generateMesh( mesh, ui );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    TA::BoundaryMesh boundaryMesh;
    base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                      meshBoundary.end(),
                                      mesh, boundaryMesh );


    // mesh size (all equal, more or less)
    const double h = base::mesh::Size<TA::Mesh::Element>::apply( mesh.elementPtr(0) );

    //--------------------------------------------------------------------------
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSetNucleus, levelSetMembrane;
    base::cut::analyticLevelSet( mesh,
                                 base::cut::Sphere<dim>( ui.RN, ui.centerN ),
                                 true, levelSetNucleus );
    base::cut::analyticLevelSet( mesh,
                                 base::cut::Sphere<dim>( ui.RM, ui.centerM ),
                                 true, levelSetMembrane );

    //--------------------------------------------------------------------------
    // make cut cell structures
    std::vector<TA::Cell> cellsNucleus, cellsMembrane, cellsCS;
    base::cut::generateCutCells( mesh, levelSetNucleus,  cellsNucleus );
    base::cut::generateCutCells( mesh, levelSetMembrane, cellsMembrane );
    
    {
        cellsCS = cellsNucleus;
        std::for_each( cellsCS.begin(), cellsCS.end(),
                       boost::bind( &TA::Cell::reverse, _1 ) );

        base::cut::generateCutCells( mesh, levelSetMembrane, cellsCS,
                                     base::cut::INTERSECT );
    }

    // intersection with the box boundary
    std::vector<TA::SurfCell> surfaceCells;
    base::cut::generateCutCells( boundaryMesh, levelSetMembrane, surfaceCells );

    // Elastic material
    Material material( mat::Lame::lambda( ui.E, ui.nu),
                       mat::Lame::mu(     ui.E, ui.nu ) );

    // Solid and fluid handlers
    typedef Solid<TA::Mesh, Material>                        Solid;
    typedef Fluid<TA::Mesh, stokesStabil>                    Fluid;
    typedef SolidFluidInterface<TA::SurfaceMesh,Solid,Fluid> SFI;
    typedef FluidFluidInterface<TA::SurfaceMesh,Fluid>       FFI;
    
    Solid nucleus(      mesh, material );
    Fluid gel(          mesh, ui.viscosityGel, ui.alpha );
    Fluid cytoSkeleton( mesh, ui.viscosityCS,  ui.alpha );

    // find geometry association for the dofs
    std::vector<std::pair<std::size_t,TA::VecDim> > doFLocationD, doFLocationU;
    base::dof::associateLocation( nucleus.getDisplacement(), doFLocationD );
    base::dof::associateLocation( gel.getVelocity(),       doFLocationU );

    // Quadrature along a surface
    TA::CutQuadrature     cutQuadratureNucleus( cellsNucleus,  true  ); 
    TA::CutQuadrature     cutQuadratureGel(     cellsMembrane, false );
    TA::CutQuadrature     cutQuadratureCS(      cellsCS,       true  );
    TA::SurfaceQuadrature surfaceQuadrature;
    TA::SurfCutQuadrature surfaceCutQuadrature( surfaceCells, true );

    //--------------------------------------------------------------------------
    // Run a zero-th step in order to write the data before computation
    // (store the previously enclosed volume)
    double prevVolumeN, initialVolumeCS;
    {
        std::cout << 0 << "  " << 0. << "  0  " << std::flush;

        writeVTKFile( "nucleus", 0, mesh,
                      nucleus.getDisplacement(),
                      cytoSkeleton.getVelocity(), cytoSkeleton.getPressure(),
                      gel.getVelocity(),          gel.getPressure(),
                      levelSetNucleus, levelSetMembrane,  ui.dt );

        // for surface field
        TA::SurfaceMesh membrane, envelope;
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cellsMembrane, membrane );
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cellsNucleus,  envelope );

        writeSurfaceVTKFile( "membrane", 0, membrane );
        writeSurfaceVTKFile( "envelope", 0, envelope );
            
        surfaceFeatures( boundaryMesh, membrane, envelope,
                         surfaceQuadrature, surfaceCutQuadrature, std::cout );
        std::cout << std::endl;

        // store the contained volume
        prevVolumeN     = surf::enclosedVolume( envelope, surfaceQuadrature );
        initialVolumeCS =
            surf::enclosedVolume( membrane, surfaceQuadrature ) +
            surf::enclosedVolume( boundaryMesh, surfaceCutQuadrature );
    }

    //--------------------------------------------------------------------------
    // * Loop over time steps *
    //--------------------------------------------------------------------------
    const double supportThreshold = std::numeric_limits<double>::min();
    for ( unsigned step = 0; step < ui.numLoadSteps; step++ ) {
        
        // for surface field
        TA::SurfaceMesh membrane, envelope;
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cellsMembrane, membrane );
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cellsNucleus,  envelope );

        // Interface handler
        SFI envelopeHandler( envelope,
                             nucleus.getDisplacement(),
                             cytoSkeleton.getVelocity(), cytoSkeleton.getPressure(),
                             ui.viscosityCS );
        
        FFI membraneHandler( membrane,
                             cytoSkeleton.getVelocity(), cytoSkeleton.getPressure(),
                             gel.getVelocity(), gel.getPressure(),
                             cellsMembrane, 
                             ui.viscosityCS, ui.viscosityGel, ui.sigma );

        std::cout << step+1 << "  " << (step+1) * ui.dt << "  " << std::flush;

        //----------------------------------------------------------------------
        // 0) Solve Lagrangian problem
        
        // compute supports
        std::vector<double> supportsD, supportsUCS, supportsPCS, supportsUGel, supportsPGel;
        base::cut::supportComputation( mesh, nucleus.getDisplacement(),
                                       cutQuadratureNucleus, supportsD );
        base::cut::supportComputation( mesh, cytoSkeleton.getVelocity(),
                                       cutQuadratureCS, supportsUCS );
        base::cut::supportComputation( mesh, cytoSkeleton.getPressure(),
                                       cutQuadratureCS, supportsPCS );
        base::cut::supportComputation( mesh, gel.getVelocity(),
                                       cutQuadratureGel, supportsUGel );
        base::cut::supportComputation( mesh, gel.getPressure(),
                                       cutQuadratureGel, supportsPGel );

        // Hard-wired Dirichlet constraints (fluid BC)
        {
            // apply to outer material = gel
            base::dof::constrainBoundary<Fluid::FEBasisU>(
                meshBoundary.begin(),
                meshBoundary.end(),
                mesh, gel.getVelocity(),
                boost::bind( &dirichletBC<dim,Fluid::DoFU>,
                             _1, _2, ui.ubar, bbox, ui.bc ) );

            // in case of boundary intersections, apply to cytoskeleton
            base::dof::constrainBoundary<Fluid::FEBasisU>(
                meshBoundary.begin(),
                meshBoundary.end(),
                mesh, cytoSkeleton.getVelocity(),
                boost::bind( &dirichletBC<dim,Fluid::DoFU>,
                             _1, _2, ui.ubar, bbox, ui.bc ) );

            // Fix first pressure dof in case of closed cavity flow
            if ( ui.bc == CAVITY ) {
                Fluid::Pressure::DoFPtrIter pIter = gel.getPressure().doFsBegin();
                (*pIter) -> constrainValue( 0, 0.0 );
            }
        }
    
        // Basis stabilisation
        base::cut::stabiliseBasis( mesh, nucleus.getDisplacement(),  supportsD,    doFLocationD );
        base::cut::stabiliseBasis( mesh, cytoSkeleton.getVelocity(), supportsUCS,  doFLocationU );
        base::cut::stabiliseBasis( mesh, cytoSkeleton.getPressure(), supportsPCS,  doFLocationD );
        base::cut::stabiliseBasis( mesh, gel.getVelocity(),          supportsUGel, doFLocationU );
        base::cut::stabiliseBasis( mesh, gel.getPressure(),          supportsPGel, doFLocationD );
        
        // number DoFs 
        const std::size_t activeDoFsD = 
            base::dof::numberDoFsConsecutively( nucleus.getDisplacement().doFsBegin(),
                                                nucleus.getDisplacement().doFsEnd() );
        const std::size_t activeDoFsUCS = 
            base::dof::numberDoFsConsecutively( cytoSkeleton.getVelocity().doFsBegin(),
                                                cytoSkeleton.getVelocity().doFsEnd(),
                                                activeDoFsD );
        const std::size_t activeDoFsPCS = 
            base::dof::numberDoFsConsecutively( cytoSkeleton.getPressure().doFsBegin(),
                                                cytoSkeleton.getPressure().doFsEnd(),
                                                activeDoFsD + activeDoFsUCS );
        const std::size_t activeDoFsUGel = 
            base::dof::numberDoFsConsecutively( gel.getVelocity().doFsBegin(),
                                                gel.getVelocity().doFsEnd(),
                                                activeDoFsD + activeDoFsUCS + activeDoFsPCS );
        const std::size_t activeDoFsPGel = 
            base::dof::numberDoFsConsecutively( gel.getPressure().doFsBegin(),
                                                gel.getPressure().doFsEnd(),
                                                activeDoFsD + activeDoFsUCS + activeDoFsPCS +
                                                activeDoFsUGel );


        // start from 0 for the increments, set rest to zero too
        base::dof::clearDoFs( nucleus.getDisplacement() );
        base::dof::clearDoFs( cytoSkeleton.getPressure() );
        base::dof::clearDoFs( cytoSkeleton.getVelocity() );
        base::dof::clearDoFs( gel.getPressure() );
        base::dof::clearDoFs( gel.getVelocity() );

        //----------------------------------------------------------------------
        // Nonlinear iterations
        unsigned iter = 0;
        for ( ; iter < ui.maxIter; iter++ ) {

            // Create a solver object
            Solver solver( activeDoFsD + activeDoFsUCS + activeDoFsPCS +
                           activeDoFsUGel + activeDoFsPGel );

            // register to solver
            nucleus.registerInSolver(      solver );
            cytoSkeleton.registerInSolver( solver );
            gel.registerInSolver(          solver );
            membraneHandler.registerInSolver(     solver );
            envelopeHandler.registerInSolver(     solver );

            // Compute as(d,dd) and Da(d,dd)[Delta d]
            nucleus.assembleBulk( cutQuadratureNucleus, solver, iter );

            // Body force
            nucleus.bodyForce( cutQuadratureNucleus, solver,
                               boost::bind( unitForce<0,dim>, _1, ui.forceValue ) );

            // Computate af(u,p; du,dp)
            cytoSkeleton.assembleBulk( cutQuadratureCS,  solver );
            gel.assembleBulk(          cutQuadratureGel, solver );

            // Penalty terms for int_Gamma (dot(d) - u) (E delta d - mu delta u) ds
            envelopeHandler.assemblePenaltyTerms( surfaceQuadrature, solver, ui.penaltyFac, ui.dt );
            membraneHandler.assemblePenaltyTerms( surfaceQuadrature, solver, ui.penaltyFac );

            // Boundary energy terms (beta = 0)  [ -int_Gamma (...) ds ]
            envelopeHandler.assembleEnergyTerms( surfaceQuadrature, solver, ui.dt );
            membraneHandler.assembleEnergyTerms( surfaceQuadrature, solver );

            // surface tension
            membraneHandler.surfaceTension( surfaceQuadrature, solver );
            

            // Finalise assembly
            solver.finishAssembly();

            // norm of residual
            const double conv1 = solver.norm(0, activeDoFsD) / ui.E;

            if ( isnan( conv1 ) ) return 1;
            
            if ( ( iter > 0 )  and ( conv1 < ui.tolerance ) ) {
                std::cout << "schon" << std::endl;
                break;
            }

            // Solve
#ifdef LOAD_PARDISO
            solver.pardisoLUSolve();
#else

#ifdef LOAD_UMFPACK
            solver.umfPackLUSolve();
#else
            solver.superLUSolve();
            //solver.biCGStabSolve();
#endif
            
#endif
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, nucleus.getDisplacement() ); //<>
            base::dof::setDoFsFromSolver(   solver, cytoSkeleton.getVelocity() );
            base::dof::setDoFsFromSolver(   solver, cytoSkeleton.getPressure() );
            base::dof::setDoFsFromSolver(   solver, gel.getVelocity() );
            base::dof::setDoFsFromSolver(   solver, gel.getPressure() );

            const double conv2 = solver.norm(0, activeDoFsD);

            if ( conv2 < ui.tolerance ) break;
            
        } // finish non-linear iteration
        std::cout << iter << "  " << std::flush;

        // write a vtk file
        writeVTKFile( "nucleus", step+1, mesh,
                      nucleus.getDisplacement(),
                      cytoSkeleton.getVelocity(), cytoSkeleton.getPressure(),
                      gel.getVelocity(),          gel.getPressure(),
                      levelSetNucleus, levelSetMembrane, ui.dt );
        writeSurfaceVTKFile( "membrane", step+1, membrane );
        writeSurfaceVTKFile( "envelope", step+1, envelope );

        //----------------------------------------------------------------------
        // 1) Pass Data to complementary domain of solid
        {
            const double extL = h; // extrapolation length

            // pass extrapolated solution to inactive DoFs
            extrapolateToFictitious( mesh, nucleus.getDisplacement(),
                                     extL, doFLocationD, levelSetNucleus, bbox );
        }

        //----------------------------------------------------------------------
        // 2) Geometry update
        // a) Nucleus
        {
            // get shape features
            const double volumeN = surf::enclosedVolume( envelope, surfaceQuadrature );
            const TA::VecDim momentN =
                surf::enclosedVolumeMoment( envelope, surfaceQuadrature );
            const TA::VecDim centroidN = momentN / volumeN;

            const double factorN = std::pow( prevVolumeN / volumeN,
                                            1./static_cast<double>( dim ) );
     
            // Move with velocity solution
            moveSurface( mesh, cytoSkeleton.getVelocity(),
                         envelope, ui.dt, factorN, centroidN );

            // compute new level set for envelope of nucleus
            base::cut::bruteForce( mesh, envelope, true, levelSetNucleus );

            // update the cut cell structure
            base::cut::generateCutCells( mesh, levelSetNucleus, cellsNucleus );


            // store the contained volume
            prevVolumeN = surf::enclosedVolume( envelope, surfaceQuadrature );
        }

        // b) CS
        {
            // get shape features
            const double volumeCS =
                surf::enclosedVolume( membrane,  surfaceQuadrature ) +
                surf::enclosedVolume( boundaryMesh, surfaceCutQuadrature );
            const TA::VecDim momentCS =
                surf::enclosedVolumeMoment( membrane, surfaceQuadrature ) +
                surf::enclosedVolumeMoment( boundaryMesh, surfaceCutQuadrature );
            const TA::VecDim centroidCS = momentCS / volumeCS;

            // Move with velocity solution
            moveSurface( mesh, gel.getVelocity(), membrane, ui.dt, 1.0, centroidCS );
            // rescale in order to preserve initial volume
            rescaleSurface( membrane, initialVolumeCS, volumeCS, centroidCS );

            // compute new level set for membrane
            base::cut::bruteForce( mesh, membrane, true, levelSetMembrane );

            // update the cut cell structure
            base::cut::generateCutCells( mesh, levelSetMembrane, cellsMembrane );
            base::cut::generateCutCells( boundaryMesh, levelSetMembrane, surfaceCells );
        }

        // new structure for cytoskeleton
        {
            cellsCS = cellsNucleus;
            std::for_each( cellsCS.begin(), cellsCS.end(),
                           boost::bind( &TA::Cell::reverse, _1 ) );
            
            base::cut::generateCutCells( mesh, levelSetMembrane, cellsCS,
                                         base::cut::INTERSECT );
        }


        //----------------------------------------------------------------------
        // 3) advect data
        {
            base::cut::supportComputation( mesh, nucleus.getDisplacement(),
                                           cutQuadratureNucleus, supportsD ); 
        
            // Find the location of the DoFs in the previous configuration
            std::vector<std::pair<std::size_t,TA::VecDim> > previousDoFLocation;
            findPreviousDoFLocations( mesh, nucleus.getDisplacement(), supportsD, 
                                      doFLocationD, previousDoFLocation, bbox, 
                                      supportThreshold, ui.findTolerance, 10 );
            
            // add previous solution to current solution: u_{n+1} = \Delta u + u_n
            {
                Solid::Displacement::DoFPtrIter dIter = nucleus.getDisplacement().doFsBegin();
                Solid::Displacement::DoFPtrIter dEnd  = nucleus.getDisplacement().doFsEnd();
                for ( ; dIter != dEnd; ++dIter ) {
                    for ( unsigned d = 0; d < dim; d++ ) {
                        const double prevU  = (*dIter) -> getHistoryValue<1>( d );
                        const double deltaU = (*dIter) -> getHistoryValue<0>( d );
                        const double currU = prevU + deltaU;
                        (*dIter) -> setValue( d, currU );
                    }
                }
            }
        
            // advect displacement field from previous to new location
            advectField( mesh, nucleus.getDisplacement(), previousDoFLocation, supportsD, 
                         supportThreshold );
        }

        //----------------------------------------------------------------------
        // push history
        base::dof::pushHistory( nucleus.getDisplacement() );

        // remove the linear constraints used in the stabilisation process
        base::dof::clearConstraints( nucleus.getDisplacement() );
        base::dof::clearConstraints( cytoSkeleton.getVelocity()     );
        base::dof::clearConstraints( cytoSkeleton.getPressure()     );
        base::dof::clearConstraints( gel.getVelocity()     );
        base::dof::clearConstraints( gel.getPressure()     );

        surfaceFeatures( boundaryMesh, membrane, envelope,
                         surfaceQuadrature, surfaceCutQuadrature, std::cout );
        std::cout << std::endl;

    } // end load-steps

    return 0;
}
