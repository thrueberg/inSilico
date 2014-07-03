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

#include <surf/Moments.hpp> 

//------------------------------------------------------------------------------
// Local includes
#include "../generateMesh.hpp"
#include "../findLocation.hpp"
#include "../moveSurface.hpp"

#include "../BoundaryConditions.hpp"
#include "../TypesAndAttributes.hpp"

#include "../Fluid.hpp"
#include "../FluidFluidInterface.hpp"

#include "../writeVTK.hpp"
#include "InputFFI.hpp"

//==============================================================================
//! Write data to VTK
template<typename MESH, typename VELOC, typename PRESS>
void writeVTKFile( const std::string& baseName,
                   const unsigned step,
                   const MESH&    mesh,
                   const VELOC&   velocity1,
                   const PRESS&   pressure1, 
                   const VELOC&   velocity2,
                   const PRESS&   pressure2, 
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet )
{
    VTKWriter<MESH> vtkWriter( mesh, baseName, step );

    vtkWriter.writeDistances( levelSet );
    vtkWriter.writeField( velocity1, "velocity1" );
    vtkWriter.writeField( pressure1, "pressure1" );
    vtkWriter.writeField( velocity2, "velocity2" );
    vtkWriter.writeField( pressure2, "pressure2" );
}

//------------------------------------------------------------------------------
//! Write the VTK file of surface mesh
template<typename SMESH,  typename VELOC, typename PRESS>
void writeSurfaceVTKFile( const std::string& name,
                          const unsigned num, 
                          const SMESH& immersedSurface,
                          const VELOC& velocity1,
                          const PRESS& pressure1,
                          const VELOC& velocity2,
                          const PRESS& pressure2,
                          const double viscosity1,
                          const double viscosity2 )
{
    VTKWriter<SMESH> vtkWriter( immersedSurface, name, num );
    
    typedef typename SMESH::Node::VecDim VecDim;

    //--------------------------------------------------------------------------
    // Evaluate displacement, velocity and pressure
    {
        // storage
        std::vector<VecDim> nodalU1, nodalU2;
        std::vector<double> nodalP1, nodalP2;
        // go through all surface elements
        typename SMESH::ElementPtrConstIter eBegin = immersedSurface.elementsBegin();
        typename SMESH::ElementPtrConstIter eEnd   = immersedSurface.elementsEnd();
        for ( ; eBegin != eEnd; ++eBegin ) {
        
            typename SMESH::Element::DomainElement* geomEp =
                (*eBegin) -> getDomainElementPointer();
        
            const std::size_t elemID = geomEp -> getID();

            // get the field elements
            typename VELOC::Element* velocEp1 = velocity1.elementPtr( elemID );
            typename PRESS::Element* pressEp1 = pressure1.elementPtr( elemID );
            typename VELOC::Element* velocEp2 = velocity2.elementPtr( elemID );
            typename PRESS::Element* pressEp2 = pressure2.elementPtr( elemID );

            // go through the local coordinates of the surface element's nodes
            typename SMESH::Element::ParamIter   pIter = (*eBegin) -> parametricBegin();
            typename SMESH::Element::ParamIter   pEnd  = (*eBegin) -> parametricEnd();
            for ( ; pIter != pEnd; ++pIter ) {
                // evaluate and store
                nodalU1.push_back( base::post::evaluateField( geomEp, velocEp1, *pIter ) );
                nodalP1.push_back( base::post::evaluateField( geomEp, pressEp1, *pIter )[0] );
                nodalU2.push_back( base::post::evaluateField( geomEp, velocEp2, *pIter ) );
                nodalP2.push_back( base::post::evaluateField( geomEp, pressEp2, *pIter )[0] );
            }
        }
        
        // pass on to writer
        vtkWriter.writePointData( nodalU1.begin(), nodalU1.end(), "Veloc1" );
        vtkWriter.writePointData( nodalP1.begin(), nodalP1.end(), "Press1" );
        vtkWriter.writePointData( nodalU2.begin(), nodalU2.end(), "Veloc2" );
        vtkWriter.writePointData( nodalP2.begin(), nodalP2.end(), "Press2" );
    }
    
    //--------------------------------------------------------------------------
    // Compute normals and areas
    
    // storage
    std::vector<VecDim> normals;
    std::vector<double> areas;
    {
        // go through all surface elements
        typename SMESH::ElementPtrConstIter eBegin = immersedSurface.elementsBegin();
        typename SMESH::ElementPtrConstIter eEnd   = immersedSurface.elementsEnd();
        for ( ; eBegin != eEnd; ++eBegin ) {

            // compute normal the element centroid
            VecDim normal;
            const double area = 
                base::SurfaceNormal<typename SMESH::Element>()(
                    *eBegin,
                    base::ShapeCentroid<SMESH::Element::shape>::apply(),
                    normal );

            // store
            normals.push_back( normal );
            areas.push_back(   area );
        }
        vtkWriter.writeCellData( normals.begin(),   normals.end(),   "normals"  );
    }

    //--------------------------------------------------------------------------
    // Compute tractions from fluid and solid
    {
        // storages
        std::vector<typename SMESH::Node::VecDim> tractions1, tractions2;

        // bind surface to fluid 1
        typedef base::asmb::SurfaceFieldBinder<const SMESH,const VELOC, const PRESS> SFB;
        SFB sfb1( immersedSurface, velocity1, pressure1 );
        SFB sfb2( immersedSurface, velocity2, pressure2 );

        // types of domain field tuples
        typedef typename SFB::template TupleBinder<1,2>::Type STBUP;
        typedef typename base::asmb::DomainFieldElementPointerTuple<typename STBUP::Tuple>::Type
            DFTUP;
        
        // traction objects
        fluid::Traction<DFTUP> traction1( viscosity1 );
        fluid::Traction<DFTUP> traction2( viscosity2 );

        // collect the sum of the forces (unused)
        VecDim sum1 = base::constantVector<SMESH::Node::dim>( 0. );
        VecDim sum2 = base::constantVector<SMESH::Node::dim>( 0. );

        // go through all surface element tuples
        typename SFB::FieldIterator sEBegin1 = sfb1.elementsBegin();
        typename SFB::FieldIterator sEBegin2 = sfb2.elementsBegin();
        typename SFB::FieldIterator sEEnd    = sfb1.elementsEnd();
        for ( std::size_t ctr = 0; sEBegin1 != sEEnd; ++sEBegin1, ++sEBegin2, ctr++ ) {

            // get local domain coordinate of the centroid of the surface element
            const typename DFTUP::GeomElement::GeomFun::VecDim xi =
                (STBUP::makeTuple( *sEBegin1 ).geomElementPtr()) -> localDomainCoordinate( 
                    base::ShapeCentroid<SMESH::Element::shape>::apply() );

            // fluid domain element tuple
            const DFTUP dftf =
                base::asmb::DomainFieldElementPointerTuple<
                    typename STBUP::Tuple>::convert( *sEBegin1 );

            // compute fluid side traction and store
            const VecDim t1 = traction1( dftf, xi, normals[ ctr ] );
            tractions1.push_back( t1 );

            const VecDim t2 = traction2( dftf, xi, normals[ ctr ] );
            tractions2.push_back( t2 );

            // add to sum (like a mid-point rule integration)
            sum1 += areas[ctr] * t1;
            sum2 += areas[ctr] * t2;
        }
    
        vtkWriter.writeCellData( tractions1.begin(), tractions1.end(), "traction1" );
        vtkWriter.writeCellData( tractions2.begin(), tractions2.end(), "traction2" );
    }
    
    return;
}

//==============================================================================
// Compute the enclosed volume, the centroid and principal directions
template<typename SURFMESH, typename FLUID, typename SQUAD, typename QUAD>
void bubbleFeatures( const SURFMESH& surfaceMesh,
                     FLUID& fluid,
                     const SQUAD& surfaceQuadrature,
                     const QUAD&  quadrature, 
                     std::ostream& out )
{
    static const unsigned dim = FLUID::dim;
    typedef typename base::Vector<dim>::Type     VecDim;
    typedef typename base::Matrix<dim,dim>::Type MatDimDim;
 
    // get shape features
    const double volume =
        surf::enclosedVolume( surfaceMesh, surfaceQuadrature );

    const VecDim moment =
        surf::enclosedVolumeMoment( surfaceMesh, surfaceQuadrature );

    const MatDimDim secondMoment =
        surf::enclosedVolumeSecondMoment( surfaceMesh, surfaceQuadrature );

    // process
    const VecDim centroid = moment / volume;
     
    MatDimDim inertia;
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            inertia(d1, d2) = secondMoment(d1,d2) - volume * centroid[d1] * centroid[d2];
        }
    }
         
    // compute principal values of inertia and the angle-axis rotation
    MatDimDim evec;
    const VecDim princVal = base::eigenPairs( inertia, evec );
    const std::pair<double,base::Vector<3>::Type> aa = base::angleAxis( evec );

    // average rise
    VecDim averageU = base::constantVector<dim>( 0. );
    base::asmb::simplyIntegrate<typename FLUID::UU>(
        quadrature, averageU,
        typename FLUID::FieldBinder( fluid.getMesh(), fluid.getVelocity(), fluid.getPressure() ),
        base::kernel::FieldIntegral<typename FLUID::UU::Tuple>() );
    averageU /= volume;
    
     
    out << volume << "  "
        << centroid.transpose() << "  "
        << averageU.transpose() << "  "
        << princVal.transpose() << "  "
        << aa.first << "  "  << (aa.second).transpose();
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
    const InputFFI<dim> ui( inputFile );

    // Simulation attributes
    const bool simplex      = false;
    const bool stokesStabil = false;
    const bool dynamic      = true; //false; 
    typedef TypesAndAttributes<dim,simplex> TA;

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


    //--------------------------------------------------------------------------
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::analyticLevelSet( mesh,
                                 base::cut::Sphere<dim>( ui.R, ui.center ),
                                 true, levelSet );

    //--------------------------------------------------------------------------
    std::vector<TA::Cell> cells;
    base::cut::generateCutCells( mesh, levelSet, cells );
    std::vector<TA::SurfCell> surfCells;
    base::cut::generateCutCells( boundaryMesh, levelSet, surfCells );

    // Solid and fluid handlers
    typedef Fluid<TA::Mesh,stokesStabil>                       Fluid;
    typedef FluidFluidInterface<TA::SurfaceMesh,Fluid>         FFI;
    
    Fluid fluid1( mesh, ui.viscosity1, ui.alpha, dynamic, ui.rho1, ui.dt  );
    Fluid fluid2( mesh, ui.viscosity2, ui.alpha, dynamic, ui.rho2, ui.dt  );

    // find geometry association for the dofs
    std::vector<std::pair<std::size_t,TA::VecDim> > doFLocationU, doFLocationP;
    base::dof::associateLocation( fluid1.getVelocity(), doFLocationU );
    base::dof::associateLocation( fluid1.getPressure(), doFLocationP );

    // Quadrature along a surface
    TA::CutQuadrature     cutQuadrature1(   cells, true  ); // in
    TA::CutQuadrature     cutQuadrature2(   cells, false ); // out
    TA::SurfaceQuadrature surfaceQuadrature;

    //--------------------------------------------------------------------------
    // Run a zero-th step in order to write the data before computation
    // (store the previously enclosed volume)
    double initialVolume;
    {
        std::cout << 0 << "  " << 0. << "  " << std::flush;

        // for surface field
        TA::SurfaceMesh immersedSurface;
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cells, immersedSurface );

        // Interface handler
        writeVTKFile( "bubble", 0, mesh, 
                      fluid1.getVelocity(), fluid1.getPressure(),
                      fluid2.getVelocity(), fluid2.getPressure(),
                      levelSet );
        writeSurfaceVTKFile( "bubbleSurf", 0, immersedSurface,
                             fluid1.getVelocity(), fluid1.getPressure(),
                             fluid2.getVelocity(), fluid2.getPressure(),
                             ui.viscosity1, ui.viscosity2 );

        bubbleFeatures( immersedSurface, fluid1, surfaceQuadrature, cutQuadrature1, std::cout );
        std::cout << std::endl;

        // store the contained volume
        initialVolume = surf::enclosedVolume( immersedSurface, surfaceQuadrature );

    }

    //--------------------------------------------------------------------------
    // * Loop over time steps *
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < ui.numLoadSteps; step++ ) {
        
        // Interface handler
        TA::SurfaceMesh immersedSurface;
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cells, immersedSurface );

        FFI ffi( immersedSurface,
                 fluid1.getVelocity(), fluid1.getPressure(),
                 fluid2.getVelocity(), fluid2.getPressure(),
                 cells, ui.viscosity1, ui.viscosity2, ui.sigma );
        
        std::cout << step+1 << "  " << (step+1) * ui.dt << "  " << std::flush;

        //----------------------------------------------------------------------
        // 0) Solve Lagrangian problem
        
        // compute supports
        std::vector<double> supportsU1, supportsP1, supportsU2, supportsP2;
        base::cut::supportComputation( mesh, fluid1.getVelocity(), cutQuadrature1, supportsU1 );
        base::cut::supportComputation( mesh, fluid1.getPressure(), cutQuadrature1, supportsP1 );
        base::cut::supportComputation( mesh, fluid2.getVelocity(), cutQuadrature2, supportsU2 );
        base::cut::supportComputation( mesh, fluid2.getPressure(), cutQuadrature2, supportsP2 );

        // Hard-wired Dirichlet constraints (fluid BC)
        {
            base::dof::constrainBoundary<Fluid::FEBasisU>(
                meshBoundary.begin(),
                meshBoundary.end(),
                mesh, fluid2.getVelocity(),
                boost::bind( &dirichletBC<dim,Fluid::DoFU>,
                             _1, _2, ui.ubar, bbox, ui.bc ) );

            // Fix first pressure dof in case of closed cavity flow
            if ( ui.bc == CAVITY ) {
                Fluid::Pressure::DoFPtrIter pIter = fluid2.getPressure().doFsBegin();
                (*pIter) -> constrainValue( 0, 0.0 );
            }

            if ( ui.bc == TANK ) {

                base::dof::constrainBoundary<Fluid::FEBasisP>(
                    meshBoundary.begin(),
                    meshBoundary.end(),
                    mesh, fluid2.getPressure(),
                    boost::bind( &tankP<dim,Fluid::DoFP>, _1, _2, bbox ) );
                
            }
        }
    
        // Basis stabilisation
        base::cut::stabiliseBasis( mesh, fluid1.getVelocity(), supportsU1, doFLocationU );
        base::cut::stabiliseBasis( mesh, fluid1.getPressure(), supportsP1, doFLocationP );
        base::cut::stabiliseBasis( mesh, fluid2.getVelocity(), supportsU2, doFLocationU );
        base::cut::stabiliseBasis( mesh, fluid2.getPressure(), supportsP2, doFLocationP );
        
        // number DoFs 
        const std::size_t activeDoFsU1= 
            base::dof::numberDoFsConsecutively( fluid1.getVelocity().doFsBegin(),
                                                fluid1.getVelocity().doFsEnd(),
                                                0 );
        const std::size_t activeDoFsP1 = 
            base::dof::numberDoFsConsecutively( fluid1.getPressure().doFsBegin(),
                                                fluid1.getPressure().doFsEnd(),
                                                activeDoFsU1 );
        const std::size_t activeDoFsU2= 
            base::dof::numberDoFsConsecutively( fluid2.getVelocity().doFsBegin(),
                                                fluid2.getVelocity().doFsEnd(),
                                                activeDoFsU1 + activeDoFsP1  );
        const std::size_t activeDoFsP2 = 
            base::dof::numberDoFsConsecutively( fluid2.getPressure().doFsBegin(),
                                                fluid2.getPressure().doFsEnd(),
                                                activeDoFsU1 + activeDoFsP1 + activeDoFsU2 );


        // start from 0 for the increments, set rest to zero too
        base::dof::clearDoFs( fluid1.getPressure() );
        base::dof::clearDoFs( fluid1.getVelocity() );
        base::dof::clearDoFs( fluid2.getPressure() );
        base::dof::clearDoFs( fluid2.getVelocity() );

        for ( unsigned iter = 0; iter < ui.fluidIter; iter++ ) {

            //------------------------------------------------------------------
            // Create a solver object
            Solver solver( activeDoFsU1 + activeDoFsP1 + activeDoFsU2 + activeDoFsP2 );

            // register to solver
            fluid1.registerInSolver( solver );
            fluid2.registerInSolver( solver );
            ffi.registerInSolver(    solver );

            // Computate af(u,p; du,dp)
            fluid1.assembleBulk( cutQuadrature1, solver, iter );

            fluid1.bodyForce( cutQuadrature1, solver,
                              boost::bind( unitForce<dim-1,dim>, _1, -0.98*ui.rho1 ) );
            
            fluid2.assembleBulk( cutQuadrature2, solver, iter );
            
            fluid2.bodyForce( cutQuadrature2, solver,
                              boost::bind( unitForce<dim-1,dim>, _1, -0.98*ui.rho2 ) );

            // Penalty terms 
            ffi.assemblePenaltyTerms( surfaceQuadrature, solver, ui.penaltyFac );

            // Boundary energy terms 
            ffi.assembleEnergyTerms( surfaceQuadrature, solver  );

            // Surface tension
            ffi.surfaceTension( surfaceQuadrature, solver );
            
            // Finalise assembly
            solver.finishAssembly();

            // Solve
#ifdef LOAD_PARDISO
            solver.pardisoLUSolve();
#else
            solver.superLUSolve();
#endif
        
            // distribute results back to dofs
            base::dof::setDoFsFromSolver(   solver, fluid1.getVelocity() );
            base::dof::setDoFsFromSolver(   solver, fluid1.getPressure() );
            base::dof::setDoFsFromSolver(   solver, fluid2.getVelocity() );
            base::dof::setDoFsFromSolver(   solver, fluid2.getPressure() );

        } // finished non-linear iterations

        
        // write a vtk file
        writeVTKFile( "bubble", step+1, mesh, 
                      fluid1.getVelocity(), fluid1.getPressure(),
                      fluid2.getVelocity(), fluid2.getPressure(),
                      levelSet );
        writeSurfaceVTKFile( "bubbleSurf", step+1, immersedSurface,
                             fluid1.getVelocity(), fluid1.getPressure(),
                             fluid2.getVelocity(), fluid2.getPressure(),
                             ui.viscosity1, ui.viscosity2 );

        //----------------------------------------------------------------------
        // 2) Geometry update
        {
            // get shape features
            const double volume =
                surf::enclosedVolume( immersedSurface, surfaceQuadrature );
            const TA::VecDim moment =
                surf::enclosedVolumeMoment( immersedSurface, surfaceQuadrature );
            const TA::VecDim centroid = moment / volume;

            const double factor = 1.0;
     
            // Move with velocity solution
            moveSurface( mesh, fluid1.getVelocity(), immersedSurface, ui.dt, factor, centroid );

            rescaleSurface( immersedSurface, initialVolume, volume, centroid );
            
            base::cut::bruteForce( mesh, immersedSurface, true, levelSet );

            // update the cut cell structure
            base::cut::generateCutCells( mesh,         levelSet, cells );
            base::cut::generateCutCells( boundaryMesh, levelSet, surfCells );


        }


        //----------------------------------------------------------------------
        // remove the linear constraints used in the stabilisation process
        base::dof::pushHistory( fluid1.getVelocity() );
        base::dof::pushHistory( fluid2.getVelocity() );

        base::dof::clearConstraints( fluid1.getVelocity()     );
        base::dof::clearConstraints( fluid1.getPressure()     );
        base::dof::clearConstraints( fluid2.getVelocity()     );
        base::dof::clearConstraints( fluid2.getPressure()     );

        bubbleFeatures( immersedSurface, fluid1, surfaceQuadrature, cutQuadrature1, std::cout );
        std::cout << std::endl;
    } // end load-steps

    return 0;
}
