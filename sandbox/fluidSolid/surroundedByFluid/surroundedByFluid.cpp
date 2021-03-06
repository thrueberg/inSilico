//#define MEMCNT
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
#include <base/cut/evaluateOnCutCells.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>

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
#include "../generateMesh.hpp"

#include "InputFSI.hpp"

#ifdef MEMCNT
#include <base/auxi/Memory.hpp>
#endif

//==============================================================================
//! Write data to VTK
template<typename MESH, typename DISP,
         typename VELOC, typename PRESS, typename MATERIAL, typename CELL>
void writeVTKFile( const std::string& baseName,
                   const unsigned step,
                   const MESH&    mesh,
                   const DISP&    displacement,
                   const VELOC&   velocity,
                   const PRESS&   pressure,
                   const MATERIAL& material, 
                   const double    viscosity,
                   const std::vector<CELL>& cells, 
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet )
{
    VTKWriter<MESH> vtkWriter( mesh, baseName, step );

    vtkWriter.writeDistances( levelSet );
    vtkWriter.writeField( displacement,         "disp" );
    vtkWriter.writeField( velocity,             "velocity" );
    vtkWriter.writeField( pressure,             "pressure" );
    vtkWriter.writePreviousField( displacement, "prevDisp" );

    // write fluid stresses
    {
        typename base::Vector<MESH::Node::dim>::Type xi =
            base::ShapeCentroid<MESH::Element::shape>::apply();
        
        std::vector<typename base::Matrix<MESH::Node::dim,MESH::Node::dim>::Type> sigma;

        typedef base::asmb::FieldBinder<const MESH,const VELOC,const PRESS> Field;
        Field field( mesh, velocity, pressure );

        typedef typename Field::template TupleBinder<1,2>::Type UP;
        fluid::Stress<typename UP::Tuple> stress( viscosity );

        typename Field::FieldIterator first = field.elementsBegin();
        typename Field::FieldIterator  last = field.elementsEnd();
        for ( std::size_t c = 0; first != last; ++first, c++ ) {

            if ( cells[c].isOutside() )
                sigma.push_back( stress( UP::makeTuple( *first ), xi ) );
            else 
                sigma.push_back( base::constantMatrix<MESH::Node::dim,
                                 MESH::Node::dim>( 0. ) );
        }
        

        vtkWriter.writeCellData( sigma.begin(), sigma.end(), "sigmaF" );

    }

    // write solid stresses
    {
        typename base::Vector<MESH::Node::dim>::Type xi =
            base::ShapeCentroid<MESH::Element::shape>::apply();
        
        std::vector<typename base::Matrix<MESH::Node::dim,MESH::Node::dim>::Type> sigma;

        typedef base::asmb::FieldBinder<const MESH,const DISP> Field;
        Field field( mesh, displacement );

        typedef typename Field::template TupleBinder<1,1>::Type UU;
        CauchyStressSolid<typename UU::Tuple,MATERIAL> stress( material );

        typename Field::FieldIterator first = field.elementsBegin();
        typename Field::FieldIterator  last = field.elementsEnd();
        for ( std::size_t c = 0; first != last; ++first, c++ ) {

            if ( cells[c].isInside() )
                sigma.push_back( stress( UU::makeTuple( *first ), xi ) );
            else 
                sigma.push_back( base::constantMatrix<MESH::Node::dim,
                                 MESH::Node::dim>( 0. ) );
        }
        
        vtkWriter.writeCellData( sigma.begin(), sigma.end(), "sigmaS" );


    }

#ifdef VTKVERBOSE
    vtkWriter.writeDoFStatus( displacement, "statusD" );
    vtkWriter.writeDoFStatus( velocity,     "statusU" );
    vtkWriter.writeDoFStatus( pressure,     "statusP" );
#endif
}

//------------------------------------------------------------------------------
template<typename MESH, typename DISP,
         typename VELOC, typename PRESS, typename MATERIAL, typename CELL>

void writeVTKFileCut( const std::string& baseName,
                      const unsigned step,
                      const MESH&    mesh,
                      const DISP&    displacement,
                      const VELOC&   velocity,
                      const PRESS&   pressure,
                      const MATERIAL& material, 
                      const double    viscosity,
                      const std::vector<CELL>&     cells,
                      const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet )
{
    typedef typename base::Vector<DISP::DegreeOfFreedom::size>::Type VecDoF;
    
    // write cut-volume files
    std::vector<double> cutCellDistance;
    base::cut::evaluateCutCellNodeDistances( mesh, levelSet, cells,
                                             cutCellDistance );
            
    std::vector<VecDoF> cutCellDisplacement;
    base::cut::evaluateFieldAtCutCellNodes( mesh, displacement, cells,
                                            cutCellDisplacement );
    std::vector<VecDoF> cutCellVelocity;
    base::cut::evaluateFieldAtCutCellNodes( mesh, velocity, cells,
                                            cutCellVelocity );

    typedef base::cut::VolumeSimplexMesh<MESH> VolumeSimplexMesh;
    VolumeSimplexMesh volumeSimplexMeshSolid, volumeSimplexMeshFluid;
    base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMeshSolid, true  );
    base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMeshFluid, false );

    const std::string vtkFileS =
        baseName + "CutS." + base::io::leadingZeros( step ) + ".vtk";
    const std::string vtkFileF =
        baseName + "CutF." + base::io::leadingZeros( step ) + ".vtk";
        
    std::ofstream vtkS( vtkFileS.c_str() );
    base::io::vtk::LegacyWriter vtkWriterS( vtkS );
    vtkWriterS.writeUnstructuredGrid( volumeSimplexMeshSolid );
    vtkWriterS.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
    vtkWriterS.writePointData( cutCellDisplacement.begin(), cutCellDisplacement.end(), "disp" );

    // write solid stresses
    {
        typedef base::asmb::FieldBinder<const MESH, const DISP> FieldBinder;
        FieldBinder fieldBinder( mesh, displacement );
        typedef typename FieldBinder::template TupleBinder<1,1>::Type UU;

        CauchyStressSolid<typename UU::Tuple,MATERIAL> cauchy( material );
        
        std::vector<typename base::Matrix<MESH::Node::dim,MESH::Node::dim>::Type> stresses;
        
        base::cut::evaluateInCutCells<UU>(
            fieldBinder, cauchy, cells, stresses, true );

        vtkWriterS.writeCellData( stresses.begin(), stresses.end(), "sigmaS" );
    }

    std::ofstream vtkF( vtkFileF.c_str() );
    base::io::vtk::LegacyWriter vtkWriterF( vtkF );
    vtkWriterF.writeUnstructuredGrid( volumeSimplexMeshFluid );
    vtkWriterF.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
    vtkWriterF.writePointData( cutCellVelocity.begin(), cutCellVelocity.end(), "velocity" );

    // write fluid stresses
    {
        typedef base::asmb::FieldBinder<const MESH, const VELOC, const PRESS> FieldBinder;
        FieldBinder fieldBinder( mesh, velocity, pressure );
        typedef typename FieldBinder::template TupleBinder<1,2>::Type UP;

        fluid::Stress<typename UP::Tuple> cauchy( viscosity );
        
        std::vector<typename base::Matrix<MESH::Node::dim,MESH::Node::dim>::Type> stresses;
        
        base::cut::evaluateInCutCells<UP>(
            fieldBinder, cauchy, cells, stresses, false );

        vtkWriterF.writeCellData( stresses.begin(), stresses.end(), "sigmaF" );
    }


}
//------------------------------------------------------------------------------
//! Write the VTK file of surface mesh
template<typename SMESH,  typename DISP, typename VELOC, typename PRESS, typename MAT>
void writeSurfaceVTKFile( const std::string& name,
                          const unsigned num, 
                          const SMESH& immersedSurface,
                          const DISP&  displacement,
                          const MAT&   material,
                          const VELOC& velocity,
                          const PRESS& pressure,
                          const double viscosity )
{
    VTKWriter<SMESH> vtkWriter( immersedSurface, name, num );
    
    typedef typename SMESH::Node::VecDim VecDim;

    //--------------------------------------------------------------------------
    // Evaluate displacement, velocity and pressure
    {
        // storage
        std::vector<VecDim> nodalD, nodalU;
        std::vector<double> nodalP;
        // go through all surface elements
        typename SMESH::ElementPtrConstIter eBegin = immersedSurface.elementsBegin();
        typename SMESH::ElementPtrConstIter eEnd   = immersedSurface.elementsEnd();
        for ( ; eBegin != eEnd; ++eBegin ) {
        
            typename SMESH::Element::DomainElement* geomEp =
                (*eBegin) -> getDomainElementPointer();
        
            const std::size_t elemID = geomEp -> getID();

            // get the field elements
            typename DISP::Element*  dispEp  = displacement.elementPtr( elemID );
            typename VELOC::Element* velocEp = velocity.elementPtr(     elemID );
            typename PRESS::Element* pressEp = pressure.elementPtr(     elemID );

            // go through the local coordinates of the surface element's nodes
            typename SMESH::Element::ParamIter   pIter = (*eBegin) -> parametricBegin();
            typename SMESH::Element::ParamIter   pEnd  = (*eBegin) -> parametricEnd();
            for ( ; pIter != pEnd; ++pIter ) {
                // evaluate and store
                nodalD.push_back( base::post::evaluateField( geomEp, dispEp,  *pIter ) );
                nodalU.push_back( base::post::evaluateField( geomEp, velocEp, *pIter ) );
                nodalP.push_back( base::post::evaluateField( geomEp, pressEp, *pIter )[0] );
            }
        }
        
        // pass on to writer
        vtkWriter.writePointData( nodalD.begin(), nodalD.end(), "Disp" );
        vtkWriter.writePointData( nodalU.begin(), nodalU.end(), "Veloc" );
        vtkWriter.writePointData( nodalP.begin(), nodalP.end(), "Press" );
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
        std::vector<typename SMESH::Node::VecDim> tractionsF, tractionsS;

        // bind surface to fluid fields
        typedef base::asmb::SurfaceFieldBinder<const SMESH,const VELOC, const PRESS> SFB1;
        SFB1 sfb1( immersedSurface, velocity, pressure );

        // bind surface to solid displacement
        typedef base::asmb::SurfaceFieldBinder<const SMESH,const DISP> SFB2;
        SFB2 sfb2( immersedSurface, displacement );

        // types of domain field tuples
        typedef typename SFB1::template TupleBinder<1,2>::Type STBUP;
        typedef typename base::asmb::DomainFieldElementPointerTuple<typename STBUP::Tuple>::Type
            DFTUP;
        
        typedef typename SFB2::template TupleBinder<1>::Type STBD;
        typedef typename base::asmb::DomainFieldElementPointerTuple<typename STBD::Tuple>::Type
            DFTD;

        // traction objects
        fluid::Traction<DFTUP> tractionF( viscosity );
        ::Traction<DFTD,MAT>   tractionS( material );

        // collect the sum of the forces (unused)
        VecDim sumF = base::constantVector<SMESH::Node::dim>( 0. );
        VecDim sumS = base::constantVector<SMESH::Node::dim>( 0. );

        // go through all surface element tuples
        typename SFB1::FieldIterator sEBegin1 = sfb1.elementsBegin();
        typename SFB2::FieldIterator sEBegin2 = sfb2.elementsBegin();
        typename SFB1::FieldIterator sEEnd    = sfb1.elementsEnd();
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
            const VecDim tf = tractionF( dftf, xi, normals[ ctr ] );
            tractionsF.push_back( tf );

            // solid domain element tuple
            const DFTD dftd =
                base::asmb::DomainFieldElementPointerTuple<
                    typename STBD::Tuple>::convert( *sEBegin2 );

            // compute solid side traction and store
            const VecDim ts = tractionS( dftd, xi, normals[ ctr ] );
            tractionsS.push_back( ts );

            // add to sum (like a mid-point rule integration)
            sumF += areas[ctr] * tf;
            sumS += areas[ctr] * ts;
        }
    
        vtkWriter.writeCellData( tractionsF.begin(), tractionsF.end(), "tractionF" );
        vtkWriter.writeCellData( tractionsS.begin(), tractionsS.end(), "tractionS" );
    }
    
    return;
}

//==============================================================================
//------------------------------------------------------------------------------
// Compute the enclosed volume, the centroid and principal directions
template<typename SURFMESH, typename BOUNDMESH, typename SQUAD, typename SCQUAD>
void surfaceFeatures( const SURFMESH&   surfaceMesh,
                      const BOUNDMESH& boundaryMesh,
                      const SQUAD&  surfaceQuadrature,
                      const SCQUAD& surfaceCutQuadrature,
                      std::ostream& out )
{
    static const unsigned dim = SURFMESH::Node::dim;
    typedef typename base::Vector<dim>::Type     VecDim;
    typedef typename base::Matrix<dim,dim>::Type MatDimDim;
 
    // get shape features
    const double volume =
        surf::enclosedVolume(  surfaceMesh, surfaceQuadrature    ) +
        surf::enclosedVolume( boundaryMesh, surfaceCutQuadrature );  

    const VecDim moment =
        surf::enclosedVolumeMoment(  surfaceMesh, surfaceQuadrature    ) + 
        surf::enclosedVolumeMoment( boundaryMesh, surfaceCutQuadrature );  

    const MatDimDim secondMoment =
        surf::enclosedVolumeSecondMoment(  surfaceMesh, surfaceQuadrature    ) + 
        surf::enclosedVolumeSecondMoment( boundaryMesh, surfaceCutQuadrature );  

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
     
    out << volume << "  " << centroid.transpose() << "  "
        << princVal.transpose() << "  "
        << aa.first << "  "  << (aa.second).transpose();
}

//==============================================================================
// MAIN
int main( int argc, char* argv[] )
{
    // spatial dimension
    const unsigned    dim = SPACEDIM; 

#ifdef MEMCNT
    const double initialMemory = base::auxi::memoryUsageInMegaBytes();
#endif
    typedef base::solver::Eigen3           Solver;
    
    // Check input arguments
    if ( argc != 2 ) {
        std::cerr << "Usage: " << argv[0] << " input.dat\n"
                  << "(Compiled for dim=" << dim << ")\n\n";
        return -1;
    }
    
    // read name of input file
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[1] );
    const InputFSI<dim> ui( inputFile );

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
    std::vector<LevelSet> levelSet;
    base::cut::analyticLevelSet( mesh,
                                 base::cut::Sphere<dim>( ui.R, ui.center ),
                                 true, levelSet );

    //--------------------------------------------------------------------------
    std::vector<TA::Cell> cells;
    base::cut::generateCutCells( mesh, levelSet, cells );
    std::vector<TA::SurfCell> surfCells;
    base::cut::generateCutCells( boundaryMesh, levelSet, surfCells );

    // Elastic material
    Material material( mat::Lame::lambda( ui.E, ui.nu),
                       mat::Lame::mu(     ui.E, ui.nu ) );

    // Solid and fluid handlers
    typedef Solid<TA::Mesh, Material>                        Solid;
    typedef Fluid<TA::Mesh, stokesStabil>                    Fluid;
    typedef SolidFluidInterface<TA::SurfaceMesh,Solid,Fluid> SFI;
    
    Solid solid( mesh, material );
    Fluid fluid( mesh, ui.viscosity, ui.alpha );

    // find geometry association for the dofs
    std::vector<std::pair<std::size_t,TA::VecDim> > doFLocationD, doFLocationU;
    base::dof::associateLocation( solid.getDisplacement(), doFLocationD );
    base::dof::associateLocation( fluid.getVelocity(),     doFLocationU );

    // Quadrature along a surface
    TA::CutQuadrature     cutQuadratureSolid(   cells, true  ); // in
    TA::CutQuadrature     cutQuadratureFluid(   cells, false ); // out
    TA::SurfaceQuadrature surfaceQuadrature;
    TA::SurfCutQuadrature surfaceCutQuadrature( surfCells, true );

    // register particles to be traced during the simulation
    base::post::ParticleTracer<dim> particleTracer;
    const std::string traces = "traces";
    {
        const boost::array<double,5> factors =
            {{ 0.5, 0.75, 0.9, 0.95, 0.975 }};
        
        for ( unsigned d = 0; d < dim; d++ ) {
            TA::VecDim x = base::constantVector<dim>( 0. );
            for ( unsigned n = 0; n < factors.size(); n++ ) {
                x[d] = ui.center[d] + factors[n] * ui.R;
                particleTracer.registerPoint( x );
            }
        }

        std::ofstream pt( traces.c_str() );
        particleTracer.writeLatest( pt );
    }


    //--------------------------------------------------------------------------
    // Run a zero-th step in order to write the data before computation
    // (store the previously enclosed volume)
    double prevVolume;
    {
        std::cout << 0 << "  " << 0. << "  0  " << std::flush;

        // for surface field
        TA::SurfaceMesh immersedSurface;
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cells, immersedSurface );

        writeVTKFile( "euler", 0, mesh, solid.getDisplacement(),
                      fluid.getVelocity(), fluid.getPressure(),
                      material, ui.viscosity, cells, levelSet );
        writeVTKFileCut( "euler", 0, mesh, solid.getDisplacement(),
                         fluid.getVelocity(), fluid.getPressure(),
                         material, ui.viscosity, cells, levelSet );
        
        writeSurfaceVTKFile( "eulerSurf", 0, immersedSurface,
                             solid.getDisplacement(), material, 
                             fluid.getVelocity(), fluid.getPressure(), ui.viscosity );

        surfaceFeatures( immersedSurface, boundaryMesh, 
                         surfaceQuadrature,  surfaceCutQuadrature,
                         std::cout );
        
        std::cout << std::endl;

        // store the contained volume
        prevVolume =
            surf::enclosedVolume( immersedSurface, surfaceQuadrature    ) + 
            surf::enclosedVolume( boundaryMesh,    surfaceCutQuadrature );  

    }

    //--------------------------------------------------------------------------
    // * Loop over time steps *
    //--------------------------------------------------------------------------
    const double supportThreshold = std::numeric_limits<double>::min();
    for ( unsigned step = 0; step < ui.numLoadSteps; step++ ) {
        
        // Interface handler
        TA::SurfaceMesh immersedSurface;
        base::cut::generateSurfaceMesh<TA::Mesh,TA::Cell>( mesh, cells, immersedSurface );
        SFI sfi( immersedSurface,
                 solid.getDisplacement(), fluid.getVelocity(),
                 fluid.getPressure(), ui.viscosity );

        std::cout << step+1 << "  " << (step+1) * ui.dt << "  " << std::flush;

        //----------------------------------------------------------------------
        // 0) Solve Lagrangian problem
        
        // compute supports
        std::vector<double> supportsD, supportsU, supportsP;
        base::cut::supportComputation( mesh, solid.getDisplacement(),
                                       cutQuadratureSolid, supportsD );
        base::cut::supportComputation( mesh, fluid.getVelocity(), cutQuadratureFluid, supportsU );
        base::cut::supportComputation( mesh, fluid.getPressure(), cutQuadratureFluid, supportsP );

        // Hard-wired Dirichlet constraints (fluid BC)
        {
            base::dof::constrainBoundary<Fluid::FEBasisU>(
                meshBoundary.begin(),
                meshBoundary.end(),
                mesh, fluid.getVelocity(),
                boost::bind( &dirichletBC<dim,Fluid::DoFU>,
                             _1, _2, ui.ubar, bbox, ui.bc ) );

            base::dof::constrainBoundary<Solid::FEBasisD>(
                meshBoundary.begin(),
                meshBoundary.end(),
                mesh, solid.getDisplacement(),
                boost::bind( &dirichletBC<dim,Solid::DoFD>,
                             _1, _2, ui.dt * ui.ubar, bbox, ui.bc ) );

            // Fix first pressure dof in case of closed cavity flow
            if ( ui.bc == CAVITY ) {
                Fluid::Pressure::DoFPtrIter pIter = fluid.getPressure().doFsBegin();
                (*pIter) -> constrainValue( 0, 0.0 );
            }
        }
    
        // Basis stabilisation
        base::cut::stabiliseBasis( mesh, solid.getDisplacement(), supportsD, doFLocationD );
        base::cut::stabiliseBasis( mesh, fluid.getVelocity(),     supportsU, doFLocationU );
        base::cut::stabiliseBasis( mesh, fluid.getPressure(),     supportsP, doFLocationD );
        
        // number DoFs 
        const std::size_t activeDoFsD = 
            base::dof::numberDoFsConsecutively( solid.getDisplacement().doFsBegin(),
                                                solid.getDisplacement().doFsEnd() );
        const std::size_t activeDoFsU = 
            base::dof::numberDoFsConsecutively( fluid.getVelocity().doFsBegin(),
                                                fluid.getVelocity().doFsEnd(),
                                                activeDoFsD );
        const std::size_t activeDoFsP = 
            base::dof::numberDoFsConsecutively( fluid.getPressure().doFsBegin(),
                                                fluid.getPressure().doFsEnd(),
                                                activeDoFsU + activeDoFsD );

        // start from 0 for the increments, set rest to zero too
        base::dof::clearDoFs( solid.getDisplacement() );
        base::dof::clearDoFs( fluid.getPressure() );
        base::dof::clearDoFs( fluid.getVelocity() );

        //----------------------------------------------------------------------
        // Nonlinear iterations
        unsigned iter = 0;
        for ( ; iter < ui.maxIter; iter++ ) {

            // Create a solver object
            Solver solver( activeDoFsD + activeDoFsU + activeDoFsP );

            // register to solver
            solid.registerInSolver( solver );
            fluid.registerInSolver( solver );
            sfi.registerInSolver(   solver );

            // Compute as(d,dd) and Da(d,dd)[Delta d]
            solid.assembleBulk( cutQuadratureSolid, solver, iter );

            solid.bodyForce( cutQuadratureSolid, solver,
                             boost::bind( unitForce<dim-1,dim>, _1, -0.98*ui.rhoS ) );

            // Computate af(u,p; du,dp)
            fluid.assembleBulk( cutQuadratureFluid, solver );

            fluid.bodyForce( cutQuadratureFluid, solver,
                             boost::bind( unitForce<dim-1,dim>, _1, -0.98 * ui.rhoF ) );
            
            // Penalty terms for int_Gamma (dot(d) - u) (E delta d - mu delta u) ds
            sfi.assemblePenaltyTerms( surfaceQuadrature, solver, ui.penaltyFac, ui.dt );

            // Boundary energy terms (beta = 0)  [ -int_Gamma (...) ds ]
            // Fluid - Mortar
            sfi.assembleEnergyTerms( surfaceQuadrature, solver, ui.dt );
            

#ifdef MEMCNT            
            double currentMem = base::auxi::memoryUsageInMegaBytes() - initialMemory;
            std::cout << "\n* before finishAssembly() " << currentMem << std::endl;
#endif
            
            // Finalise assembly
            solver.finishAssembly();

#ifdef MEMCNT            
            currentMem = base::auxi::memoryUsageInMegaBytes() - initialMemory;
            std::cout << "\n* after finishAssembly() " << currentMem << std::endl;
#endif

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
            base::dof::addToDoFsFromSolver( solver, solid.getDisplacement() ); //<>
            base::dof::setDoFsFromSolver(   solver, fluid.getVelocity() );
            base::dof::setDoFsFromSolver(   solver, fluid.getPressure() );

            const double conv2 = solver.norm(0, activeDoFsD);

            if ( conv2 < ui.tolerance ) break;
            
        } // finish non-linear iteration
        std::cout << iter << "  " << std::flush;

        // write a vtk file
        writeVTKFile( "euler", step + 1, 
                      mesh, solid.getDisplacement(),
                      fluid.getVelocity(), fluid.getPressure(),
                      material, ui.viscosity, cells, levelSet );
        writeVTKFileCut( "euler", step + 1,
                         mesh, solid.getDisplacement(),
                         fluid.getVelocity(), fluid.getPressure(),
                         material, ui.viscosity, cells, levelSet );
        writeSurfaceVTKFile( "eulerSurf", step + 1, immersedSurface,
                             solid.getDisplacement(), material, 
                             fluid.getVelocity(), fluid.getPressure(), ui.viscosity );

        {
            particleTracer.update( mesh, solid.getDisplacement(), ui.findTolerance, 10 );
            std::ofstream pt( traces.c_str(), std::ios_base::app );
            particleTracer.writeLatest( pt );
        }

        //----------------------------------------------------------------------
        // 1) Pass Data to complementary domain
        {
            const double extL = h; // extrapolation length

            // pass extrapolated solution to inactive DoFs
            extrapolateToFictitious( mesh, solid.getDisplacement(),
                                     extL, doFLocationD, levelSet, bbox );
        }

        //----------------------------------------------------------------------
        // 2) Geometry update
        {
            // get shape features
            const double volume =
                surf::enclosedVolume( immersedSurface, surfaceQuadrature    ) + 
                surf::enclosedVolume( boundaryMesh,    surfaceCutQuadrature );  

            const TA::VecDim moment =
                surf::enclosedVolumeMoment( immersedSurface, surfaceQuadrature    ) + 
                surf::enclosedVolumeMoment( boundaryMesh,    surfaceCutQuadrature );  
            

            const TA::VecDim centroid = moment / volume;

            const double factor = std::pow( prevVolume / volume, 
                                            1./static_cast<double>( dim ) );
     
            // Move with velocity solution
            moveSurface( mesh, fluid.getVelocity(), immersedSurface, ui.dt, factor, centroid );

            
            base::cut::bruteForce( mesh, immersedSurface, true, levelSet );

            // update the cut cell structure
            base::cut::generateCutCells( mesh,         levelSet, cells );
            base::cut::generateCutCells( boundaryMesh, levelSet, surfCells );


            // store the contained volume
            prevVolume =
                surf::enclosedVolume( immersedSurface, surfaceQuadrature    ) + 
                surf::enclosedVolume( boundaryMesh,    surfaceCutQuadrature );  
        }

        //----------------------------------------------------------------------
        // 3) advect data
        {
            base::cut::supportComputation( mesh, solid.getDisplacement(),
                                           cutQuadratureSolid, supportsD ); 
        
            // Find the location of the DoFs in the previous configuration
            std::vector<std::pair<std::size_t,TA::VecDim> > previousDoFLocation;
            findPreviousDoFLocations( mesh, solid.getDisplacement(), supportsD, 
                                      doFLocationD, previousDoFLocation, bbox, 
                                      supportThreshold, ui.findTolerance, 10 );
            
            // add previous solution to current solution: u_{n+1} = \Delta u + u_n
            {
                Solid::Displacement::DoFPtrIter dIter = solid.getDisplacement().doFsBegin();
                Solid::Displacement::DoFPtrIter dEnd  = solid.getDisplacement().doFsEnd();
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
            advectField( mesh, solid.getDisplacement(), previousDoFLocation, supportsD, 
                         supportThreshold );
        }

        //----------------------------------------------------------------------
        // push history
        base::dof::pushHistory( solid.getDisplacement() );

        // remove the linear constraints used in the stabilisation process
        base::dof::clearConstraints( solid.getDisplacement() );
        base::dof::clearConstraints( fluid.getVelocity()     );
        base::dof::clearConstraints( fluid.getPressure()     );

        surfaceFeatures( immersedSurface, boundaryMesh, 
                         surfaceQuadrature,  surfaceCutQuadrature, std::cout );

        std::cout << std::endl;

    } // end load-steps

    return 0;
}
