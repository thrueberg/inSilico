//#define WRITEMAT


#include <string>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <base/dof/location.hpp>

#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/Field.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/mesh/Size.hpp>

#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>
#include <base/io/List.hpp>

#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/ScaledField.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>
#include <base/cut/stabiliseBasis.hpp>
#include <base/cut/evaluateOnCutCells.hpp>
#include <base/cut/bruteForce.hpp>

#include <base/fe/Basis.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/location.hpp>
#include <base/dof/constrainBoundary.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <heat/Laplace.hpp>

#include <solid/HyperElastic.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include <base/auxi/FundamentalSolution.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>
#include <base/post/ErrorNorm.hpp>

#include "generateMesh.hpp"

const double coordTol = 1.e-6;

//------------------------------------------------------------------------------
template<typename FUN, typename DOF>
void dirichletBC( const typename FUN::arg1_type x, DOF* doFPtr, FUN fun,
                  const bool scale, const std::vector<double>& supportAreas )
{
    typename FUN::result_type value = fun( x );

    if ( scale ) {
        const double support = supportAreas[ doFPtr -> getID() ];
        value *= support;
    }

    // fix the box face with x_1 = 0
    if ( base::auxi::almostEqualNumbers( x[0], 0. ) ) {

        for ( unsigned d = 0; d < DOF::size; d++ ) {
            if ( doFPtr -> isActive(d) ) {
                doFPtr -> constrainValue( d, value[d] );
            }
        }
    }
    return;
}

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
void writeVTKFile( const std::string& baseName,
                   const MESH&    mesh,
                   const FIELD&   field,
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet )
{
    const std::string vtkFile = baseName + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );

    std::vector<double> distances;
    std::transform( levelSet.begin(), levelSet.end(),
                    std::back_inserter( distances ),
                    boost::bind( &base::cut::LevelSet<MESH::Node::dim>::getSignedDistance, _1 ) );

    vtkWriter.writeUnstructuredGrid( mesh );
    vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
    base::io::vtk::writePointData( vtkWriter, mesh, field, "u   " );

    vtk.close();
}



//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // spatial dimension
    const unsigned    dim = 3;

    const bool isLaplace = true;
    
    if ( argc != 3 ) {
        std::cerr << "Usage: " << argv[0] << " N  input.dat\n"
                  << "(Compiled for dim=" << dim << ")\n\n";
        return -1;
    }

    // read name of input file
    const unsigned    numElements = boost::lexical_cast<unsigned>(    argv[1] );
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double xmax, E, nu, f, tolerance, cutThreshold, sX, sY, sZ;
    base::io::List<std::string> surfaceMeshes;
    unsigned numLoadSteps, maxIter;
    std::string meshFile;
    unsigned stabilise; // 0 - not, 1 - Höllig, 2 - scaling
    bool compute;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "xmax",         xmax );
        prop.registerPropertiesVar( "E",            E );
        prop.registerPropertiesVar( "nu",           nu );
        prop.registerPropertiesVar( "f",            f );
        prop.registerPropertiesVar( "numLoadSteps", numLoadSteps );
        prop.registerPropertiesVar( "maxIter",      maxIter );
        prop.registerPropertiesVar( "tolerance",    tolerance );

        prop.registerPropertiesVar( "surfaceMeshes", surfaceMeshes );

        prop.registerPropertiesVar( "sX",           sX );
        prop.registerPropertiesVar( "sY",           sY );
        prop.registerPropertiesVar( "sZ",           sZ );

        prop.registerPropertiesVar( "stabilise",    stabilise );
        prop.registerPropertiesVar( "meshFile",     meshFile );
        prop.registerPropertiesVar( "compute",      compute );
        prop.registerPropertiesVar( "cutThreshold", cutThreshold );
        
        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );
    }

    // in case of no computation, do not stabilise
    if ( not compute ) stabilise = 0;
    
    // number of immersed surfaces
    const unsigned numSurfaces = static_cast<unsigned>( surfaceMeshes.size() );
    
    // source point of the fundamental solution
    base::Vector<dim>::Type sourcePoint
        = base::constantVector<dim>( 0. );
    sourcePoint[0] = sX;
    if ( dim > 1 ) sourcePoint[1] = sY;
    if ( dim > 2 ) sourcePoint[2] = sZ;

    
    // basic attributes of the computation
    const unsigned             geomDeg  = 1;
    const unsigned             fieldDeg = 1;
    const base::Shape             shape = //base::SimplexShape<dim>::value;
        base::HyperCubeShape<dim>::value;
    const base::Shape surfShape = base::SimplexShape<dim-1>::value;

    const unsigned    kernelDegEstimate = 3;
    const unsigned              doFSize = (isLaplace ? 1 : dim);
    const unsigned              nHist   = 1;
    const bool        isSigned  = true;

    // Bulk mesh
    typedef base::Unstructured<shape,geomDeg>  Mesh;
    Mesh mesh;
    std::string baseName;
    if ( numElements > 0 ) {
        const double left = 0;
        generateMesh( mesh, numElements, left, xmax );

        baseName = "nut." + base::io::leadingZeros( numElements );
    }
    else{
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );

        baseName = base::io::baseName( meshFile, ".smf" );
    }

    // Boundary mesh
    typedef base::mesh::BoundaryMeshBinder<Mesh,true>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );
    base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                      meshBoundary.end(),
                                      mesh, boundaryMesh );
        
    typedef Mesh::Node::VecDim VecDim;

    // Cell structures
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells;

    typedef base::cut::Cell<surfShape> SurfCell;
    std::vector<SurfCell> surfCells;

    // union of all level sets
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSetUnion;

    // Surface meshes
    typedef base::Unstructured<surfShape,1,dim>    SurfMesh;

    //--------------------------------------------------------------------------
    // go through immersed surfaces
    std::cout << numSurfaces << "  ";
    for ( unsigned s = 0; s < numSurfaces; s++ ) {
    
        // surface mesh from input
        const std::string surfMeshFileName = surfaceMeshes[s];

        SurfMesh surfMesh;
        {
            std::ifstream smf( surfMeshFileName.c_str() );
            base::io::smf::readMesh( smf, surfMesh );
            smf.close();
        }

        //----------------------------------------------------------------------
        // Compute the level set data
        std::vector<LevelSet> levelSet;
        base::cut::bruteForce( mesh, surfMesh, isSigned, levelSet );

        //----------------------------------------------------------------------
        base::cut::generateCutCells( mesh, levelSet, cells, s > 0 );

        // Make cut cell structure (surface)
        base::cut::generateCutCells( boundaryMesh, levelSet, surfCells, s>0 );

        // merge level set
        if ( s == 0 ) levelSetUnion = levelSet;
        else 
            std::transform( levelSetUnion.begin(), levelSetUnion.end(),
                            levelSet.begin(),
                            levelSetUnion.begin(),
                            boost::bind( base::cut::setUnion<dim>, _1, _2 ) );

        std::for_each( cells.begin(), cells.end(),
                       boost::bind( &Cell::compress, _1, cutThreshold ) );

    }


    typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
    SurfaceMesh surfaceMesh;
    base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, surfaceMesh );

    //--------------------------------------------------------------------------
    // FE
    typedef base::fe::Basis<shape,fieldDeg>               FEBasis;
    typedef base::cut::ScaledField<FEBasis,doFSize,nHist> Field;
    Field field;
    base::dof::generate<FEBasis>( mesh, field );

    // for domain field
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type UU;

    // for surface field
    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder boundaryFieldBinder( boundaryMesh, field );
    SurfaceFieldBinder  surfaceFieldBinder(  surfaceMesh, field );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SU;

    //--------------------------------------------------------------------------
    // Quadratures
    typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
    CutQuadrature cutQuadrature( cells, true );

    typedef base::Quadrature<kernelDegEstimate,surfShape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    typedef base::cut::Quadrature<kernelDegEstimate,surfShape> SurfaceCutQuadrature;
    SurfaceCutQuadrature surfaceCutQuadrature( surfCells, true );

    // compute supports, scale basis
    const std::size_t numDoFs = std::distance( field.doFsBegin(),
                                               field.doFsEnd() );
    std::vector<double> supports;
    supports.resize(  numDoFs );
    
    base::cut::supportComputation( mesh, field, cutQuadrature,  supports );

    if ( stabilise==1 ) {
        std::vector<std::pair<std::size_t,VecDim> > doFLocation;
        base::dof::associateLocation( field, doFLocation );
        base::cut::stabiliseBasis( mesh, field, supports, doFLocation );
    }
    else if ( stabilise==2 ) {
        field.scaleAndTagBasis(  supports,  1.e-10 );
    }
    else field.tagBasis( supports,  1.e-10 );

    typedef base::Vector<1>::Type VecDoF;

    typedef base::auxi::FundSolLaplace<dim> FSol;
    FSol fSol;

    typedef boost::function< VecDoF( const VecDim& ) > FFun;
    FFun fFun =
        boost::bind( &FSol::fun, &fSol, _1, sourcePoint );
    

    // apply dirichlet constraints
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, field, 
                                           boost::bind( &dirichletBC<FFun,
                                                        Field::DegreeOfFreedom>,
                                                        _1, _2, fFun,
                                                        (stabilise==2),
                                                        boost::ref(supports) ) );
    // number DoFs
    const std::size_t activeDoFsU = 
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );

    typedef mat::hypel::NeoHookeanCompressible Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );
    typedef solid::HyperElastic<Material,UU::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    typedef heat::Laplace<UU::Tuple> Laplace;
    Laplace laplace( 1. );

    //const unsigned effIter = (isLaplace ? 1 : maxIter );

    unsigned numCGIter;
    if ( compute ) {
    //for ( unsigned iter = 0; iter < effIter; iter++ ) {

        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( activeDoFsU );

        typedef boost::function< VecDoF( const VecDim&,
                                         const VecDim&) > FFun2;
        FFun2 fFun2 = boost::bind( &FSol::coNormal, &fSol, _1, sourcePoint, _2 );

        // Neumann boundary condition -- box boundary
        base::asmb::neumannForceComputation<SU>( //surfaceQuadrature,
                                                 surfaceCutQuadrature,
                                                 solver, boundaryFieldBinder,
                                                 fFun2 );

        // Neumann boundary condition -- immersed surface
        base::asmb::neumannForceComputation<SU>( surfaceQuadrature,
                                                 solver, surfaceFieldBinder,
                                                 fFun2 );


        if ( isLaplace ) {
            base::asmb::stiffnessMatrixComputation<UU>( cutQuadrature, solver,
                                                        fieldBinder, laplace, 
                                                        false ); //iter > 0 );
        }
        else {

            base::asmb::computeResidualForces<UU>( cutQuadrature, solver,
                                                   fieldBinder, hyperElastic );
            
            base::asmb::stiffnessMatrixComputation<UU>( cutQuadrature, solver,
                                                        fieldBinder, hyperElastic,
                                                        false ); //iter > 0 );

        }

        // Finalise assembly
        solver.finishAssembly();

#ifdef WRITEMAT
        {
            std::string matName = baseName + ".mat";
            std::ofstream mat( matName.c_str() );
            solver.debugLHS( mat );

            std::string vecName = baseName + ".vec";
            std::ofstream vec( vecName.c_str() );
            solver.debugRHS( vec );
        }
#endif
        
        // norm of residual
        //const double conv1 = solver.norm() / E;
        //std::cout << "* " << iter << " " << conv1 << " ";

        //if ( conv1 < tolerance ) {
        //    std::cout << std::endl;
        //    break;
        //}

        // Solve
        //solver.superLUSolve();
        numCGIter = solver.cgSolve();

        
        // distribute results back to dofs
        base::dof::addToDoFsFromSolver( solver, field );

        //const double conv2 = solver.norm();
        //std::cout << conv2 << std::endl;
        //if ( conv2 < tolerance ) break;

    }

    
    {
        std::vector<double> cutCellDistance;
        base::cut::evaluateCutCellNodeDistances( mesh, levelSetUnion, cells,
                                                 cutCellDistance );

        std::vector<VecDoF> cutCellSolution;
        base::cut::evaluateFieldAtCutCellNodes( mesh, field, cells,
                                                cutCellSolution );

        typedef base::cut::VolumeSimplexMesh<Mesh> VolumeSimplexMesh;
        VolumeSimplexMesh volumeSimplexMeshIn;
        base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMeshIn,
                                                  true );

        const std::string vtkFileIn = baseName + ".volin.vtk";
        std::ofstream vtkIn( vtkFileIn.c_str() );
        base::io::vtk::LegacyWriter vtkWriterIn( vtkIn );
        vtkWriterIn.writeUnstructuredGrid( volumeSimplexMeshIn );
        vtkWriterIn.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
        vtkWriterIn.writePointData( cutCellSolution.begin(), cutCellSolution.end(), "u" );

        VolumeSimplexMesh volumeSimplexMeshOut;
        base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMeshOut,
                                                  false );
        
        const std::string vtkFileOut = baseName + ".volout.vtk";
        std::ofstream vtkOut( vtkFileOut.c_str() );
        base::io::vtk::LegacyWriter vtkWriterOut( vtkOut );
        vtkWriterOut.writeUnstructuredGrid( volumeSimplexMeshOut );
        vtkWriterOut.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
        vtkWriterOut.writePointData( cutCellSolution.begin(), cutCellSolution.end(), "u" );

        typedef base::cut::SurfaceSimplexMesh<Mesh> SurfaceSimplexMesh;
        SurfaceSimplexMesh surfaceSimplexMesh;
        base::cut::extractSurfaceMeshFromCutCells( mesh, cells, surfaceSimplexMesh );

        const std::string vtkFileSurf = baseName + ".surf.vtk";
        std::ofstream vtkSurf( vtkFileSurf.c_str() );
        base::io::vtk::LegacyWriter vtkWriterSurf( vtkSurf );

        vtkWriterSurf.writeUnstructuredGrid( surfaceSimplexMesh );
        vtkWriterSurf.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
        vtkWriterSurf.writePointData( cutCellSolution.begin(), cutCellSolution.end(), "u" );

        //const std::string vtkFileBound = baseName + ".bound.vtk";
        //std::ofstream vtkBound( vtkFileBound.c_str() );
        //base::io::vtk::LegacyWriter vtkWriterBound( vtkBound );
        //
        //std::vector<double> distances;
        //std::transform( levelSetUnion.begin(), levelSetUnion.end(),
        //                std::back_inserter( distances ),
        //                boost::bind( &LevelSet::getSignedDistance, _1 ) );
        //
        //vtkWriterBound.writeUnstructuredGrid( boundaryMesh );
        //vtkWriterBound.writePointData( distances.begin(), distances.end(), "distances" );

    }

    //--------------------------------------------------------------------------
    {
        typedef base::cut::VolumeSimplexMesh<BoundaryMesh> BoundarySimplexMesh;
        BoundarySimplexMesh boundarySimplexMeshIn;
        base::cut::extractVolumeMeshFromCutCells( boundaryMesh,
                                                  surfCells,
                                                  boundarySimplexMeshIn, true );
        
        const std::string smfBoundIn = baseName + ".boundIn.smf";
        std::ofstream smfIn( smfBoundIn.c_str() );
        base::io::smf::writeMesh( boundarySimplexMeshIn, smfIn );
        smfIn.close();

        BoundarySimplexMesh boundarySimplexMeshOut;
        base::cut::extractVolumeMeshFromCutCells( boundaryMesh,
                                                  surfCells,
                                                  boundarySimplexMeshOut, false );
        
        const std::string smfBoundOut = baseName + ".boundOut.smf";
        std::ofstream smfOut( smfBoundOut.c_str() );
        base::io::smf::writeMesh( boundarySimplexMeshOut, smfOut );
        smfOut.close();
    }

    // write a vtk file
    writeVTKFile( baseName, mesh, field, levelSetUnion );

    double hmin, hmax;
    {
        std::vector<double> h;
        Mesh::ElementPtrConstIter eIter = mesh.elementsBegin();
        Mesh::ElementPtrConstIter eEnd  = mesh.elementsEnd();
        for ( ; eIter != eEnd; ++eIter )  {
            h.push_back( base::mesh::elementSize( *eIter ) );
        }
     
        hmin = *std::min_element( h.begin(), h.end() );
        hmax = *std::max_element( h.begin(), h.end() );
    }

    std::cout << hmin << "  " << hmax << "  ";
    //std::cout << numElements << "  ";
    
    //--------------------------------------------------------------------------
    // compute L2-error
    std::cout << base::post::errorComputation<0>(
                  cutQuadrature, mesh, field,
                  fFun )
              << "  ";

    //--------------------------------------------------------------------------
    // compute H1-error
    boost::function<FSol::Grad( const VecDim& )> derivative
        = boost::bind( &FSol::grad, &fSol, _1, sourcePoint );
    
    std::cout << base::post::errorComputation<1>(
                  cutQuadrature, mesh, field,
                  derivative )
              << "  ";
    
#if 1
    //--------------------------------------------------------------------------
    // Compute mesh volume
    double volume = 0.;
    base::asmb::simplyIntegrate<UU>( cutQuadrature, volume, fieldBinder,
                                     base::kernel::Measure<UU::Tuple>() );
    std::cout << volume << "  ";

    //--------------------------------------------------------------------------
    // Compute surface area
    double area = 0.;
    base::asmb::simplyIntegrate<SU>( surfaceQuadrature, area, surfaceFieldBinder,
                                     base::kernel::Measure<SU::Tuple>() );

    base::asmb::simplyIntegrate<SU>( surfaceCutQuadrature, area, boundaryFieldBinder,
                                     base::kernel::Measure<SU::Tuple>() );


    std::cout << area << "  ";
    
#endif

    //std::cout << numCGIter;
    std::cout << std::endl;

    return 0;
}
