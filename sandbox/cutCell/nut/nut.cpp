//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   cutCell/nut/nut.cpp
//! @author Thomas Rueberg
//! @date   2014

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
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>
#include <base/cut/stabiliseBasis.hpp>
#include <base/cut/evaluateOnCutCells.hpp>

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

#include "../../generateMesh.hpp"
#include "../Laplace.hpp"
#include "../SurfaceField.hpp"
#include "../ImplicitGeometry.hpp"

namespace cutCell{
    
    const double coordTol = 1.e-6;
    
    //--------------------------------------------------------------------------
    template<typename FUN, typename DOF>
    void dirichlet( const typename FUN::arg1_type x, DOF* doFPtr, FUN fun )
    {
        typename FUN::result_type value = fun( x );
        
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

    int nut( int argc, char* argv[] );
}


//------------------------------------------------------------------------------
int cutCell::nut( int argc, char* argv[] )
{
    // spatial dimension
    const unsigned    dim = 3;

    if ( argc != 3 ) {
        std::cerr << "Usage: " << argv[0] << " N  input.dat\n"
                  << "(Compiled for dim=" << dim << ")\n\n";
        return -1;
    }

    // read name of input file
    const unsigned    numElements = boost::lexical_cast<unsigned>(    argv[1] );
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double xmax, cutThreshold, sX, sY, sZ;
    base::io::List<std::string> surfaceMeshes;
    std::string meshFile;
    unsigned stabilise; // 0 - not, 1 - Höllig
    bool compute;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "xmax",         xmax );
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
    const unsigned    kernelDegEstimate = 5;
    
    // Bulk mesh
    typedef base::Unstructured<shape,geomDeg>  Mesh;
    typedef Mesh::Node::VecDim VecDim;
    Mesh mesh;
    std::string baseName;
    if ( numElements > 0 ) {
        base::Vector<dim,unsigned>::Type N;
        VecDim a, b;
        for ( unsigned d = 0; d < dim; d++ ) {
            N[d] = numElements;
            a[d] = 0.;
            b[d] = xmax;
        }
        generateMesh<dim>( mesh, N, a, b );
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

    cutCell::ImplicitGeometry<Mesh> geometry( mesh, boundaryMesh );
    
    //--------------------------------------------------------------------------
    // go through immersed surfaces
    for ( unsigned s = 0; s < numSurfaces; s++ ) {
    
        // surface mesh from input
        const std::string surfMeshFileName = surfaceMeshes[s];

        // Surface mesh
        typedef base::Unstructured<surfShape,1,dim>    SurfMesh;
        SurfMesh surfMesh;
        {
            std::ifstream smf( surfMeshFileName.c_str() );
            base::io::smf::readMesh( smf, surfMesh );
            smf.close();
        }

        geometry.intersectSurfaceMesh( surfMesh, cutThreshold );

    }
    // Cell structures
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells = geometry.getCells();

    typedef base::cut::Cell<surfShape> SurfCell;
    std::vector<SurfCell> surfCells = geometry.getSurfCells();

    // intersection of all level sets
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSetIntersection = geometry.getLevelSet();

    typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
    SurfaceMesh surfaceMesh;
    base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, surfaceMesh );

    //--------------------------------------------------------------------------
    // FE
    typedef Laplace<Mesh,fieldDeg> Laplace;
    Laplace laplace( mesh, 1.0 );

    typedef SurfaceField<SurfaceMesh,Laplace::Field> SurfaceField;
    SurfaceField surfaceField( surfaceMesh,   laplace.getField() );
    SurfaceField boundaryField( boundaryMesh, laplace.getField() );

    
    //--------------------------------------------------------------------------
    // Quadratures
    typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
    CutQuadrature cutQuadrature( cells, true );

    typedef base::Quadrature<kernelDegEstimate,surfShape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    typedef base::cut::Quadrature<kernelDegEstimate,surfShape> SurfaceCutQuadrature;
    SurfaceCutQuadrature surfaceCutQuadrature( surfCells, true );

    // compute supports, scale basis
    const std::size_t numDoFs = std::distance( laplace.getField().doFsBegin(),
                                               laplace.getField().doFsEnd() );
    std::vector<double> supports;
    supports.resize(  numDoFs );
    
    base::cut::supportComputation( mesh, laplace.getField(), cutQuadrature,  supports );
    std::vector<std::pair<std::size_t,VecDim> > doFLocation;
    base::dof::associateLocation( laplace.getField(), doFLocation );

    typedef base::auxi::FundSolLaplace<dim> FSol;
    FSol fSol;

    typedef boost::function< Laplace::VecDoF( const VecDim& ) > FFun;
    FFun fFun = boost::bind( &FSol::fun, &fSol, _1, sourcePoint );

    // apply dirichlet constraints
    base::dof::constrainBoundary<Laplace::FEBasis>( meshBoundary.begin(),
                                                    meshBoundary.end(),
                                                    mesh, laplace.getField(),
                                                    boost::bind( &dirichlet<FFun,
                                                                 Laplace::DoF>,
                                                                 _1, _2, fFun ) );

    base::cut::stabiliseBasis( mesh, laplace.getField(), supports, doFLocation );

    // number DoFs
    const std::size_t activeDoFsU = 
        base::dof::numberDoFsConsecutively( laplace.getField().doFsBegin(),
                                            laplace.getField().doFsEnd() );

    if ( compute ) {

        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( activeDoFsU );

        typedef boost::function< Laplace::VecDoF( const VecDim&,
                                                  const VecDim&) > FFun2;
        FFun2 fFun2 = boost::bind( &FSol::coNormal, &fSol, _1, sourcePoint, _2 );

        // Neumann boundary condition -- box boundary
        boundaryField.applyNeumannBoundaryConditions( surfaceCutQuadrature,
                                                      solver, fFun2 );
        
        // Neumann boundary condition -- immersed surface
        surfaceField.applyNeumannBoundaryConditions( surfaceQuadrature,
                                                     solver, fFun2 );

        laplace.assembleBulk( cutQuadrature, solver );

        // Finalise assembly
        solver.finishAssembly();

#ifdef WRITEMAT
        {
            std::string matName = baseName + ".mat";
            std::ofstream mat( matName.c_str() );
            solver.debugLHS( mat );
        }
#endif
        // Solve
        //solver.superLUSolve();
        solver.cgSolve();

        
        // distribute results back to dofs
        base::dof::addToDoFsFromSolver( solver, laplace.getField() );

    }

    laplace.writeVTKFile(    baseName, levelSetIntersection );
    laplace.writeVTKFileCut( baseName, levelSetIntersection, cells );
    

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

    double hmin;
    {
        std::vector<double> h;
        Mesh::ElementPtrConstIter eIter = mesh.elementsBegin();
        Mesh::ElementPtrConstIter eEnd  = mesh.elementsEnd();
        for ( ; eIter != eEnd; ++eIter )  {
            h.push_back( base::mesh::elementSize( *eIter ) );
        }
     
        hmin = *std::min_element( h.begin(), h.end() );
    }

    std::cout << hmin << "  ";
    
    //--------------------------------------------------------------------------
    boost::function<FSol::Grad( const VecDim& )> derivative
        = boost::bind( &FSol::grad, &fSol, _1, sourcePoint );
    const std::pair<double,double> error =
        laplace.computeErrors( cutQuadrature,
                               fFun, derivative );
    std::cout << error.first << "  " << error.second << "  ";


    //--------------------------------------------------------------------------
    // Compute mesh volume
    const double volume = laplace.computeVolume( cutQuadrature );
    std::cout << volume << "  ";

    //--------------------------------------------------------------------------
    // Compute surface area
    double area = surfaceField.computeArea( surfaceQuadrature );
    area += boundaryField.computeArea(      surfaceCutQuadrature );

    std::cout << area << "  ";

    //std::cout << numCGIter;
    std::cout << std::endl;

    return 0;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    return cutCell::nut( argc, argv );
}
