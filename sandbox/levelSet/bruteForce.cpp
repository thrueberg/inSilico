#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/io/Format.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/cut/LevelSet.hpp>
#include <base/cut/bruteForce.hpp>
#include <base/cut/Cell.hpp>

#include <base/cut/generateCutCells.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>
#include <base/cut/evaluateOnCutCells.hpp>

#include <base/cut/Quadrature.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/kernel/Measure.hpp>


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    // User input
    if ( argc < 3 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf surf.smf [surf2.smf] \n\n";
        return -1;
    }

    // volume mesh and basename
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    // number of surface
    const unsigned numSurfaces = static_cast<unsigned>( argc ) - 2;

    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
    const unsigned    dim       = 3;
    const bool        isSigned  = true;
    
    const base::Shape shape     = base::HyperCubeShape<dim>::value;
    const base::Shape surfShape = base::SimplexShape<dim-1>::value;

    //--------------------------------------------------------------------------
    // Domain mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;

    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Boundary mesh
    typedef base::mesh::BoundaryMeshBinder<Mesh,true>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
        
    }

    // Cell structures
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells;

    typedef base::cut::Cell<surfShape> SurfCell;
    std::vector<SurfCell> surfCells;

    // intersection of all level sets
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSetIntersection;

    //--------------------------------------------------------------------------
    // Surface meshes
    typedef base::Unstructured<surfShape,1,dim>    SurfMesh;


    for ( unsigned s = 0; s < numSurfaces; s++ ) {

        // surface mesh from input
        const std::string surfMeshFileName =
            boost::lexical_cast<std::string>( argv[2+s] );
    
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

        //--------------------------------------------------------------------------
        // Make cut cell structure (volume)
        base::cut::generateCutCells( mesh, levelSet, cells, s>0 );

        // Make cut cell structure (surface)
        base::cut::generateCutCells( boundaryMesh, levelSet, surfCells, s>0 );

        // merge level set
        if ( s == 0 ) levelSetIntersection = levelSet;
        else 
            std::transform( levelSetIntersection.begin(), levelSetIntersection.end(),
                            levelSet.begin(),
                            levelSetIntersection.begin(),
                            boost::bind( base::cut::setIntersection<dim>, _1, _2 ) );
        

    }

    //--------------------------------------------------------------------------
    // Extract distances, closestPoints and location flags from level set data
    std::vector<double>           distances;
    std::vector<LevelSet::VecDim> closestPoints;
    std::vector<bool>             location;
    std::vector<std::size_t>      closestElements;
    std::vector<double>           distancesToPlane;
    {
        std::transform( levelSetIntersection.begin(), levelSetIntersection.end(),
                        std::back_inserter( distances ),
                        boost::bind( &LevelSet::getSignedDistance, _1 ) );

        std::transform( levelSetIntersection.begin(), levelSetIntersection.end(),
                        std::back_inserter( closestPoints ),
                        boost::bind( &LevelSet::getClosestPoint, _1 ) );

        std::transform( levelSetIntersection.begin(), levelSetIntersection.end(),
                        std::back_inserter( location ),
                        boost::bind( &LevelSet::isInterior, _1 ) );

        std::transform( levelSetIntersection.begin(), levelSetIntersection.end(),
                        std::back_inserter( closestElements ),
                        boost::bind( &LevelSet::getClosestElement, _1 ) );

        std::transform( levelSetIntersection.begin(), levelSetIntersection.end(),
                        std::back_inserter( distancesToPlane ),
                        boost::bind( &LevelSet::getDistanceToPlane, _1 ) );
    }

    //--------------------------------------------------------------------------
    // Evaluate distances at cut cell nodes
    std::vector<double> cutCellDistance;
    base::cut::evaluateCutCellNodeDistances( mesh, levelSetIntersection, cells,
                                             cutCellDistance );
    //--------------------------------------------------------------------------
    // VTK file -- bulk mesh
    {
        // compute element areas in parameter coordinates (in- and outside)
        std::vector<double> areaIn, areaOut;
        std::transform( cells.begin(), cells.end(), std::back_inserter( areaIn ),
                        boost::bind( &Cell::parameterArea, _1, true ) );
        std::transform( cells.begin(), cells.end(), std::back_inserter( areaOut ),
                        boost::bind( &Cell::parameterArea, _1, false ) );
    
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( mesh );
        vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
        vtkWriter.writePointData( closestPoints.begin(), closestPoints.end(), "points" );
        vtkWriter.writePointData( location.begin(), location.end(), "location" );
        vtkWriter.writePointData( closestElements.begin(), closestElements.end(), "elemID" );
        vtkWriter.writePointData( distancesToPlane.begin(), distancesToPlane.end(), "d2p" );


        vtkWriter.writeCellData( areaIn.begin(),  areaIn.end(),  "areaIn" );
        vtkWriter.writeCellData( areaOut.begin(), areaOut.end(), "areaOut" );
        
        vtk.close();
    }

    //--------------------------------------------------------------------------
    // VTK file -- extracted surface
    {
        typedef base::cut::SurfaceSimplexMesh<Mesh>::Type SurfaceSimplexMesh;
        SurfaceSimplexMesh surfaceSimplexMesh;
        base::cut::extractSurfaceMeshFromCutCells( mesh, cells, surfaceSimplexMesh );

        const std::string vtkFile = baseName + ".surf.vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( surfaceSimplexMesh );
        vtkWriter.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
    }

    //--------------------------------------------------------------------------
    // VTK file -- extracted volumina
    {
        typedef base::cut::VolumeSimplexMesh<Mesh>::Type VolumeSimplexMesh;
        VolumeSimplexMesh volumeSimplexMeshIn;
        base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMeshIn,
                                                  true );

        const std::string vtkFileIn = baseName + ".volin.vtk";
        std::ofstream vtkIn( vtkFileIn.c_str() );
        base::io::vtk::LegacyWriter vtkWriterIn( vtkIn );
        vtkWriterIn.writeUnstructuredGrid( volumeSimplexMeshIn );
        vtkWriterIn.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );

        VolumeSimplexMesh volumeSimplexMeshOut;
        base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMeshOut,
                                                  false );
        
        const std::string vtkFileOut = baseName + ".volout.vtk";
        std::ofstream vtkOut( vtkFileOut.c_str() );
        base::io::vtk::LegacyWriter vtkWriterOut( vtkOut );
        vtkWriterOut.writeUnstructuredGrid( volumeSimplexMeshOut );
        vtkWriterOut.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
    }

    //--------------------------------------------------------------------------
    {
        typedef base::cut::VolumeSimplexMesh<BoundaryMesh>::Type BoundarySimplexMesh;
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
    
    //--------------------------------------------------------------------------
    // integrate over volumina
    base::cut::Quadrature<1,shape> quadrature( cells, true );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Mesh> FieldBinder;
    FieldBinder fieldBinder( mesh, mesh );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // integrate
    double volumeIn = 0.;
    base::asmb::simplyIntegrate<FTB>( quadrature, volumeIn, fieldBinder,
                                      base::kernel::Measure<FTB::Tuple>() );


    double volumeOut = 0.;
    quadrature.flipInside();
    base::asmb::simplyIntegrate<FTB>( quadrature, volumeOut, fieldBinder,
                                      base::kernel::Measure<FTB::Tuple>() );
    
    std::cout << "Volume of mesh: " << volumeIn << " + " << volumeOut
              << " = " << volumeIn + volumeOut
              << '\n';

    


    return 0;
}
