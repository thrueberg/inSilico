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

#include <base/cut/Quadrature.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/kernel/Measure.hpp>

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    // User input
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf surf.smf \n\n";
        return -1;
    }
        
    const std::string smfFile          = boost::lexical_cast<std::string>( argv[1] );
    const std::string surfMeshFileName = boost::lexical_cast<std::string>( argv[2] );
    const std::string baseName         = base::io::baseName( smfFile, ".smf" );

    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
    const unsigned    dim       = 2;
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

    //--------------------------------------------------------------------------
    // Surface mesh
    typedef base::Unstructured<surfShape,1,dim>    SurfMesh;

    SurfMesh surfMesh;
    {
        std::ifstream smf( surfMeshFileName.c_str() );
        base::io::smf::readMesh( smf, surfMesh );
        smf.close();
    }

    //--------------------------------------------------------------------------
    // Compute the level set data
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::bruteForce( mesh, surfMesh, isSigned, levelSet );

    //--------------------------------------------------------------------------
    // Make cut cell structure
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells;
    base::cut::generateCutCells( mesh, levelSet, cells );

    //--------------------------------------------------------------------------
    {
        typedef base::cut::SurfaceSimplexMesh<Mesh>::Type SurfaceSimplexMesh;
        SurfaceSimplexMesh surfaceSimplexMesh;
        base::cut::extractSurfaceMeshFromCutCells( mesh, cells, surfaceSimplexMesh );

        const std::string smfSurf = baseName + ".surf.smf";
        std::ofstream smf( smfSurf.c_str() );
        base::io::smf::writeMesh( surfaceSimplexMesh, smf );
        smf.close();
    }
    
    //--------------------------------------------------------------------------
    {
        typedef base::cut::VolumeSimplexMesh<Mesh>::Type VolumeSimplexMesh;
        VolumeSimplexMesh volumeSimplexMesh;
        base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMesh,
                                                  true );

        const std::string smfVol = baseName + ".volin.smf";
        std::ofstream smf( smfVol.c_str() );
        base::io::smf::writeMesh( volumeSimplexMesh, smf );
        smf.close();
    }

    //--------------------------------------------------------------------------
    {
        typedef base::cut::VolumeSimplexMesh<Mesh>::Type VolumeSimplexMesh;
        VolumeSimplexMesh volumeSimplexMesh;
        base::cut::extractVolumeMeshFromCutCells( mesh, cells, volumeSimplexMesh,
                                                  false );

        const std::string smfVol = baseName + ".volout.smf";
        std::ofstream smf( smfVol.c_str() );
        base::io::smf::writeMesh( volumeSimplexMesh, smf );
        smf.close();
    }

    //--------------------------------------------------------------------------
    // Extract distances, closestPoints and location flags from level set data
    std::vector<double>           distances;
    std::vector<LevelSet::VecDim> closestPoints;
    std::vector<bool>             location;
    std::vector<std::size_t>      closestElements;
    std::vector<double>           distancesToPlane;
    {
        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( distances ),
                        boost::bind( &LevelSet::getSignedDistance, _1 ) );

        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( closestPoints ),
                        boost::bind( &LevelSet::getClosestPoint, _1 ) );

        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( location ),
                        boost::bind( &LevelSet::isInterior, _1 ) );

        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( closestElements ),
                        boost::bind( &LevelSet::getClosestElement, _1 ) );

        std::transform( levelSet.begin(), levelSet.end(),
                        std::back_inserter( distancesToPlane ),
                        boost::bind( &LevelSet::getDistanceToPlane, _1 ) );
    }

    // compute element areas in parameter coordinates (in- and outside)
    std::vector<double> areaIn, areaOut;
    std::transform( cells.begin(), cells.end(), std::back_inserter( areaIn ),
                    boost::bind( &Cell::parameterArea, _1, true ) );
    std::transform( cells.begin(), cells.end(), std::back_inserter( areaOut ),
                    boost::bind( &Cell::parameterArea, _1, false ) );
    
    //--------------------------------------------------------------------------
    // output to a VTK file
    {
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
