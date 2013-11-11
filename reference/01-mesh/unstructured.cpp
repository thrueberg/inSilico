//[std]{ std includes
#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
//[std]}

//[boost] boost includes
#include <boost/lexical_cast.hpp> // lexical cast between objects

//------------------------------------------------------------------------------
//[headerMesh]{
// base::Shape and the traits object
#include <base/shape.hpp>
// Unstructured mesh
#include <base/Unstructured.hpp>
// Mesh boundary definition
#include <base/mesh/MeshBoundary.hpp>
// Creation of a boundary mesh
#include <base/mesh/generateBoundaryMesh.hpp>
//[headerMesh]}
//[headerIO]{
// IO helper
#include <base/io/Format.hpp>
// SMF Reader
#include <base/io/smf/Reader.hpp>
// Write a VTK file
#include <base/io/vtk/LegacyWriter.hpp>
// Write an SMF file
#include <base/io/smf/Writer.hpp>
//[headerIO]}

STATIC_ASSERT_MSG( (SPACEDIM>1) and (SPACEDIM<4), "Inapt choice of dimension" );

//------------------------------------------------------------------------------
// Message to the user if called without arguments
template<base::Shape SHAPE, unsigned GEOMDEG>
void usageMessage( char* arg )
{
    std::cout << "Usage: " << arg << " file.smf \n\n"
              << "[Compiled for shape=" << base::ShapeName<SHAPE>::apply()
              << " and geometry degree=" << GEOMDEG << "]\n";
}

//------------------------------------------------------------------------------
//[main]{
int main( int argc, char* argv[] )
{
    //[main]}

    //[attributes]{
    //--------------------------------------------------------------------------
    // Definition of the spatial dimension
    const unsigned dim = SPACEDIM;
    // Defintion of an element shape
    const base::Shape elemShape = base::SimplexShape<dim>::value;
    // Definition of the polynomial degree of the geometry approximation
    const unsigned geomDeg = 1;
    //[attributes]}

    //[userInput]{
    // Check the number of input arguments
    if ( argc != 2 ) { // note argv[0] is the program name itself
        usageMessage<elemShape,geomDeg>( argv[0] );
        return 0;
    }
    // convert input argument to a string object
    const std::string smfFileName = boost::lexical_cast<std::string>( argv[1] );
    // extract the basename of the file
    const std::string baseName = base::io::baseName( smfFileName, ".smf" );
    //[userInput]}

    //[mesh]{
    typedef base::Unstructured<elemShape,geomDeg>            Mesh;

    // Create an empty mesh object
    Mesh mesh;
    //[mesh]}

    //--------------------------------------------------------------------------
    {
        //[smfInp]{  Input from an smf file
        // input file stream from smf file
        std::ifstream smf( smfFileName.c_str() );
        // read smf mesh
        base::io::smf::readMesh( smf, mesh );
        // close the stream
        smf.close();
        //[smfInp]}
    }

    //--------------------------------------------------------------------------
    //[boundary]{ Boundary creation
    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create list of <Element,faceNo> pairs
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
    }
    //[boundary]}

    //--------------------------------------------------------------------------
    // Output functions
    {

        // VTK file of the input mesh
        {
            //[vtk]{
            // file name for vtk mesh file
            const std::string vtkMeshFileName = baseName + ".vtk";
            // output file stream 
            std::ofstream vtk( vtkMeshFileName.c_str() );
            // create a VTK writer
            base::io::vtk::LegacyWriter vtkWriter( vtk );
            // write the mesh
            vtkWriter.writeUnstructuredGrid( mesh );
            // close the stream
            vtk.close();
            //[vtk]}
        }

        // SMF file of the boundary mesh
        {
            //[smfOut]{
            // file name for boundary mesh output
            const std::string smfBoundaryFileName = baseName + "_boundary.smf";
            // output stream
            std::ofstream smf( smfBoundaryFileName.c_str() );
            // write to smf
            base::io::smf::writeMesh( boundaryMesh, smf );
            // close stream
            smf.close();
            //[smfOut]}
        }
        
    }


    return 0;
}
