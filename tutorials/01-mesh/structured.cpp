//[std]{ std includes
#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
//[std]}

//[boost] boost includes
#include <boost/lexical_cast.hpp> // lexical cast between objects

//------------------------------------------------------------------------------
//[headerMesh]{
// Structured mesh
#include <base/Structured.hpp>
// Mesh boundary definition
#include <base/mesh/MeshBoundary.hpp>
// Creation of a boundary mesh
#include <base/mesh/generateBoundaryMesh.hpp>
// sampling, e.g. refinement
#include <base/mesh/sampleStructured.hpp>
//[headerMesh]}
//[headerIO]{
// IO helper
#include <base/io/Format.hpp>
// SMF Reader
#include <base/io/sgf/Reader.hpp>
// Write a VTK file
#include <base/io/vtk/LegacyWriter.hpp>
// Write an SMF file
#include <base/io/smf/Writer.hpp>
//[headerIO]}

STATIC_ASSERT_MSG( (SPACEDIM>1) and (SPACEDIM<4), "Inapt choice of dimension" );

//------------------------------------------------------------------------------
// Message to the user if called without arguments
template<unsigned DIM, unsigned GEOMDEG>
void usageMessage( char* arg )
{
    std::cout << "Usage: " << arg << " file.sgf \n\n"
              << "[Compiled for dim=" << DIM
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
    // Definition of the polynomial degree of the geometry approximation
    const unsigned geomDeg = 1;
    // The output can be with a higher resolution
    const unsigned resolution = 2;
    // Boundary mesh shall be composed of simplex elements
    const bool boundarySimplex = true;
    //[attributes]}

    //[userInput]{
    // Check the number of input arguments
    if ( argc != 2 ) { // note argv[0] is the program name itself
        usageMessage<dim,geomDeg>( argv[0] );
        return 0;
    }
    // convert input argument to a string object
    const std::string sgfFileName = boost::lexical_cast<std::string>( argv[1] );
    // extract the basename of the file
    const std::string baseName = base::io::baseName( sgfFileName, ".sgf" );
    //[userInput]}

    //[grid]{
    typedef base::Structured<dim,geomDeg>          Grid;
    // Create an empty mesh object
    Grid grid;
    //[grid]}
    
    //--------------------------------------------------------------------------
    {
        //[sgfInp]{  Input from an smf file
        // input file stream from smf file
        std::ifstream sgf( sgfFileName.c_str() );
        // read sgf mesh
        base::io::sgf::readGrid( sgf, grid );
        // close the stream
        sgf.close();
        //[sgfInp]}
    }

    //--------------------------------------------------------------------------
    //[boundary]{ Boundary creation
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary gridBoundary;
    gridBoundary.create<dim>( grid.gridSizes() );

    typedef base::mesh::BoundaryMeshBinder<Grid::Element,
                                           boundarySimplex>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( gridBoundary.boundaryBegin(),
                                          gridBoundary.boundaryEnd(),
                                          grid, boundaryMesh );
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
            base::io::vtk::LegacyWriter vtkWriter( vtk, resolution );
            // write the mesh
            vtkWriter.writeStructuredGrid( grid );
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
