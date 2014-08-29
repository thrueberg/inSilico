#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
#include <boost/lexical_cast.hpp> // lexical cast between objects

//------------------------------------------------------------------------------
// Structured mesh
#include <base/Structured.hpp>
// Mesh boundary definition
#include <base/mesh/MeshBoundary.hpp>
// Creation of a boundary mesh
#include <base/mesh/generateBoundaryMesh.hpp>
// sampling, e.g. refinement
#include <base/mesh/sampleStructured.hpp>
// IO helper
#include <base/io/Format.hpp>
// SMF Reader
#include <base/io/sgf/Reader.hpp>
// Write a VTK file
#include <base/io/vtk/LegacyWriter.hpp>
// Write an SMF file
#include <base/io/smf/Writer.hpp>

STATIC_ASSERT_MSG( (SPACEDIM>1) and (SPACEDIM<4), "Inapt choice of dimension" );

//------------------------------------------------------------------------------
// Message to the user if called without arguments
namespace ref01{
    
     //! Create a message how to call this program
    template<unsigned DIM, unsigned GEOMDEG>
    void usageMessage2( char* arg )
    {
        std::cout << "Usage: " << arg << " file.sgf \n\n"
                  << "[Compiled for dim=" << DIM
                  << " and geometry degree=" << GEOMDEG << "]\n";
    }

    int structured( int argc, char* argv[] );
}

//------------------------------------------------------------------------------
/** In- and output of a structured mesh.
 *  Does the same job as ref01::unstructured, but is using the SGF format
 *  (base::io::sgf) for structured meshes. These are meshes in which every
 *  element is a hypercube, every interior node has 2 x DIM direct neighbors
 *  (as connected by edges) and a lexicographic ordering of the indices is used.
 *  Structured meshes have many convenience advantages and have therefore been
 *  introduced to inSilico. Nevertheless, almost all applications are based
 *  on unstructured meshes for the sake of generality and the power of structured
 *  meshes has not been fully exploited.
 *
 *  \param[in] argc Number of command line arguments
 *  \param[in] argv Values of command line arguments
 */
int ref01::structured( int argc, char* argv[] )
{
    //--------------------------------------------------------------------------
    // Definition of the spatial dimension
    const unsigned dim = SPACEDIM;
    // Definition of the polynomial degree of the geometry approximation
    const unsigned geomDeg = 1;
    // The output can be with a higher resolution
    const unsigned resolution = 2;
    // Boundary mesh shall be composed of simplex elements
    const bool boundarySimplex = true;
    // Check the number of input arguments
    if ( argc != 2 ) { // note argv[0] is the program name itself
        usageMessage2<dim,geomDeg>( argv[0] );
        return 0;
    }
    // convert input argument to a string object
    const std::string sgfFileName = boost::lexical_cast<std::string>( argv[1] );
    // extract the basename of the file
    const std::string baseName = base::io::baseName( sgfFileName, ".sgf" );

    typedef base::Structured<dim,geomDeg>          Grid;
    // Create an empty mesh object
    Grid grid;
    
    //--------------------------------------------------------------------------
    {
        // input file stream from smf file
        std::ifstream sgf( sgfFileName.c_str() );
        // read sgf mesh
        base::io::sgf::readGrid( sgf, grid );
        // close the stream
        sgf.close();
    }

    //--------------------------------------------------------------------------
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary gridBoundary;
    gridBoundary.create<dim>( grid.gridSizes() );

    typedef base::mesh::BoundaryMeshBinder<Grid::Element,
                                           boundarySimplex>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( gridBoundary.begin(),
                                          gridBoundary.end(),
                                          grid, boundaryMesh );
    }

    //--------------------------------------------------------------------------
    // Output functions
    {

        // VTK file of the input mesh
        {
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
        }

        // SMF file of the boundary mesh
        {
            // file name for boundary mesh output
            const std::string smfBoundaryFileName = baseName + "_boundary.smf";
            // output stream
            std::ofstream smf( smfBoundaryFileName.c_str() );
            // write to smf
            base::io::smf::writeMesh( boundaryMesh, smf );
            // close stream
            smf.close();
        }
        
    }


    return 0;
}

//------------------------------------------------------------------------------
// Delegate the call 
int main( int argc, char* argv[] )
{
    return ref01::structured( argc, argv );
}
