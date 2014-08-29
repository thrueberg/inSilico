// std includes
#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
// boost includes
#include <boost/lexical_cast.hpp> // lexical cast between objects
// base::Shape and the traits object
#include <base/shape.hpp>
// Unstructured mesh
#include <base/Unstructured.hpp>
// Mesh boundary definition
#include <base/mesh/MeshBoundary.hpp>
// Creation of a boundary mesh
#include <base/mesh/generateBoundaryMesh.hpp>
// IO helper
#include <base/io/Format.hpp>
// SMF Reader
#include <base/io/smf/Reader.hpp>
// Write a VTK file
#include <base/io/vtk/LegacyWriter.hpp>
// Write an SMF file
#include <base/io/smf/Writer.hpp>

// Make sure that the dimensions are correct
STATIC_ASSERT_MSG( (SPACEDIM>1) and (SPACEDIM<4), "Inapt choice of dimension" );

//------------------------------------------------------------------------------
// Message to the user if called without arguments
namespace ref01{

    //! Create a message how to call this program
    template<base::Shape SHAPE, unsigned GEOMDEG>
    void usageMessage1( char* arg )
    {
        std::cout << "Usage: " << arg << " file.smf \n\n"
                  << "[Compiled for shape=" << base::ShapeName<SHAPE>::apply()
                  << " and geometry degree=" << GEOMDEG << "]\n";
    }

    int unstructured( int argc, char* argv[] );

}

//------------------------------------------------------------------------------
/** In- and output of an unstructured mesh.
 *  This function demonstrates
 *     - how to define and create a mesh
 *     - read a  mesh from a file (base::io::smf)
 *     - extract the boundary from the mesh
 *     - write output in vtk and SMF format
 *
 *  The image shows the reference meshes in this application
 *  \image html simplexMeshes.png "Triangle (left) and tetrahedron mesh (right)"
 *  and the boundary meshes are given in the next picture
 *  \image html boundaryMeshes.png "Line (left) and triangle mesh(right)"
 *
 *  \param[in] argc Number of command line arguments
 *  \param[in] argv Values of command line arguments
 */
int ref01::unstructured( int argc, char* argv[] )
{
    //--------------------------------------------------------------------------
    // Definition of the spatial dimension
    const unsigned dim = SPACEDIM;
    // Defintion of an element shape
    const base::Shape elemShape = base::SimplexShape<dim>::value;
    // Definition of the polynomial degree of the geometry approximation
    const unsigned geomDeg = 1;

    // Check the number of input arguments
    if ( argc != 2 ) { // note argv[0] is the program name itself
        usageMessage1<elemShape,geomDeg>( argv[0] );
        return 0;
    }
    // convert input argument to a string object
    const std::string smfFileName = boost::lexical_cast<std::string>( argv[1] );
    // extract the basename of the file
    const std::string baseName = base::io::baseName( smfFileName, ".smf" );

    typedef base::Unstructured<elemShape,geomDeg>            Mesh;

    // Create an empty mesh object
    Mesh mesh;

    //--------------------------------------------------------------------------
    {
        // input file stream from smf file
        std::ifstream smf( smfFileName.c_str() );
        // read smf mesh
        base::io::smf::readMesh( smf, mesh );
        // close the stream
        smf.close();
    }

    //--------------------------------------------------------------------------
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
            base::io::vtk::LegacyWriter vtkWriter( vtk );
            // write the mesh
            vtkWriter.writeUnstructuredGrid( mesh );
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
    return ref01::unstructured( argc, argv );
}
