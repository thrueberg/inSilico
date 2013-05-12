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
// Lagrangian shape functions
#include <base/LagrangeShapeFun.hpp>
// Node definition
#include <base/mesh/Node.hpp>
// Element definition
#include <base/mesh/Element.hpp>
// Unstructured mesh
#include <base/mesh/Unstructured.hpp>
// Mesh boundary definition
#include <base/mesh/MeshBoundary.hpp>
// Creation of a boundary mesh
#include <base/mesh/CreateBoundaryMesh.hpp>
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
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//[main]{
int main( int argc, char* argv[] )
{
    // Check the number of input arguments
    if ( argc != 2 ) { // note argv[0] is the program name itself
        std::cout << "Usage:  " << argv[0] << " file.smf \n";
        return 0;
    }
    
    // convert input argument to a string object
    const std::string smfFileName = boost::lexical_cast<std::string>( argv[1] );
    // extract the basename of the file
    const std::string baseName = base::io::baseName( smfFileName, ".smf" );
    //[main]}

    //[attributes]{
    //--------------------------------------------------------------------------
    // Defintion of an element shape
    const base::Shape elemShape = base::TRI;
    // Definition of the polynomial degree of the geometry approximation
    const unsigned geomDeg = 1;

    //--------------------------------------------------------------------------
    // spatial dimension dictated by element shape (not necessarily)
    const unsigned dim = base::ShapeDim<elemShape>::value;
    //[attributes]}
    //[mesh]{
    // type of a geometry node
    typedef base::mesh::Node<dim>                        Node;
    // type of shape function for geometry approximation
    typedef base::LagrangeShapeFun<geomDeg,elemShape>    GeomFun;
    // type of element 
    typedef base::mesh::Element<Node,GeomFun>            Element;
    // type of mesh
    typedef base::mesh::Unstructured<Element>            Mesh;

    // Create an empty mesh object
    Mesh mesh;
    //[mesh]}

    //--------------------------------------------------------------------------
    {
        //[smfInp]{  Input from an smf file
        // input file stream from smf file
        std::ifstream smf( smfFileName.c_str() );
        // smf-reader object
        base::io::smf::Reader<Mesh> smfReader;
        // read into mesh
        smfReader( mesh, smf );
        // close the stream
        smf.close();
        //[smfInp]}
    }

    //--------------------------------------------------------------------------
    //[boundary]{ Boundary creation
    typedef base::mesh::CreateBoundaryMesh<Element> CreateBoundaryMesh;
    CreateBoundaryMesh::BoundaryMesh boundaryMesh;
    {
        // Create list of <Element,faceNo> pairs
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        // Create a real mesh object from this list
        CreateBoundaryMesh::apply( meshBoundary.boundaryBegin(),
                                   meshBoundary.boundaryEnd(),
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
            // smf-writer object
            base::io::smf::Writer<CreateBoundaryMesh::BoundaryMesh> smfWriter;
            // write mesh to stream
            smfWriter( boundaryMesh, smf );
            // close stream
            smf.close();
            //[smfOut]}
        }
        
    }


    return 0;
}
