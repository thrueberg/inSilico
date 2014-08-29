//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smf2smf.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <map>
// boost includes
#include <boost/lexical_cast.hpp>
// base includes
#include <base/verify.hpp>
#include <base/numbers.hpp>
#include <base/Unstructured.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
// tools includes
#include <tools/converter/smf2xx/SMFHeader.hpp>
#include <tools/converter/smf2xx/Conversion.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace smf2smf{
        
            template<base::Shape SHAPE,unsigned DEGREE>
            struct SMFClean;
            
        }
    }
}

//------------------------------------------------------------------------------
/** Converter writes mesh without orphaned nodes.
    
 *  If a mesh has nodes which are node part of the element connectivity (e.g.,
 *  from the gmesh2smf converter), this can lead to problems in the computation.
 *  An isoparametric mesh would simply generate for every node a degree of
 *  freedom. Degrees of freedom of unconnected nodes stray in the system  and
 *  make it singular. Therefore, this converter is to be used prior to
 *  computation and the orphaned nodes get removed.
 *
 *  The process of removing nodes is done as follows:
 *   1.  Read and create mesh with possible orphans
 *   2.  Collect from mesh all IDs and pointers of nodes which are connected to
 *       an element; at the same time generate a map between give node IDs and
 *       a new consecutive ID number
 *   3.  Set new ID (from the map of 2) to all nodes connected to an element
 *       (and an invalid number to the orphans)
 *   4.  Create a new mesh from the collect node pointers and the connectivity
 *       of the old mesh (but now with new node IDs)
 *   5.  Write output
 *
 *   \tparam SHAPE   The shape of the elements
 *   \tparam DEGREE The degree of the geometry approximation
 */
template<base::Shape SHAPE,unsigned DEGREE>
struct tools::converter::smf2smf::SMFClean
{
    static void apply( std::istream& smfIn, std::ostream& smfOut )
    {
        // Attributes of the mesh
        static const base::Shape shape   = SHAPE;
        static const unsigned degree     = DEGREE;
        static const unsigned    dim     = 3;

        // Mesh type and object
        typedef base::Unstructured<shape,degree,dim>  Mesh;
        Mesh mesh;

        // SMF input
        base::io::smf::readMesh( smfIn, mesh );

        // Collect pointers to all nodes connect to an element
        // establish index map
        std::vector<typename Mesh::Node*> goodNodes;
        std::map<std::size_t,std::size_t> nodeNumberMap;
        typename Mesh::ElementPtrConstIter eIter = mesh.elementsBegin();
        typename Mesh::ElementPtrConstIter eLast = mesh.elementsEnd();
        std::size_t newNodeID = 0;
        for ( ; eIter != eLast; ++eIter ) { // go through all elements

            typename Mesh::Element::NodePtrConstIter nIter = (*eIter) -> nodesBegin();
            typename Mesh::Element::NodePtrConstIter nLast = (*eIter) -> nodesEnd();
            for ( ; nIter != nLast; ++nIter ) { // go through element nodes

                // old ID of the existing mesh
                const std::size_t oldID = (*nIter) -> getID();

                // find this ID in number-map
                std::map<std::size_t,std::size_t>::iterator idIter
                    = nodeNumberMap.find( oldID );

                // if not found assign a new ID to the old one
                if (idIter == nodeNumberMap.end() ) {
                    nodeNumberMap[ oldID ] = newNodeID;
                    newNodeID++;
                    // store good node
                    goodNodes.push_back( *nIter );
                }
            }
        }

        //----------------------------------------------------------------------
        // assign new IDs to nodes
        typename Mesh::NodePtrIter nIter = mesh.nodesBegin();
        typename Mesh::NodePtrIter nLast = mesh.nodesEnd();
        for ( ; nIter != nLast; ++nIter ) { // go through all nodes

            // get old node ID
            const std::size_t oldID = (*nIter) -> getID();

            // find old ID in the number-map
            std::map<std::size_t,std::size_t>::iterator idIter
                = nodeNumberMap.find( oldID );

            // if not found, the node is an orphan :(
            if (idIter == nodeNumberMap.end() ) {
                (*nIter) -> setID( base::invalidInt );
            }
            else { // found, give new ID
                const std::size_t newID = idIter -> second;
                (*nIter) -> setID( newID );
            }
        }

        //----------------------------------------------------------------------
        // Create a new mesh from good nodes and connectivity
        Mesh newMesh;
        newMesh.append( goodNodes.begin(), goodNodes.end(),
                        mesh.elementsBegin(), mesh.elementsEnd() );

        // SMF output
        base::io::smf::writeMesh( newMesh, smfOut );

    }
};

//------------------------------------------------------------------------------
/** Read smf formatted file, create a temporary mesh and write transformed mesh
 */
int main( int argc, char * argv[] )
{
    namespace smf2xx   = tools::converter::smf2xx;
    
    // Sanity check of the number of input arguments
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0]
                  << " fileIn.smf fileOut.smf \n\n";
        return 0;
    }

    // Name of smf input file, its basename and the data output file name
    const std::string smfFileIn  = boost::lexical_cast<std::string>( argv[1] );
    const std::string smfFileOut = boost::lexical_cast<std::string>( argv[2] );

    // Element attributes
    base::Shape elementShape;
    unsigned    elementNumPoints;
    
    {
        // extract data from header
        std::ifstream smf( smfFileIn.c_str() );
        smf2xx::readSMFHeader( smf, elementShape, elementNumPoints );
        smf.close();
    }
    
    // Input and output file streams
    std::ifstream smfIn(   smfFileIn.c_str() );
    std::ofstream smfOut( smfFileOut.c_str() );

    // write to file for traceback
    smfOut << "# Generated by smf2smf \n";

    // Call generic conversion helper
    smf2xx::Conversion<tools::converter::SMFClean>::apply( elementShape,
                                                           elementNumPoints,
                                                           smfIn, smfOut );

    // Close the streams
    smfIn.close();
    smfOut.close();
    
    return 0;
}
