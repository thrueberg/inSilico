//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smfReverse.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <fstream>
#include <string>
// boost includes
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
// base includes
#include <base/verify.hpp>
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
// base/fe includes
#include <base/fe/LagrangeElement.hpp>
// base/io includes
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
// tools/converter/smf2xx includes
#include <tools/converter/smf2xx/SMFHeader.hpp>
#include <tools/converter/smf2xx/Conversion.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace smfReverse{

            //------------------------------------------------------------------
            /** Reverse a given array of node pointers which represents the
             *  connectivity of an element. Only surface elements shall be
             *  considered.  Since the reorientation of nodes located in the
             *  interior of faces is rather complicated, this implementation
             *  only considers polynomial degrees of QUADs until quadratic
             *  and of TRIs until cubic. In these extreme cases, only one node
             *  node for the face exists and does not need to be re-oriented.
             *  \tparam NODE    Type of node
             *  \tparam DEGREE  Polynomial degree of the geometry representation
             *  \tparam SHAPE   Geometric shape of the element
             */
            template<typename NODE, unsigned DEGREE, base::Shape SHAPE>
            struct ReverseConnectivity
            {

                typedef base::fe::LagrangeElement<SHAPE,DEGREE> LagrangeElement;

                static const unsigned numNodes       = LagrangeElement::numTotalDoFs;
                static const int      numVertexNodes = LagrangeElement::numVertexDoFs;
                static const int      numEdgeNodes   = LagrangeElement::numEdgeDoFs;

                // apply to array of node pointers
                static void apply( boost::array<NODE*,numNodes>& connec )
                {
                    boost::array<NODE*,numNodes> tmp = connec;

                    // re-orient the vertices
                    for ( int v = 0; v < numVertexNodes; v++ ){
                        connec[v] = tmp[numVertexNodes - v - 1];
                    }

                    // re-orient the edges
                    for ( int e = 0; e < numEdgeNodes; e++ ) {
                        connec[numVertexNodes + e] =
                            tmp[numVertexNodes + numEdgeNodes - e - 1];
                    }

                    // The re-orientation of the face nodes is not so straightforward.
                    // Therefore this function does not handle more than one face
                    // node, which obviously needs no re-orientation
                    if ( SHAPE == base::TRI ) 
                        VERIFY_MSG( (DEGREE < 4),
                                    "Reversal of face nodes not implemented" );

                    if ( SHAPE == base::QUAD )
                        VERIFY_MSG( (DEGREE < 3),
                                    "Reversal of face nodes not implemented" );
                    
                }
            };

            //------------------------------------------------------------------
            /** Read smf file, reverse the element connectivities and write
             *  \tparam SHAPE  Type of element shape
             *  \tparam DEGREE Polynomial degree of the element
             */
            template<base::Shape SHAPE,unsigned DEGREE>
            struct Converter
            {
                static void apply( std::istream& smfIn,
                                   std::ostream& smfOut )
                {
                    // It does not make sense to reverse 3D elements
                    VERIFY_MSG( (base::ShapeDim<SHAPE>::value < 3),
                                "Do not reverse 3D elements" );

                    // Attributes of the mesh
                    static const base::Shape shape   = SHAPE;
                    static const unsigned degree     = DEGREE;
                    static const unsigned    dim     = 3;

                    // Mesh type and object
                    typedef base::Unstructured<shape,degree,dim>  Mesh;
                    Mesh mesh;

                    typedef typename Mesh::Element Element;
                    typedef typename Mesh::Node    Node;
                    static const unsigned numNodes = Element::numNodes;

                    // SMF input
                    base::io::smf::readMesh( smfIn, mesh );

                    // go through elements and modify connectivities
                    typename Mesh::ElementPtrIter eIter = mesh.elementsBegin();
                    typename Mesh::ElementPtrIter eEnd  = mesh.elementsEnd();
                    for (; eIter != eEnd; ++eIter ) {
                        // extract connectivity
                        boost::array<Node*,numNodes> connec;
                        typename Element::NodePtrIter nIter = (*eIter) -> nodesBegin();
                        typename Element::NodePtrIter nEnd  = (*eIter) -> nodesEnd();
                        for ( unsigned n = 0; nIter != nEnd; ++nIter )
                            connec[n++] = (*nIter);

                        // reverse the connectivity
                        ReverseConnectivity<Node,degree,shape>::apply( connec );

                        // set connectivity
                        nIter = (*eIter) -> nodesBegin();
                        for ( unsigned n = 0; nIter != nEnd; ++nIter )
                            (*nIter) = connec[n++];
                    }

                    // SMF output
                    base::io::smf::writeMesh( mesh, smfOut );
                }

            };

    
        } // namespace smfReverse
    } // namespace converter
} // namespace tools

//------------------------------------------------------------------------------
/** Read smf formatted file, create a temporary mesh and write gnuplot data file
 */
int main( int argc, char * argv[] )
{
    namespace smfReverse = tools::converter::smfReverse;
    namespace smf2xx    = tools::converter::smf2xx;
    
    // Sanity check of the number of input arguments
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0]
                  << " file.smf \n\n";
        return 0;
    }

    // Name of smf input file, its basename and the data output file name
    const std::string smfFileIn   = boost::lexical_cast<std::string>( argv[1] );
    const std::string base        = smfFileIn.substr(0, smfFileIn.find( ".smf") );
    const std::string smfFileOut  = base + ".reverse.smf";

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
    smfOut << "# Generated by smfReverse \n";

    // Call generic conversion helper
    smf2xx::Conversion< smfReverse::Converter >::apply( elementShape,
                                                        elementNumPoints,
                                                        smfIn, smfOut );

    // Close the streams
    smfIn.close();
    smfOut.close();
    
    return 0;
}
