//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ExtractElementFaces.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_extractelementfaces_hpp
#define base_mesh_extractelementfaces_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/shape.hpp>
// base/mesh includes
#include <base/mesh/ElementFaces.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename ELEMENT>
        struct ExtractElementVertices;

        template<base::Shape SHAPE, base::NFace NFACE>
        struct ExtractElementFaces;
    }
}

//------------------------------------------------------------------------------
/** Extract the IDs of the element nodes which fall on vertices.
 *  Note that a hierarchic node ordering is assumed and, therefore, the first
 *  \f$ n_V \f$ nodes are the vertex nodes, \f$ n_V \f$ being the number of
 *  vertices associated with the shape of the element.
 *  \tparam ELEMENT Type of element to extract from
 */
template<typename ELEMENT>
struct base::mesh::ExtractElementVertices
{
    //! Number of vertices for given shape
    static const unsigned numVertices =
        base::NumNFaces<ELEMENT::shape,base::VERTEX>::value;

    //! Read number of vertices node IDs from element
    static void apply( const ELEMENT * ep, 
                       boost::array<std::size_t,numVertices>& vertexIDs )
    {
        typename ELEMENT::NodePtrConstIter node = ep -> nodesBegin();
        typename ELEMENT::NodePtrConstIter last = ep -> nodesEnd();
                    
        for ( unsigned n = 0;
              (node != last) and (n < numVertices); ++node, n++ ) {
            vertexIDs[n] = (*node) -> getID();
        }
    }

};
        

//------------------------------------------------------------------------------
/** Extract the connectivities of faces from a given element.
 *  For instance, the quadrilateral faces of a HEX or the edges of a TET might
 *  be of interest. In that case, this object does the following steps:
 *   
 *  1. Extract the IDs of nodes which fall on vertices 
 *  2. Use base::mesh::ElementFaces in order to extract the connectivity
 *     of the faces
 *  3. Permutate the connectivity in order to have a unique face identification
 *
 *  \tparam SHAPE    Type of shape to extract from (QUAD, TET, ... )
 *  \tparam NFACE    Type of face to extract (EDGE, FACE, ... )
 */
template<base::Shape SHAPE, base::NFace NFACE>
struct base::mesh::ExtractElementFaces
{
    //! Traits for array sizes
    typedef base::mesh::ElementFaceTraits<SHAPE,NFACE,std::size_t> ElementFaceTraits;

    //! Object for face extraction
    typedef base::mesh::ElementFaces<SHAPE,NFACE> ElementFaces;
            
    //! Number of vertices on that face
    static const unsigned numFaceVertices = ElementFaceTraits::numVertices;
                
    //! Define array of indices as face
    typedef typename ElementFaceTraits::Type Face;

    //! Number of vertices of the element shape
    static const unsigned numElementVertices =
        base::NumNFaces<SHAPE,base::VERTEX>::value;
    
    //--------------------------------------------------------------------------
    /** Main action: extract face from given element.
     *  \tparam     NUMNODES   Number of nodes of the element
     *  \param[in]  vertexIDs  Array of global node IDs of the element
     *  \param[out] nFaces     Vector with face  connectivities
     */
    static void apply( const boost::array<std::size_t,
                                          numElementVertices>& vertexIDs,
                       std::vector<Face>& nFaces )
    {
        //! Go through all faces of interest
        for ( unsigned f = 0; f < ElementFaces::numFaces; f ++ ) {

            //! A face
            Face face;

            //! Go through all vertices of the face
            for ( unsigned v = 0; v < numFaceVertices; v ++ ) {
                const std::size_t vertexID = vertexIDs[ ElementFaces::index( f, v ) ];
                ElementFaceTraits::assignIndex( face, v, vertexID );
            }
                    
            // store in provided container
            nFaces.push_back( face );

        }

    }
};


#endif
