//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   GenerateMeshFromFaces.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_generatemeshfromfaces_hpp
#define base_mesh_generatemeshfromfaces_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
#include <iterator>
// base/mesh includes
#include <base/mesh/Unstructured.hpp>
#include <base/mesh/ElementFaces.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename MESH, base::NFace NFACE>
        class GenerateMeshFromFaces;
        
    }
}

//------------------------------------------------------------------------------
template<typename MESH, base::NFace NFACE>
class base::mesh::GenerateMeshFromFaces
{
public:
    //! Element face traits object
    typedef base::mesh::ElementFaceTraits<MESH::Element::shape,
                                          NFACE,std::size_t>  ElementFaceTraits;
    //! Type of face storage
    typedef typename ElementFaceTraits::Type  Face;

    //! Number of vertices of a face
    static const unsigned numFaceVertices = ElementFaceTraits::numVertices;

    //! @name Typedefs for face mesh
    //@{
    typedef typename MESH::Node                                Node;
    typedef base::LagrangeShapeFun<1,ElementFaceTraits::shape> GeomFun;
    typedef base::mesh::Element<Node,GeomFun>                  Element;
    typedef base::mesh::Unstructured<Element>                  FaceMesh;
    //@}

    //! Get required faces of every element and store the unqiue ones
    static void apply( const MESH & mesh,
                       const std::vector<Face> & faces,
                       FaceMesh & faceMesh )
    {
        // 1. Allocate space for nodes and face elements
        const std::size_t numNodes = std::distance( mesh.nodesBegin(),
                                                    mesh.nodesEnd() );
        const std::size_t numElements = faces.size();
        faceMesh.allocate( numNodes, numElements );

        // 2. Copy nodes (not just the pointers!!!)
        typename MESH::NodePtrConstIter sourceNode = mesh.nodesBegin();
        typename MESH::NodePtrConstIter finalNode  = mesh.nodesEnd();
        typename FaceMesh::NodePtrIter  copyNode   = faceMesh.nodesBegin();
        for ( ; sourceNode != finalNode; ++sourceNode, ++copyNode ) {

            // Copy the global ID
            (*copyNode) -> setID( (*sourceNode) -> getID() );

            // Copy the coordinates
            std::vector<double> x( Node::dim );
            (*sourceNode) -> getX( x.begin() );
            (*copyNode)   -> setX( x.begin() );
        }
        

        // 3. Generate new elements from faces
        typename FaceMesh::ElementPtrIter elemIter = faceMesh.elementsBegin();
        typename FaceMesh::ElementPtrIter elemEnd  = faceMesh.elementsEnd();
        typename std::vector<Face>::const_iterator faceIter = faces.begin();
        for ( ; elemIter != elemEnd; ++elemIter, ++faceIter ) {

            // Extract face's node IDs
            Face nodeIDs;
            for ( unsigned v = 0; v < nodeIDs.size(); v ++ )
                nodeIDs[v] = faceIter -> at( v );

            // Set face element connectivity
            typename FaceMesh::Element::NodePtrIter node =
                (*elemIter) -> nodesBegin();
            for ( unsigned v = 0; v < Face::size(); v++, ++node ) {

                (*node) = mesh.nodePtr( nodeIDs[v] );
            }
            
        }
        return;
    }
    
};
#endif
