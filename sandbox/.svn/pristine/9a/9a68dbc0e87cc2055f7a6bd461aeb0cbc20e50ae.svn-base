//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   CollapseMesh.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef tfm_collapsemesh_hpp
#define tfm_collapsemesh_hpp

//------------------------------------------------------------------------------
#include <iterator>
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace tfm{
    template<typename MESH>
    void collapseMesh( const MESH& in, MESH& out, const double tolerance );
}

//------------------------------------------------------------------------------
template<typename MESH>
void tfm::collapseMesh( const MESH& in, MESH& out, const double tolerance )
{
    //--------------------------------------------------------------------------
    // fill with IDs of the input mesh
    typename MESH::NodePtrConstIter nIter = in.nodesBegin();
    typename MESH::NodePtrConstIter nEnd  = in.nodesEnd();

    // do not act on an empty mesh
    if ( nIter == nEnd ) return;

    std::vector<std::size_t> nodeIDs;
    std::vector<bool>        isUnique;


    for ( std::size_t n = 0; nIter != nEnd; ++nIter, n++ ) {

        // make sure the IDs are consecutive
        VERIFY_MSG( ( (*nIter) -> getID() ) == n,
                    "Works only for consecutive node IDs" );
        
        nodeIDs.push_back(  n    );
        isUnique.push_back( true );

    }

    //--------------------------------------------------------------------------
    // check nodes for geometric coincidence
    nIter = in.nodesBegin();
    nEnd  = in.nodesEnd();
    typename MESH::NodePtrConstIter nEnd2  = nEnd;
    nEnd--;
    
    for ( std::size_t m = 0; nIter != nEnd; ++nIter, m++ ) {

        if ( isUnique[m] ) {
        
            typename MESH::Node::VecDim x1 = (*nIter) -> getX();

            typename MESH::NodePtrConstIter nIter2 = nIter;
            nIter2++;
            for ( std::size_t n = m+1; nIter2 != nEnd2; ++nIter2, n++) {

                // distance between nodes
                typename MESH::Node::VecDim diff =
                    ( (*nIter2) -> getX() ) - x1;

                if ( base::norm( diff ) < tolerance ) {
                    nodeIDs[  n ] = m;
                    isUnique[ n ] = false;
                }
            }
        }
        
    }

    //--------------------------------------------------------------------------
    // create new mesh
    
    // number of unique nodes
    const std::size_t numUniqueNodes = std::count( isUnique.begin(),
                                                   isUnique.end(), true );

    // allocate space for new mesh
    out.allocate( numUniqueNodes,
                  std::distance( in.elementsBegin(), in.elementsEnd() ) );

    // copy nodes
    nIter = in.nodesBegin();
    nEnd  = in.nodesEnd();
    typename MESH::NodePtrIter nIterNew = out.nodesBegin();

    std::vector<std::size_t> newNodeIDs( nodeIDs.size(), base::invalidInt );
    
    std::size_t nodeID = 0;
    for ( std::size_t n = 0;  nIter != nEnd; nIter++, n++ ) {

        if ( isUnique[n] ) {

            const typename MESH::Node::VecDim x = (*nIter) -> getX();

            (*nIterNew) -> setX(  x );
            (*nIterNew) -> setID( nodeID );

            newNodeIDs[ n ] = nodeID;
            
            nIterNew++;
            nodeID++;
        }
    }

    // copy elements
    typename MESH::ElementPtrConstIter eIter = in.elementsBegin();
    typename MESH::ElementPtrConstIter eEnd  = in.elementsEnd();
    typename MESH::ElementPtrIter      eIterNew = out.elementsBegin();
    for ( ; eIter != eEnd; ++eIter, ++eIterNew ) {

        // copy ID
        (*eIterNew) -> setID( (*eIter) -> getID() );

        typename MESH::Element::NodePtrConstIter eNIter = (*eIter) -> nodesBegin();
        typename MESH::Element::NodePtrConstIter eNEnd  = (*eIter) -> nodesEnd();
        typename MESH::Element::NodePtrIter      eNIterNew = (*eIterNew) -> nodesBegin();

        for ( ; eNIter != eNEnd; ++eNIter, ++eNIterNew ) {

            const std::size_t oldNodeID = (*eNIter) -> getID();
            const std::size_t newNodeID = newNodeIDs[ nodeIDs[ oldNodeID ] ];

            *eNIterNew = out.nodePtr( newNodeID );
        }

    }

    return;
    
}
#endif
