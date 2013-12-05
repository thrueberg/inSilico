//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   bruteForce.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_bruteForce_hpp
#define base_cut_bruteForce_hpp

//------------------------------------------------------------------------------
#include <base/cut/LevelSet.hpp>
#include <base/cut/distanceToElement.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename DOMAINMESH, typename SURFMESH>
        void bruteForce( const DOMAINMESH& domainMesh,
                         const SURFMESH&   surfaceMesh,
                         const bool        isSigned, 
                         std::vector< base::cut::LevelSet<DOMAINMESH::Node::dim> >&
                         levelSet );
    }
}

//------------------------------------------------------------------------------
/** Brute-force computation of the distance function between domain and surface
 *  meshes.
 *  For every surface element, its signed distance to every node in the mesh is
 *  computed.
 */
template<typename DOMAINMESH, typename SURFMESH>
void base::cut::bruteForce( const DOMAINMESH& domainMesh,
                            const SURFMESH&   surfaceMesh,
                            const bool        isSigned, 
                            std::vector< base::cut::LevelSet<DOMAINMESH::Node::dim> >&
                            levelSets )
{
    typedef base::cut::LevelSet<DOMAINMESH::Node::dim> LevelSet;

    
    
    // Fill level set with node data
    typename DOMAINMESH::NodePtrConstIter nIter = domainMesh.nodesBegin();
    typename DOMAINMESH::NodePtrConstIter nEnd  = domainMesh.nodesEnd();
    for ( ; nIter != nEnd; ++nIter ) {

        typename LevelSet::VecDim x;
        (*nIter) -> getX( &(x[0]) );

        LevelSet ls( x );
        levelSets.push_back( ls );
    }

    // Go through the surface elements
    typename SURFMESH::ElementPtrConstIter seIter = surfaceMesh.elementsBegin();
    typename SURFMESH::ElementPtrConstIter seEnd  = surfaceMesh.elementsEnd();
    for ( ; seIter != seEnd; ++seIter ) {

        // get pointer to surface element
        const typename SURFMESH::Element* surfEp = *seIter;

        // go through all level set points
        for ( unsigned el = 0; el < levelSets.size(); el++ ) {

            // get level set datum
            LevelSet ls = levelSets[el];

            // compute distance data
            base::cut::distanceToElement( surfEp, isSigned, ls );

            // set level set (assignment does not work?)
            {
                if ( ls.isInterior() ) levelSets[el].setInterior();
                else                   levelSets[el].setExterior();

                levelSets[el].setClosestPoint(   ls.getClosestPoint() );
                levelSets[el].setClosestElement( ls.getClosestElement() );
            }
        }

    }

    return;
}

#endif
