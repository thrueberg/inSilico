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
// std includes
#include <vector>
// base/mesh includes
#include <base/mesh/Size.hpp>
// base/cut includes
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
                         levelSet,
                         const double pointIdentityTolerance = 
                         std::sqrt( std::numeric_limits<double>::epsilon() ) );
    }
}

//------------------------------------------------------------------------------
/** Brute-force computation of the distance function between domain
 *  and surface meshes.
 *  For every surface element, its distance to every node in the mesh is
 *  computed. The distance is either signed or unsigned, depending on a given
 *  flag. Moreover, the closest point on the surface and the index of the
 *  surface element, it lies on, are stored.
 *
 *  \warning This function does not rely on any structure of the domain mesh nor
 *           does it restrict its computation to the vicinity of the surface.
 *           Therefore, it bears the maximal (worst!) complexity by checking for
 *           every domain mesh point every surface element. There is a variety
 *           of methods which accelerate this compuation but are not used here.
 *
 *  \warning Using the nodes of the domain mesh only works properly for either
 *           Lagrangian geometry representations or linear B-splines. Higher
 *           order B-splines do not really have nodes and therefore the geometry
 *           needs to be evaluated at the vertices of a reference element.
 *
 *  \tparam DOMAINMESH  Type of domain mesh (no requirement on structure)
 *  \tparam SURFMESH    Type of surface mesh
 *  \param[in]   domainMesh  Access to the domain mesh
 *  \param[in]   surfaceMesh Access to the surface mesh
 *  \param[in]   isSigned    Flag for signed/unsigned distance function
 *  \param[out]  levelSet    Vector with level set data
 *  \param[in]   pointIdentityTolerance Tolerance for point comparison
 */
template<typename DOMAINMESH, typename SURFMESH>
void base::cut::bruteForce( const DOMAINMESH& domainMesh,
                            const SURFMESH&   surfaceMesh,
                            const bool        isSigned, 
                            std::vector< base::cut::LevelSet<DOMAINMESH::Node::dim> >&
                            levelSet,
                            const double pointIdentityTolerance )
{
    // convenience typedef
    typedef base::cut::LevelSet<DOMAINMESH::Node::dim> LevelSet;
    
    // Fill level set data with node coordinates
    levelSet.resize( 0 );
    typename DOMAINMESH::NodePtrConstIter nIter = domainMesh.nodesBegin();
    typename DOMAINMESH::NodePtrConstIter nEnd  = domainMesh.nodesEnd();
    for ( ; nIter != nEnd; ++nIter ) {

        typename LevelSet::VecDim x;
        (*nIter) -> getX( &(x[0]) );

        LevelSet ls( x );
        levelSet.push_back( ls );
    }

    typename SURFMESH::ElementPtrConstIter seIter = surfaceMesh.elementsBegin();
    typename SURFMESH::ElementPtrConstIter seEnd  = surfaceMesh.elementsEnd();
    for ( ; seIter != seEnd; ++seIter ) {

        // get pointer to surface element
        const typename SURFMESH::Element* surfEp = *seIter;

        // go through all level set points
        for ( unsigned el = 0; el < levelSet.size(); el++ ) {

            // get level set datum
            LevelSet ls = levelSet[el];

            // compute distance data
            base::cut::distanceToElement( surfEp, isSigned, ls,
                                          pointIdentityTolerance );

            // set level set
            levelSet[el] = ls;
        }

    }

    return;
}

#endif
