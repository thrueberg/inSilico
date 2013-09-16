//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   analyticLevelSet.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_analyticalevelset_hpp
#define base_cut_analyticalevelset_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// base/cut includes
#include <base/cut/LevelSet.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename DOMAINMESH, typename SURFACE>
        void analyticLevelSet( const DOMAINMESH& domainMesh,
                               const SURFACE&    surface, 
                               const bool        isSigned, 
                               std::vector< base::cut::LevelSet<DOMAINMESH::Node::dim> >&
                               levelSet );

        //! Define type of an analytic surface representation
        template<unsigned DIM>
        struct AnalyticSurface
        {
            typedef boost::function<bool( const typename base::Vector<DIM>::Type&,
                                          typename base::Vector<DIM>::Type& )> Type;
        };
         
    }
}

//------------------------------------------------------------------------------
/** Compute the level set data for a given analytic surface expression.
 *  By passing the functor \em surface  of type \em SURFACE with the signature
 *  \code{.cpp}
 *  bool surface( const base::Vector<DIM>::Type& x,
 *                base::Vector<DIM>::Type& xClosest )
 *  \endcode
 *  (as defined by base::cut::AnalyticSurface ), the level set data are directly
 *  determined at every point. The functor is given \p x (the location of the
 *  level set datum) and  \p xClosest (the point on surface closest to \p x and
 *  modified as reference). It returns a flag if \p x is inside or outside of
 *  the surface.
 *  \tparam DOMAINMESH The domain mesh which gives the nodal coordinates
 *  \tparam SURFACE    Type of the functor describing the surface implicitly
 *  \param[in]  domainMesh The domain mesh
 *  \param[in]  surface    The surface representation
 *  \param[in]  isSigned   Flag if the level set is a signed distance function
 *  \param[out] levelSet   Vector of level set data, filled on return
 */
template<typename DOMAINMESH, typename SURFACE>
void base::cut::analyticLevelSet( const DOMAINMESH& domainMesh,
                                  const SURFACE&    surface, 
                                  const bool        isSigned, 
                                  std::vector< base::cut::LevelSet<DOMAINMESH::Node::dim> >&
                                  levelSet )
{
    // convenience typedef
    typedef base::cut::LevelSet<DOMAINMESH::Node::dim> LevelSet;
    
    // Fill level set data with node coordinates
    levelSet.resize( 0 );
    typename DOMAINMESH::NodePtrConstIter nIter = domainMesh.nodesBegin();
    typename DOMAINMESH::NodePtrConstIter nEnd  = domainMesh.nodesEnd();
    for ( ; nIter != nEnd; ++nIter ) {

        typename LevelSet::VecDim x, xClosest;
        (*nIter) -> getX( &(x[0]) );

        // store coordinate
        LevelSet ls( x );

        if ( isSigned ) { 
            // get closest point on the surface and side flag
            const bool isInterior = surface( x, xClosest );
            // set flag
            if ( isInterior ) ls.setInterior();
            else              ls.setExterior();
        }
        else {
            // get closest point on surface
            surface( x, xClosest );
        }

        // set closest point
        ls.setClosestPoint( xClosest );
            
        // store level set datum
        levelSet.push_back( ls );
    }

    return;
}

#endif
