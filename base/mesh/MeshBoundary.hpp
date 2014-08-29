//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   MeshBoundary.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_meshboundary_hpp
#define base_mesh_meshboundary_hpp

//------------------------------------------------------------------------------
// std  includes
#include <utility>
#include <vector>
// boost includes
#include <boost/utility.hpp>
// base includes
#include <base/MultiIndex.hpp>
// base/mesh includes
#include <base/mesh/createBoundaryFromUnstructured.hpp>
#include <base/mesh/createBoundaryFromStructuredFace.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        class MeshBoundary;
    }
}

//------------------------------------------------------------------------------
/** Storage of the boundary of a mesh.
 *  Note that the boundary of a mesh is just a list of index pairs; the first
 *  number refers to the index of the element that lies adjacent to the boundary
 *  and the second number refers to the face number of this element which forms
 *  the boundary, see ElementFaces.
 *
 *  The action to create this list is delegated to either
 *  createBoundaryFromUnstructured() in case of an Unstructured mesh or to
 *  createBoundaryFromStructuredFace() in case of a structured mesh.
 *
 *  After creation of this list of index pairs, the function
 *  generateBoundaryMesh() allows to create a true surface mesh.
 *
 */
class base::mesh::MeshBoundary
    : boost::noncopyable
{
public:
    //--------------------------------------------------------------------------
    /** Create storage of (element pointer, face number) pairs of boundary.
     *  \param[in] first, last  Range of elements to extract the boundary from
     *  \tparam EITER Type of element iterator
     */
    template<typename EITER>
    void create( EITER first, EITER last )
    {
        base::mesh::createBoundaryFromUnstructured( first, last,
                                                    boundaryElements_ );
    }

    template<unsigned DIM>
    void create( const typename base::MultiIndex<DIM>::Type& gridSizes )
    {
        for ( unsigned d = 0; d < DIM; d ++ ) {
            this -> createFromFace<DIM>( gridSizes, d, 0 );
            this -> createFromFace<DIM>( gridSizes, d, 1 );
        }
    }

    template<unsigned DIM>
    void createFromFace( const typename base::MultiIndex<DIM>::Type& gridSizes,
                         const unsigned direction,
                         const unsigned side )
    {
        base::mesh::createBoundaryFromStructuredFace<DIM>( gridSizes,
                                                           direction,
                                                           side, 
                                                           boundaryElements_ );
    }

public:
    //--------------------------------------------------------------------------
    /** Storage of element IDs along the boundary together with the
     *  corresponding number of the element face which coincides with the
     *  boundary.
     */
    typedef std::vector< std::pair<std::size_t,unsigned> > BoundaryElementContainer;

    //! @name Access to the storage of the boundary definition
    //@{
    typedef BoundaryElementContainer::const_iterator BoundConstIter;
    BoundConstIter begin() const { return boundaryElements_.begin(); }
    BoundConstIter end()   const { return boundaryElements_.end();   }
    //@}
    
private:
    //! Local storage of the boundary defintion
    BoundaryElementContainer boundaryElements_;
};

#endif
