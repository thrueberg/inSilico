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
#include <map>
#include <vector>
// boost includes
#include <boost/utility.hpp>
// base includes
#include <base/shape.hpp>
#include <base/types.hpp>
// base/aux
#include <base/aux/SortArray.hpp>
// base/mesh includes
#include <base/mesh/ElementFaces.hpp>
#include <base/mesh/FaceIterator.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        class MeshBoundary;

    }
}

//------------------------------------------------------------------------------
/** Extract and store from a range of elements the boundary of that range.
 *  The boundary of a range of elements is defined by the surfaces of these
 *  elements, which appear only once. In order to have access to the full
 *  geometry, pairs of element pointers and face numbers are stored.
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
        // Template parameter: type of element
        typedef typename base::TypeReduction<typename EITER::value_type>::Type Element;

        // Deduction of geometric entities
        static const base::Shape elemShape = Element::shape;
        static const base::NFace surface   = base::ShapeSurface<elemShape>::value;

        // Iterator over all faces of the element range
        typedef base::mesh::FaceIterator<EITER,surface> FaceIter;

        // Type of face
        typedef typename FaceIter::Face Face;

        // Temporary map for storage of the boundary faces
        typedef std::map< Face, std::pair<std::size_t, unsigned> > BoundaryMap;
        BoundaryMap boundaryMap;

        // Go through all faces of the range
        FaceIter faceIter = FaceIter( first );
        FaceIter faceEnd  = FaceIter( last  );

        for ( ; faceIter != faceEnd; ++faceIter ) {

            // Get face be dereferencing the iterator
            Face faceSorted = *faceIter;

            // Sort face array in order to have a unique key
            base::aux::SortArray<Face>::apply( faceSorted );

            // Try to find sorted face in the map
            typename BoundaryMap::iterator check = boundaryMap.find( faceSorted );

            // If not found, insert new pair of Element pointer and face number
            if ( check == boundaryMap.end() ) {
                const std::size_t elemID = (*(faceIter.elementIterator() )) -> getID();
                const unsigned nf = faceIter.faceNum();
                boundaryMap.insert( std::make_pair( faceSorted,
                                                    std::make_pair( elemID, nf ) ) );
            }
            else // Otherwise, face is interior to mesh and to be removed
                boundaryMap.erase( check );
        }

        // Pass element,faceNumber pairs to storage
        typename BoundaryMap::iterator bMapIter = boundaryMap.begin();
        typename BoundaryMap::iterator bMapEnd  = boundaryMap.end();
        for( ; bMapIter != bMapEnd; ++bMapIter )
            boundaryElements_.push_back( bMapIter -> second );
        
    }
    
public:
    //--------------------------------------------------------------------------
    /** Storage of element pointers along the boundary together with the
     *  corresponding number of the element face which coincides with the
     *  boundary.
     */
    typedef std::vector< std::pair<std::size_t,unsigned> > BoundaryElementContainer;

    //! @name Access to the storage of the boundary definition
    //@{
    typedef typename BoundaryElementContainer::const_iterator BoundConstIter;
    BoundConstIter boundaryBegin() const { return boundaryElements_.begin(); }
    BoundConstIter boundaryEnd()   const { return boundaryElements_.end();   }
    //@}
    
private:
    //! Local storage of the boundary defintion
    BoundaryElementContainer boundaryElements_;
};

#endif
