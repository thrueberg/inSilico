//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FaceIterator.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_faceiterator_hpp
#define base_mesh_faceiterator_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/array.hpp>
#include <boost/iterator/iterator_facade.hpp>
// base includes
#include <base/shape.hpp>
#include <base/types.hpp>
// base/mesh includes
#include <base/mesh/ElementFaces.hpp>
#include <base/mesh/ExtractElementFaces.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename EITER, base::NFace NFACE>
        class FaceIterator;

        namespace detail_{
            //! Reduce the type name sizes
            template<typename EITER, base::NFace NFACE>
            struct FaceType
            {
                typedef typename base::TypeReduction<
                    typename EITER::value_type>::Type       Element;
                typedef typename ElementFaceTraits<Element::shape,
                                                   NFACE,
                                                   std::size_t>::Type Type;
            };

            //! Workaround for the vertex case
            template<typename ARRAY>
            void assignIndex( ARRAY & array, const unsigned num,
                              const std::size_t index )
            {
                array[ num ] = index;
            }

            void assignIndex( std::size_t & result, const unsigned dummy,
                              const std::size_t index )
            {
                result = index;
            }
        }

    }
}

//--------------------------------------------------------------------------
/** Iterator over the n-faces of a elements of a triangulation.
 *  Using this object, one can easily iterate over the vertices, edges,
 *  faces or cells of a mesh. The iterator holds an element iterator and
 *  a face counter, both of which are incremented appropriately.
 *
 *  \note There is no uniqueness in the iteration, i.e. if the same
 *        n-Face is contained by many elements it will be touched
 *        by this iterator each time an element contains it.
 *
 *  Currently, this iterator is implemented as a random access iterator. One
 *  could think of a forward iterator only which remembers the already
 *  considered faces and thus guarantees uniqueness.
 *
 *  This iterator fulfils all requirements of a 'Random access iterator'
 *  and is implemented by using Boost's iterator facade.
 *
 *  \tparam EITER  Underlying element iterator
 *  \tparam NFACE  Type of n-face to iterate over
 */
template<typename EITER, base::NFace NFACE>
class base::mesh::FaceIterator
    : public boost::iterator_facade<base::mesh::FaceIterator<EITER,
                                                             NFACE>,
                                    const typename
                                    base::mesh::detail_::FaceType<EITER,NFACE>::Type,
                                    boost::random_access_traversal_tag,
                                    const typename
                                    base::mesh::detail_::FaceType<EITER,NFACE>::Type>
{

public:
    //! @name template parameter
    //@{
    typedef EITER ElementIterator;
    static  const base::NFace nFace = NFACE;
    //@}

    //! For better legibility
    typedef FaceIterator<ElementIterator,nFace> SelfType;

    //! Helper based on Element iterator
    typedef base::mesh::detail_::FaceType<ElementIterator,nFace> FaceType;

    //! Element type
    typedef typename FaceType::Element Element;

    //! Type of face storage
    typedef typename FaceType::Type    Face;

    /** Boost's iterator facace
     *  Template parameters: - FaceIterator: Used forCRTP
     *                       - const Face:   The value type of this
     *                       - tag:          Iterator category
     *                       - const Face:   Reference type
     *  Note that the latter deviates from the default 'const Simplex&'
     *  in order to have the 'dereference' function const as required.
     */
    typedef boost::iterator_facade<SelfType,const Face,
                                   boost::random_access_traversal_tag,
                                   const Face> Facade;

    //! Element face traits
    typedef base::mesh::ElementFaceTraits<Element::shape,nFace> ElementFaceTraits;

    //! Lookup object
    typedef typename base::mesh::ElementFaces<Element::shape,nFace> ElementFaces;

    //! @name Number of vertices involved
    //@{
    static const unsigned numElemVertices =
        base::NumNFaces<Element::shape,base::VERTEX>::value;
    static const unsigned numFaceVertices = ElementFaceTraits::numVertices;
    //@}

    //! Number of faces
    static const unsigned numFaces = ElementFaces::numFaces;

    //--------------------------------------------------------------------------
    //! No default construction shall be allowed
    FaceIterator() { VERIFY_MSG( false, "Must not call empty constructor" ); }

    //! Constructor with an element iterator
    FaceIterator( ElementIterator elementIterator )
    : elementIterator_( elementIterator ),
      faceNum_( 0 )
    {
        //! empty
    }

    //! Copy constructor mandatory
    FaceIterator( const SelfType & si )
    : elementIterator_( si.elementIterator_ ),
      faceNum_( si.faceNum_ )
    {
        //! empty
    }

public:
    //--------------------------------------------------------------------------
    //! @name Accessor methods to private data
    //@{
    unsigned        faceNum()         const { return faceNum_; }
    ElementIterator elementIterator() const { return elementIterator_; }
    //@}
    
private:
    //--------------------------------------------------------------------------
    //! @name Implementation of iterator's core behaviour
    //@{

    //! Open back-door for the facade to acces these functions
    friend class boost::iterator_core_access;

    /** This iterator goes over all faces of an element and goes to the next
     *  element after the last face. Therefore, the incrementation has two
     *  cases (i) incrementation of the face counter and (ii) if the face
     *  counter reaches its maximum, reset it and increment the element
     *  reader.
     */
    void increment()
    {
        // Increment face number (i)
        faceNum_++;

        // If maximum is reached, reset face number and increment element (ii)
        if ( faceNum_ == numFaces ) {
            faceNum_ = 0;       // reset face number
            ++elementIterator_; // increment element iterator
        }
        
    }

    //! The decrementation inverses the logic of #increment()
    void decrement()
    {
        if ( faceNum_ == 0 ) {
            faceNum_ = numFaces-1;
            --elementIterator_;
        }
        else faceNum_--;
    }

    /** Two iterators are equal if
     *  (i)  the face numbers and
     *  (ii) the element iterators are equal
     */
    bool equal( const SelfType & other ) const
    {
        return ( (other.faceNum_ == faceNum_) and
                 (other.elementIterator_ == elementIterator_ ) );
    }

    //! Return a face. This requires to extract a face based on the current state.
    const Face dereference() const 
    {
        Face face;

        boost::array<std::size_t,numElemVertices> vertexIDs;
        
        base::mesh::ExtractElementVertices<Element>::apply( *elementIterator_,
                                                            vertexIDs );

        for ( unsigned v = 0; v < numFaceVertices; v ++ )
            ElementFaceTraits::assignIndex( face, v,
                                            vertexIDs[ ElementFaces::index( faceNum_, v ) ]
                );

        return face;
    }

    //! Move forward (backward) n steps
    void advance( const typename Facade::difference_type n )
    {
        // Deduce increments in element and face from n
        typename Facade::difference_type elemInc = n / numFaces;
        typename Facade::difference_type faceInc = n - elemInc * numFaces;

        // If face increment is too large, move to next element
        if ( (faceNum_ + faceInc) >= numFaces ) {
            elemInc++;
            faceInc -= numFaces;
        }
        else if ( (faceNum_ + faceInc) < 0 ) {
            // If the face increment is too small, move to previous element
            elemInc--;
            faceInc += numFaces;
        }

        // Change element iterator
        std::advance( elementIterator_, elemInc );

        // Change face number
        faceNum_ += numFaces;
    }

    //! Compute the difference between two iterators
    typename Facade::difference_type distance_to( const SelfType & other ) const
    {
        typename Facade::distance_type elemDiff =
            std::distance( elementIterator_, other.elementIterator_ );
        typename Facade::distance_type faceDiff =
            other.faceNum_ - faceNum_;

        return elemDiff * numFaces + faceDiff;
    }
    //@}

private:
    ElementIterator elementIterator_; //!< Iterator to mesh element
    unsigned        faceNum_;         //!< Number of current face of element
};


#endif
