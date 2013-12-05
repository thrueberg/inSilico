//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FieldIterator.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_fielditerator_hpp
#define base_asmb_fielditerator_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/tuple/tuple.hpp>
#include <boost/iterator/iterator_facade.hpp>
// base/asmb/ includes
#include <base/asmb/FieldElementPointerTuple.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        //----------------------------------------------------------------------
        namespace detail_ {
            
            //! Dummy object as template placeholder for iterators
            struct DummyIterator
                : public boost::iterator_facade<DummyIterator,
                                                detail_::DummyElementPtr,
                                                boost::random_access_traversal_tag,
                                                detail_::DummyElementPtr>
            {
                //! Boost's facade will provide the proper behaviour
                typedef boost::iterator_facade<DummyIterator,
                                               const detail_::DummyElementPtr,
                                               boost::random_access_traversal_tag,
                                               const detail_::DummyElementPtr> Facade;
                
            private:
                //! Open back-door for the facade to acces these functions
                friend class boost::iterator_core_access;

                void increment() {}
                void decrement() {}
                void advance( const Facade::difference_type n ) {}

                Facade::difference_type distance_to( const DummyIterator & other ) const
                {
                    return 0;
                }

                bool equal( const DummyIterator & other ) const
                {
                    return true;
                }
                
                Facade::reference dereference() const
                {
                    return detail_::makeDummyElementPtr();
                }
            };

            //! Helper for construction of a dummy interator
            inline DummyIterator makeDummyIterator()
            { return DummyIterator(); }
            
        } // namespace detail_

        template<typename GITER,
                 typename FITER1,
                 typename FITER2 = detail_::DummyIterator,
                 typename FITER3 = detail_::DummyIterator,
                 typename FITER4 = detail_::DummyIterator,
                 typename FITER5 = detail_::DummyIterator>
        class FieldIterator;

        
    } // namespace asmb
} // namespace base

//------------------------------------------------------------------------------
/** Iterator for combining a mesh with fields.
 *  When assembling system matrices and vectors, it is required to iterate over
 *  the elements of the mesh and a certain number of fields. This iterator
 *  is a random access iterator based on boost's iterator_facade. Note that
 *  contrary to boost's zip_iterator, this iterator dereferences to a tuple of
 *  the value types of the iterators, i.e. pointers to elements. This can only
 *  be achieved by changing the iterators reference type to a value type.
 *
 *  \tparam GITER           Iterator over the mesh elements (mandatory)
 *  \tparam FITER1          Iterator over first  field      (mandatory)
 *  \tparam FITER2, FITER3,
 *          FITER4, FITER5  Optional iterators over other fields.
 */
template<typename GITER,
         typename FITER1, typename FITER2, typename FITER3,
         typename FITER4, typename FITER5>
class base::asmb::FieldIterator
    : public boost::iterator_facade< base::asmb::FieldIterator<GITER,  FITER1,
                                                               FITER2, FITER3, 
                                                               FITER4, FITER5>,
                                     const typename
                                     base::asmb::detail_::
                                     TupleTypeBinder<GITER, FITER1,
                                                     FITER2,FITER3,
                                                     FITER4, FITER5>::Type,
                                     boost::random_access_traversal_tag,
                                     const typename // important: overload of
                                     base::asmb::detail_:: // reference type !!
                                     TupleTypeBinder<GITER, FITER1,
                                                     FITER2,FITER3,
                                                     FITER4, FITER5>::Type>
{
public:
    //! Return type of the dereference operation: a tuple of pointers
    typedef typename base::asmb::detail_::TupleTypeBinder<GITER, FITER1,
                                                          FITER2,FITER3,
                                                          FITER4, FITER5>::Type
    ElementPtrTuple;

    //! Tuple of iterators
    typedef boost::tuple<GITER,FITER1,FITER2,FITER3,FITER4,FITER5> IterTuple;

    //! For better legibility
    typedef base::asmb::FieldIterator<GITER,  FITER1,
                                      FITER2, FITER3, 
                                      FITER4, FITER5> SelfType;

    //! Boost's facade will provide the proper behaviour
    typedef boost::iterator_facade<SelfType,
                                   const ElementPtrTuple,
                                   boost::random_access_traversal_tag,
                                   const ElementPtrTuple> Facade;

    //! Construction
    FieldIterator( GITER  gIter,
                   FITER1 fIter1,
                   FITER2 fIter2 = FITER2(),
                   FITER3 fIter3 = FITER3(),
                   FITER4 fIter4 = FITER4(),
                   FITER5 fIter5 = FITER5() )
        : iterTuple_( gIter, fIter1, fIter2, fIter3, fIter4, fIter5 )
    {
        // emtpy
    }

    //! Copy constructor mandatory
    FieldIterator( const SelfType & other )
    : iterTuple_( other.iterTuple_ )
    {
        // empty
    }
   
private:
    //--------------------------------------------------------------------------
    //! @name Implementation of iterator's core behaviour
    //@{

    //! Open back-door for the facade to acces these functions
    friend class boost::iterator_core_access;

    /** Increment every member iterator
     */
    void increment()
    {
        ++(iterTuple_.template get<0>());
        ++(iterTuple_.template get<1>());
        ++(iterTuple_.template get<2>());
        ++(iterTuple_.template get<3>());
        ++(iterTuple_.template get<4>());
        ++(iterTuple_.template get<5>());
    }

    /** Decrement every memeber iterator
     */
    void decrement()
    {
        --(iterTuple_.template get<0>());
        --(iterTuple_.template get<1>());
        --(iterTuple_.template get<2>());
        --(iterTuple_.template get<3>());
        --(iterTuple_.template get<4>());
        --(iterTuple_.template get<5>());
    }

    //! Move forward (backward) n steps
    void advance( const typename Facade::difference_type n )
    {
        std::advance(  iterTuple_.template get<0>(), n);
        std::advance(  iterTuple_.template get<1>(), n);
        std::advance(  iterTuple_.template get<2>(), n);
        std::advance(  iterTuple_.template get<3>(), n);
        std::advance(  iterTuple_.template get<4>(), n);
        std::advance(  iterTuple_.template get<5>(), n);
    }

    //! Compute the difference between two iterators
    typename Facade::difference_type distance_to( const SelfType & other ) const
    {
        // Distance between the mesh iterators is sufficient
        return std::distance( iterTuple_.template       get<0>(),
                              other.iterTuple_.template get<0>() );
    }

    /** Two iterators are equal if the tuples they hold are equal
     *  Sufficient to compare the mesh iterators
     */
    bool equal( const SelfType & other ) const
    {
        return ( iterTuple_.template       get<0>( ) ==
                 other.iterTuple_.template get<0>( ) );
    }

    //! Return a tuple containing the values of the iterators
    const  ElementPtrTuple dereference() const 
    {
        ElementPtrTuple ept( *(iterTuple_.template get<0>() ),
                             *(iterTuple_.template get<1>() ),
                             *(iterTuple_.template get<2>() ),
                             *(iterTuple_.template get<3>() ),
                             *(iterTuple_.template get<4>() ),
                             *(iterTuple_.template get<5>() ) );
        return ept;
    }
    //@}

private:
    IterTuple iterTuple_; //!< Holds a tuple of iterators
};

#endif
