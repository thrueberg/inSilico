//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FieldElementPointerTuple.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_fieldelementpointertuple_hpp
#define base_asmb_fieldelementpointertuple_hpp

//------------------------------------------------------------------------------
// std includes
#include <iterator>
// boost includes
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_same.hpp>
// base includes
#include <base/types.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        //----------------------------------------------------------------------
        namespace detail_{

            //! A dummy structure for placeholder purposes
            struct DummyElementPtr{};

            //! Function constructing the dummy
            inline DummyElementPtr makeDummyElementPtr()
            { return DummyElementPtr(); }
            
        }

        //----------------------------------------------------------------------
        template<typename GEOMEPTR,
                 typename FIELD1EPTR = detail_::DummyElementPtr,
                 typename FIELD2EPTR = detail_::DummyElementPtr,
                 typename FIELD3EPTR = detail_::DummyElementPtr,
                 typename FIELD4EPTR = detail_::DummyElementPtr,
                 typename FIELD5EPTR = detail_::DummyElementPtr>
        class FieldElementPointerTuple;

        //----------------------------------------------------------------------
        
        //! Generate a tuple from a surface field tuple
        template<typename SFEPT>
        struct DomainFieldElementPointerTuple
        {
            typedef FieldElementPointerTuple<
                typename
                base::TypeReduction<typename
                                    SFEPT::GeomElementPtr>::Type::DomainElement*,
                typename SFEPT::TestElementPtr,
                typename SFEPT::TrialElementPtr,
                typename SFEPT::AuxField1ElementPtr,
                typename SFEPT::AuxField2ElementPtr,
                typename SFEPT::AuxField3ElementPtr>
            Type;

            static Type convert( const SFEPT & sfept )
            {
                return Type( sfept.geomElementPtr() -> getDomainElementPointer(),
                             sfept.template get<1>(),
                             sfept.template get<2>(),
                             sfept.template get<3>(),
                             sfept.template get<4>(),
                             sfept.template get<5>() );
            }
        };


        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            // Helper to bind the type of FieldElementPointerTuple given 
            // iterators as types
            template<typename GITER,  typename FITER1,
                     typename FITER2, typename FITER3,
                     typename FITER4, typename FITER5>
            struct TupleTypeBinder
            {
                typedef
                base::asmb::FieldElementPointerTuple<
                    typename std::iterator_traits<GITER >::value_type,
                    typename std::iterator_traits<FITER1>::value_type,
                    typename std::iterator_traits<FITER2>::value_type,
                    typename std::iterator_traits<FITER3>::value_type,
                    typename std::iterator_traits<FITER4>::value_type,
                    typename std::iterator_traits<FITER5>::value_type>
                Type;
            };

            //------------------------------------------------------------------
            // Number of active fields
            template<typename FIELD1EPTR,
                     typename FIELD2EPTR, typename FIELD3EPTR,
                     typename FIELD4EPTR, typename FIELD5EPTR>
            struct NumFields
            {
                static const unsigned value =
                    ( boost::is_same<FIELD1EPTR,DummyElementPtr>::value ? 0 :
                      ( boost::is_same<FIELD2EPTR,DummyElementPtr>::value ? 1 :
                        ( boost::is_same<FIELD3EPTR,DummyElementPtr>::value ? 2 :
                          ( boost::is_same<FIELD4EPTR,DummyElementPtr>::value ? 3 :
                            ( boost::is_same<FIELD5EPTR,DummyElementPtr>::value ? 4 :
                              5 ) ) ) ) );
            };

            //------------------------------------------------------------------
            template<typename TUPLE,int N>
            struct TupleElementType
            {
                typedef typename boost::tuples::element<N,TUPLE>::type Type;
            };

            template<typename TUPLE>
            struct TupleElementType<TUPLE,-1>
            {
                typedef DummyElementPtr Type;
            };

            template<typename TUPLE,int N>
            struct GetTupleElement
            {
                static typename TupleElementType<TUPLE,N>::Type apply( const TUPLE& t )
                {
                    return t.template get<N>();
                }
            };

            template<typename TUPLE>
            struct GetTupleElement<TUPLE,-1>
            {
                static DummyElementPtr apply( const TUPLE& t )
                {
                    return makeDummyElementPtr();
                }
            };


        } // namespace detail_

    } // namespace asmb
} // namespace base


//------------------------------------------------------------------------------
/** A tuple of element  pointers.
 *  On the lowest level of iterating over a field compound, an object of this
 *  class is returned which allows access to the individual element pointers.
 *  This object is just a wrapper around boost::tuple for the two reason:
 *  1.  Future implementations might go to a std::tuple or other structs
 *  2.  For the caller this wrapper provides human-readable types and accessors
 *
 *  \tparam GEOMEPTR     Geometry element pointer
 *  \tparam FIELD1EPTR   Type of pointer to an element of the first field
 *  ...
 */
template<typename GEOMEPTR,   typename FIELD1EPTR,
         typename FIELD2EPTR, typename FIELD3EPTR,
         typename FIELD4EPTR, typename FIELD5EPTR>
class base::asmb::FieldElementPointerTuple
{
    //! The basic data type of object: a tuple of pointers
    typedef boost::tuple<GEOMEPTR,  FIELD1EPTR,FIELD2EPTR,
                         FIELD3EPTR,FIELD4EPTR,FIELD5EPTR> Tuple;
    
public:
    //! Maximal storage of fields as defined by the implementation
    static const unsigned maxNumFields = 5; 
    
    //! Number of fields which are not a dummy
    static const unsigned numFields = detail_::NumFields<FIELD1EPTR,FIELD2EPTR,
                                                         FIELD3EPTR,FIELD4EPTR,
                                                         FIELD5EPTR>::value;
    
    //! Type via index
    template<int N>
    struct Binder
    {
        typedef typename detail_::TupleElementType<Tuple,N>::Type Type;
    };

    //! @name Human readable names of the template parameter
    //@{
    typedef GEOMEPTR   GeomElementPtr;
    typedef FIELD1EPTR TestElementPtr;
    typedef FIELD2EPTR TrialElementPtr;
    typedef FIELD3EPTR AuxField1ElementPtr;
    typedef FIELD4EPTR AuxField2ElementPtr;
    typedef FIELD5EPTR AuxField3ElementPtr;
    //@}

    //! @name Extraction of the element types
    //@{
    typedef typename base::TypeReduction<GeomElementPtr     >::Type  GeomElement;
    typedef typename base::TypeReduction<TestElementPtr     >::Type  TestElement;
    typedef typename base::TypeReduction<TrialElementPtr    >::Type  TrialElement;
    typedef typename base::TypeReduction<AuxField1ElementPtr>::Type  AuxField1Element;
    typedef typename base::TypeReduction<AuxField2ElementPtr>::Type  AuxField2Element;
    typedef typename base::TypeReduction<AuxField3ElementPtr>::Type  AuxField3Element;
    //@}

    
    //! Constructor with all pointers
    FieldElementPointerTuple(
        GeomElementPtr      geomElementPtr,     
        TestElementPtr      testElementPtr      = detail_::makeDummyElementPtr(),
        TrialElementPtr     trialElementPtr     = detail_::makeDummyElementPtr(),
        AuxField1ElementPtr auxField1ElementPtr = detail_::makeDummyElementPtr(),
        AuxField2ElementPtr auxField2ElementPtr = detail_::makeDummyElementPtr(),
        AuxField3ElementPtr auxField3ElementPtr = detail_::makeDummyElementPtr() )
        : tuple_( geomElementPtr,      testElementPtr,      trialElementPtr,
                  auxField1ElementPtr, auxField2ElementPtr, auxField3ElementPtr )
    {
        // empty
    }
    
    //! Value access via index
    template<int N>
    const typename Binder<N>::Type get() const
    {
        return detail_::GetTupleElement<Tuple,N>::apply( tuple_ );
    }

    //! @name Access with names
    //@{
    GeomElementPtr           geomElementPtr() const { return this ->template get<0>(); }
    TestElementPtr           testElementPtr() const { return this ->template get<1>(); }
    TrialElementPtr         trialElementPtr() const { return this ->template get<2>(); }
    AuxField1ElementPtr auxField1ElementPtr() const { return this ->template get<3>(); }
    AuxField2ElementPtr auxField2ElementPtr() const { return this ->template get<4>(); }
    AuxField3ElementPtr auxField3ElementPtr() const { return this ->template get<5>(); }
    //@}

    //! @name Transpose test and trial fields
    //@{
    typedef FieldElementPointerTuple<GeomElementPtr,
                                     TrialElementPtr,TestElementPtr, // transposed
                                     AuxField1ElementPtr,
                                     AuxField2ElementPtr,
                                     AuxField3ElementPtr> TransposedTuple;
    TransposedTuple transpose() const
    {
        return TransposedTuple( tuple_.template get<0>(),
                                tuple_.template get<2>(), tuple_.template get<1>(), // T
                                tuple_.template get<3>(),
                                tuple_.template get<4>(),
                                tuple_.template get<5>() );
    }
    //@}

        
private:
    //! A tuple of element pointers (and dummy place holders)
    Tuple tuple_;
};


#endif
