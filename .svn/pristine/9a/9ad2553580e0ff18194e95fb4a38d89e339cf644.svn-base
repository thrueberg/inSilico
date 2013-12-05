//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FieldBinder.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_fieldbinder_hpp
#define base_fieldbinder_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
// base/asmb includes
#include <base/asmb/FieldIterator.hpp>
#include <base/asmb/FieldTupleBinder.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        namespace detail_{

            //------------------------------------------------------------------
            /** Template placeholder and provider of dummy iterators.
             *  In order to allow for a large number of fields bound to the
             *  mesh, the non-used parameters are occupied by this object.
             *  For a smooth connection with base::asmb::FieldIterator it
             *  also provides dummy element iterators.
             */
            struct DummyField
            {
                detail_::DummyIterator elementsBegin() const
                {
                    return detail_::makeDummyIterator();
                }

                detail_::DummyIterator elementsEnd() const
                {
                    return detail_::makeDummyIterator();
                }

                detail_::DummyElementPtr elementPtr( const std::size_t n ) const
                {
                    return detail_::makeDummyElementPtr();
                }
                
            };

            //! Give a constant dummy field for the default parameters
            inline const DummyField makeDummyField() { return DummyField(); }
        }

        //----------------------------------------------------------------------
        template<typename MESH,
                 typename FIELD1,
                 typename FIELD2 = const detail_::DummyField,
                 typename FIELD3 = const detail_::DummyField,
                 typename FIELD4 = const detail_::DummyField,
                 typename FIELD5 = const detail_::DummyField>
        class FieldBinder;

        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            // Element pointer iterator types in dependence of a const qualifier
            template<typename FIELD>
            struct ElementPtrIterType
            {
                typedef typename FIELD::ElementPtrIter Type;
            };

            template<typename FIELD>
            struct ElementPtrIterType<const FIELD>
            {
                typedef typename FIELD::ElementPtrConstIter Type;
            };

            // The type is the dummy iterator
            template<>
            struct ElementPtrIterType<const DummyField>
            {
                typedef detail_::DummyIterator Type;
            };

        } // namespace detail_
        
    } // namespace asmb
} // namespace base


//------------------------------------------------------------------------------
/** Bind a mesh with up to five different fields.
 *
 *  \tparam MESH            Type of mesh to which a field is/fields are bound
 *  \tparam FIELD1          Type of field bound to the mesh
 *  \tparam FIELD2, FIELD3,
 *          FIELD4, FIELD5  Optional fields bound to the mesh
 */
template<typename MESH,   typename FIELD1,
         typename FIELD2, typename FIELD3,
         typename FIELD4, typename FIELD5>
class base::asmb::FieldBinder
    : public boost::noncopyable
{
public:
    //! Tuple of references to the fields
    typedef boost::tuple<FIELD1&, FIELD2&, FIELD3&, FIELD4&, FIELD5&> FieldTuple;

    //! The compound iterator over the bound fields
    typedef
    base::asmb::FieldIterator<typename detail_::ElementPtrIterType<MESH  >::Type,
                              typename detail_::ElementPtrIterType<FIELD1>::Type,
                              typename detail_::ElementPtrIterType<FIELD2>::Type,
                              typename detail_::ElementPtrIterType<FIELD3>::Type,
                              typename detail_::ElementPtrIterType<FIELD4>::Type,
                              typename detail_::ElementPtrIterType<FIELD5>::Type>
    FieldIterator;

    //! For instrospection
    typedef typename FieldIterator::ElementPtrTuple ElementPtrTuple;

    //--------------------------------------------------------------------------
    //! For convenience, delegate the type specification from here
    template<int I, int J=-1, int K=-1, int L=-1, int M=-1>
    struct TupleBinder
    {
        // Define the type of the tuple-binder
        typedef typename
        base::asmb::FieldTupleBinder<ElementPtrTuple,I,J,K,L,M> Type;
    };

    //! Constructor with mesh, a field and four optional fields
    FieldBinder( MESH& mesh,
                 FIELD1& field1,
                 FIELD2& field2 = detail_::makeDummyField(),
                 FIELD3& field3 = detail_::makeDummyField(),
                 FIELD4& field4 = detail_::makeDummyField(),
                 FIELD5& field5 = detail_::makeDummyField() )
        : mesh_( mesh ),
          fieldTuple_( field1, field2, field3, field4, field5 )
    { }

    //! Begin of elements of compound
    FieldIterator elementsBegin() const
    {
        return
            FieldIterator( mesh_.elementsBegin(),
                           (fieldTuple_.template get<0>()).elementsBegin(),
                           (fieldTuple_.template get<1>()).elementsBegin(),
                           (fieldTuple_.template get<2>()).elementsBegin(),
                           (fieldTuple_.template get<3>()).elementsBegin(),
                           (fieldTuple_.template get<4>()).elementsBegin() );
    }

    //! End of elements of compound
    FieldIterator elementsEnd() const
    {
        return
            FieldIterator( mesh_.elementsEnd(),
                           (fieldTuple_.template get<0>()).elementsEnd(),
                           (fieldTuple_.template get<1>()).elementsEnd(),
                           (fieldTuple_.template get<2>()).elementsEnd(),
                           (fieldTuple_.template get<3>()).elementsEnd(),
                           (fieldTuple_.template get<4>()).elementsEnd() );
    }

    //! Random-access to the pointer tuples
    ElementPtrTuple elementPtr( const std::size_t& e ) const
    {
        FieldIterator iter = this -> elementsBegin();
        std::advance( iter, e );
        return *iter;
    }

private:
    //! Mesh reference
    MESH&      mesh_;
    //! Tuple of field references
    FieldTuple fieldTuple_;
};

#endif
