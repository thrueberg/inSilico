//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SurfaceFieldBinder.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_surfacefieldbinder_hpp
#define base_surfacefieldbinder_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>
// base/asmb includes
#include <base/asmb/FieldElementPointerTuple.hpp>
#include <base/asmb/FieldBinder.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{


        //----------------------------------------------------------------------
        template<typename SURFMESH,
                 typename FIELD1,
                 typename FIELD2 = const detail_::DummyField,
                 typename FIELD3 = const detail_::DummyField,
                 typename FIELD4 = const detail_::DummyField,
                 typename FIELD5 = const detail_::DummyField>
        class SurfaceFieldBinder;

        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            // Element pointer iterator types in dependence of a const qualifier
            template<typename FIELD>
            struct ElementPtrType
            {
                typedef typename FIELD::Element* Type;
            };

            // The type is the dummy iterator
            template<>
            struct ElementPtrType<const DummyField>
            {
                typedef detail_::DummyElementPtr Type;
            };

        } // namespace detail_
        
    } // namespace asmb
} // namespace base


//------------------------------------------------------------------------------
/** Bind a mesh with up to five different fields.
 *
 *  \tparam SURFMESH        Type of surface mesh to which a field is/fields are bound
 *  \tparam FIELD1          Type of field bound to the mesh
 *  \tparam FIELD2, FIELD3,
 *          FIELD4, FIELD5  Optional fields bound to the mesh
 */
template<typename SURFMESH, typename FIELD1,
         typename FIELD2,   typename FIELD3,
         typename FIELD4,   typename FIELD5>
class base::asmb::SurfaceFieldBinder
    : public boost::noncopyable
{
public:

    //! Directly define the element pointer tuple
    typedef FieldElementPointerTuple<typename SURFMESH::Element*,
                                     typename detail_::ElementPtrType<FIELD1>::Type,
                                     typename detail_::ElementPtrType<FIELD2>::Type,
                                     typename detail_::ElementPtrType<FIELD3>::Type,
                                     typename detail_::ElementPtrType<FIELD4>::Type,
                                     typename detail_::ElementPtrType<FIELD5>::Type>
    ElementPtrTuple;

    //! Iterator of local container
    typedef typename std::vector<ElementPtrTuple>::const_iterator FieldIterator;
    

    //! Constructor with mesh, a field and four optional fields
    SurfaceFieldBinder( SURFMESH& surfMesh,
                        FIELD1& field1,
                        FIELD2& field2 = detail_::makeDummyField(),
                        FIELD3& field3 = detail_::makeDummyField(),
                        FIELD4& field4 = detail_::makeDummyField(),
                        FIELD5& field5 = detail_::makeDummyField() )
    {
        // iterate over the surface mesh
        for ( typename SURFMESH::ElementPtrConstIter bElemIter =
                  surfMesh.elementsBegin();
              bElemIter != surfMesh.elementsEnd(); ++bElemIter ) {

            const std::size_t elemID = (*bElemIter) -> getID();

            ElementPtrTuple ept =
                ElementPtrTuple( *bElemIter,
                                 field1.elementPtr( elemID ),
                                 field2.elementPtr( elemID ),
                                 field3.elementPtr( elemID ),
                                 field4.elementPtr( elemID ),
                                 field5.elementPtr( elemID ) );
                                 

            elementPtrs_.push_back( ept );
        }

    }

    //! Begin of elements of compound
    FieldIterator elementsBegin() const
    {
        return elementPtrs_.begin();
    }

    //! End of elements of compound
    FieldIterator elementsEnd() const
    {
        return elementPtrs_.end();
    }

private:
    //! Container of element pointer tuples
    std::vector<ElementPtrTuple> elementPtrs_;

};

#endif
