//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FieldTupleBinder.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_fieldtuplebinder_hpp
#define base_asmb_fieldtuplebinder_hpp

#include <base/types.hpp>
#include <base/verify.hpp>
#include <base/asmb/FieldElementPointerTuple.hpp>


//------------------------------------------------------------------------------

namespace base{
    namespace asmb{

        template<typename FIELDTUPLE,
                 int I=-1, int J=-1, int K=-1, int L=-1, int M=-1>
        struct FieldTupleBinder;

        namespace detail_{

            //! Compile-time assertion of the index values
            template<int I, int MAX>
            struct RangeCheck
            {
                STATIC_ASSERT_MSG( (I >= -1) and (I != 0) and (I <= MAX),
                                   "Index is out of range" );
            };
        }
        
        
    }
}

//------------------------------------------------------------------------------
/** Re-constructs tuples of element pointers according to given indices.
 *  In many applications, e.g. mixed methods, various solution fields are used
 *  in different orders. For instance in case of a standard discretisation of
 *  Stokes' system, the used bilinear forms expands as
 *  \f[
 *       a(u,v) + b(p,v) + b(q,u)
 *  \f]
 *  Ignoring boundary conditions, \f$ v \f$ is commonly taken from the same
 *  space as \f$ u \f$ and \f$ q \f$ from the space of \f$ p \f$. Therefore,
 *  the field composition \f$ ( u, p ) \f$ is used in various combinations
 *  where each constituent can serve as a test or a trial field.
 *  In order to avoid numerous redeclarations of these combinations, this class
 *  generates the desired combination of test, trial and auxiliary fields based
 *  on a given field composition.
 *  Assuming that all five indices are used, we would get the following picture
 *
 *   |    Index  |  Physical meaning                  |
 *   | --------- | ---------------------------------- |
 *   |      I    |  Test field                        |
 *   |      J    |  Trial field                       |
 *   |  K,L,M    |  1st, 2nd and 3rd auxiliary fields |
 *
 *
 *  Example
 *  -------
 *  For above example of Stokes' system, the following piece of code shall
 *  illustrate the functionality of this class.
 *  \code{.cpp}
 *  // Fields u and p of types U and P are defined
 *  base::asmb::FieldBinder<U,P> fieldBinder( u, p );
 *  // Legibility shortcut
 *  typedef base::asmb::FieldBinder<U,P>::ElementPtrTuple EPT;
 *  // generate the combinations for the three matrix blocks
 *  //
 *  //     ( A   B ) (u) = ...
 *  //     ( Bt  0 ) (p) = ...
 *  //
 *  typedef base::asmb::FieldTupleBinder<EPT,1,1> A;
 *  typedef base::asmb::FieldTupleBinder<EPT,1,2> B;
 *  typedef base::asmb::FieldTupleBinder<EPT,2,1> Bt;
 *  \endcode
 *  The types defined in the last 3 lines each provide a \a Tuple type and
 *  are themselves used to call the routine for the stiffness matrix
 *  computation.
 *  Note that up to 3 auxiliary fields can be bound into the tuple for multi-
 *  physical situations.
 *
 *  \tparam FIELDTUPLE Type of originally given field tuple, representing the
 *                     the field composition.
 *  \tparam I,J,K,L,M  Indices describing how a new field tuple shall be
 *                     constructed 
 */
template<typename FIELDTUPLE, int I, int J, int K, int L, int M>
struct base::asmb::FieldTupleBinder
{
    //! Type of field tuple to operate on
    typedef FIELDTUPLE FieldTuple;

    //! @name Sanity checks
    //@{
    detail_::RangeCheck<I,FieldTuple::maxNumFields> rc1;
    detail_::RangeCheck<J,FieldTuple::maxNumFields> rc2;
    detail_::RangeCheck<K,FieldTuple::maxNumFields> rc3;
    detail_::RangeCheck<L,FieldTuple::maxNumFields> rc4;
    detail_::RangeCheck<M,FieldTuple::maxNumFields> rc5;
    //@}

    //--------------------------------------------------------------------------
    //! @name Tuple construction
    //@{

    //! Type of the constructed tuple
    typedef
    base::asmb::FieldElementPointerTuple<typename FieldTuple::GeomElementPtr,
                                         typename FieldTuple::template Binder<I>::Type,
                                         typename FieldTuple::template Binder<J>::Type,
                                         typename FieldTuple::template Binder<K>::Type,
                                         typename FieldTuple::template Binder<L>::Type,
                                         typename FieldTuple::template Binder<M>::Type>
    Tuple;

    //! Generate a tuple from a given tuple based on the index ordering
    static Tuple makeTuple( const FieldTuple& ft )
    {
        return Tuple( ft.geomElementPtr(),
                      ft.template get<I>(),
                      ft.template get<J>(),
                      ft.template get<K>(),
                      ft.template get<L>(),
                      ft.template get<M>() );
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Transposed tuple construction
    //@{

    //! Type of the constructed tuple
    typedef
    base::asmb::FieldElementPointerTuple<typename FieldTuple::GeomElementPtr,
                                         typename FieldTuple::template Binder<J>::Type,
                                         typename FieldTuple::template Binder<I>::Type,
                                         typename FieldTuple::template Binder<K>::Type,
                                         typename FieldTuple::template Binder<L>::Type,
                                         typename FieldTuple::template Binder<M>::Type>
    TransposedTuple;

    //! Generate a tuple from a given tuple based on the index ordering
    static TransposedTuple makeTransposedTuple( const FieldTuple& ft )
    {
        return TransposedTuple( ft.geomElementPtr(),
                                ft.template get<J>(),
                                ft.template get<I>(),
                                ft.template get<K>(),
                                ft.template get<L>(),
                                ft.template get<M>() );
    }
    //@}

};

#endif
