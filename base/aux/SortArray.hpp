//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SortArray.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_aux_sortarray_hpp
#define base_aux_sortarray_hpp

//------------------------------------------------------------------------------
// std includes
#include <algorithm>
//
#include <boost/type_traits.hpp>
//------------------------------------------------------------------------------
namespace base{
    namespace aux{

        namespace detail_{
            
            template<typename NOTANARRAY>
            struct DoNotDoAnything
            {
                static void apply( NOTANARRAY & notAnArray ) {  }
            };

            template<typename ARRAY>
            struct SortAnArray
            {
                static void apply( ARRAY & array )
                {
                    std::sort( array.begin(), array.end() );
                }
            };
        }

        //----------------------------------------------------------------------
        /** Sort a given array of numbers if it really is an array.
         *  In the implementation faces are stored as arrays of vertex numbers.
         *  In order to have unique identification of faces, a sorted version
         *  of such an array is used as identification key (as in maps or sets).
         *  But, for being generic, a 0-face is a vertex and stored as a single
         *  number only. The application of a sort algorithm to that number is
         *  not feasible. For this reason, the following struct inherits from
         *  a sorting struct or a do-nothing struct in dependence of the
         *  condition if the given type is really an array.
         *
         *  \note This implementation checks it the type is a fundamental
         *        type (int, float, etc.). This dose \em not imply that
         *        in the contrary case, the type would be an array.
         *        Unfortunately, boost's is_array only works for standard
         *        C-arrays, but not its own array class.
         *        Better solution to be found
         *
         *  \tparam MAYBEANARRAY  Type which might be an array
         */
        template<typename MAYBEANARRAY>
        struct SortArray
            : public base::IfElse< boost::is_fundamental<MAYBEANARRAY>::value,
                                   detail_::DoNotDoAnything<MAYBEANARRAY>,
                                   detail_::SortAnArray<MAYBEANARRAY> >::Type
        { };

    }
}
#endif
