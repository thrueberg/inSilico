//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   algorithms.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_aux_algorithms_hpp
#define base_aux_algorithms_hpp

namespace base{
    namespace aux{

        //! Apply op to two iterators' values simultaneously
        template<typename IT1, typename IT2, typename OP>
        OP forEach2( IT1 first1, IT1 last1, IT2 first2, OP op )
        {
            while( first1 != last1 ) {
                op( *first1, *first2 );
                ++first1;
                ++first2;
            }
            return op;
        }

        //! Apply op to three iterators' values simultaneously
        template<typename IT1, typename IT2, typename IT3, typename OP>
        OP forEach3( IT1 first1, IT1 last1, IT2 first2, IT3 first3, OP op )
        {
            while( first1 != last1 ) {
                op( *first1, *first2, *first3 );
                ++first1;
                ++first2;
                ++first3;
            }
            return op;
        }
        

    }
}

#endif

