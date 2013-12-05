//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   EqualPointers.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_auxi_equalpointers_hpp
#define base_auxi_equalpointers_hpp

//------------------------------------------------------------------------------
namespace base{
    namespace auxi{

        //! To pointers of different types must not be equal
        template<typename A, typename B>
        struct EqualPointers
        {
            static bool apply( const A* a, const B* b )
            {
                return false;
            }
        };

        //! Compare the addresses of two pointers of same type
        template<typename A>
        struct EqualPointers<A,A>
        {
            static bool apply( const A* a, const A* b )
            {
                return a==b;
            }
        };
        
    }
}

#endif
