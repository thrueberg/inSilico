//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   compareNumbers.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_auxi_comparenumbers_hpp
#define base_auxi_comparenumbers_hpp

// std includes
#include <limits>
// base includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace auxi{

        //----------------------------------------------------------------------
        /** Tolerance comparison of two doubles.
         *  This version of comparison is taken from Ericson 'Real-Time
         *  Collision Detection' and tries to avoid all pitfalls.
         *  If not provided by the caller, the tolerance is set to
         *  the square-root of the machine epsilon
         *  (Difference between 1 and the least value greater than 1 that is
         *   representable).
         *  The comparison for \f$ a == b \f$ becomes
         *  \f[
         *      |a-b| <= eps (|a| + |b| + 1)
         *  \f]
         *                                      
         *  \param[in] a,b  The values to compare
         *  \param[in] eps  Tolerance
         *  \return    True if values are equal within the tolerance
         */
        bool almostEqualNumbers( const double a, const double b,
                                 const double eps =
                                 std::sqrt( std::numeric_limits<double>::epsilon() ) )
        {
            return
                std::abs( a-b ) <=
                eps * (std::abs(a) + std::abs(b) + 1.0 );
        }

        //----------------------------------------------------------------------
        /** Tolerance comparison of two vectors of size DIM.
         *  Compares all components of the given vectors with
         *  almostEqualNumbers.
         *  \param[in] a,b  Vectors to compare
         *  \param[in] eps  Tolerance
         *  \return    True if vectors are equal within the tolerance
         */
        template<unsigned DIM>
        bool almostEqualVectors( const typename base::Vector<DIM>::Type& a,
                                 const typename base::Vector<DIM>::Type& b,
                                 const double eps =
                                 std::sqrt( std::numeric_limits<double>::epsilon() ) )
        {
            for ( unsigned d = 0; d < DIM; d++ ) 
                if ( not almostEqualNumbers( a[d], b[d], eps ) ) return false;

            return true;
        }
        
    }
}
        

#endif
