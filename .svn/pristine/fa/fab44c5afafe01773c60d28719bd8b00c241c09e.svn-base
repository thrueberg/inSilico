//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   numbers.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_numbers_hpp
#define base_numbers_hpp

//------------------------------------------------------------------------------
// std includes
#include <limits>
#include <complex>

//------------------------------------------------------------------------------
namespace base{

    //! The largest available unsigned is used to represent an invalid number
    static const unsigned invalidInt = static_cast<unsigned>( -1 );

    //! The largest representable real number
    inline double invalidReal()
    {
        return  std::numeric_limits<double>::max();
    }

    //! Number type for degrees of freedom
#ifdef INSILICO_COMPLEX
    typedef std::complex<double>  number;

    //! The largest representable number
    static number invalidNumber()
    {
        return number( std::numeric_limits<double>::max() );
    }

#else
    
    typedef double        number;

    //! The largest representable number
    inline number invalidNumber()
    {
        return  invalidReal();
    }
    
#endif


}

//------------------------------------------------------------------------------
#endif
