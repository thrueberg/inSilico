//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Lame.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef mat_lame_hpp
#define mat_lame_hpp

//------------------------------------------------------------------------------
// std  includes
#include <utility>
// base includes
#include <base/verify.hpp>
//------------------------------------------------------------------------------

namespace mat{

    //! First Lame parameter
    inline double lambda( const double E,
                          const double nu )
    {
        VERIFY_MSG( E > 0., "Negative Young's modulus" );
        VERIFY_MSG( (nu > -1.) and (nu < 0.5), "Poisson's ratio out of bounds" );
        return E * nu / (1. + nu) / (1. - 2.*nu);
    }

    //! Second Lame parameter
    inline double mu( const double E,
                      const double nu )
    {
        VERIFY_MSG( E > 0., "Negative Young's modulus" );
        VERIFY_MSG( (nu > -1.), "Poisson's ratio out of bound" );
        return E / 2. / (1. + nu );
    }

    //! Combine both Lame parameters
    inline std::pair<double,double> lame( const double E,
                                          const double nu )
    {
        return std::make_pair( lambda( E, nu ), mu( E, nu ) );
    }

    //! Compute bulk modulus
    inline double bulk( const double E,
                        const double nu )
    {
        VERIFY_MSG( E > 0., "Negative Young's modulus" );
        VERIFY_MSG( (nu > -1.) and (nu < 0.5), "Poisson's ratio out of bounds" );
        return E / 3. / (1. - 2.*nu);
    }
}

#endif
