//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   NeoHookean.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef surf_neohookean_hpp
#define surf_neohookean_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/verify.hpp>

//------------------------------------------------------------------------------
namespace surf{
    
    class NeoHookean;
}

//------------------------------------------------------------------------------
/** Hyperelastic membrane energy a la neo-Hooke.
 *
 */
class surf::NeoHookean
{
public:
    
    //! Construct with Lame parameters
    NeoHookean( const double A )
        : A_( A )
    {
        VERIFY_MSG( (A_ >  0.0 ), "Paramter has to be positive" );
    }


public:
    //--------------------------------------------------------------------------
    /** The hyperelastic energy 
     *  \f[
     *      W = \frac{A}{4} \left( I_1 -1 + \frac{1}{I_2 + 1} \right)
     *  \f]
     *  \param[in] I1,I2 Strain invariants
     *  \return       Value of elastic energy
     */
    double energy( const double I1, const double I2 ) const
    {
        // give value of energy
        return A_/4. * ( I1 -1. + 1. /(I2 + 1.) );
    }

public:
    
    //! @name Partial derivatives of energy w.r.t invariants
    //@{
    double dWdI1( const double I1, const double I2 ) const
    {
        return A_/2.;
    }

    double dWdI2( const double I1, const double I2 ) const
    {
        return -A_/2. * 1./( (I2+1.) * (I2+1.) );
    }

    double d2WdI1dI1( const double I1, const double I2 ) const
    {
        return 0.;
    }
    
    double d2WdI2dI2( const double I1, const double I2 ) const
    {
        return A_ * 1./( (I2+1.) * (I2+1.) * (I2+1.) );
    }

    double d2WdI1dI2( const double I1, const double I2 ) const
    {
        return 0.;
    }
    //@}

private:
    //! Material parameter
    const double A_;
};

#endif
