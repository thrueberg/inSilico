//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Skalak.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef surf_skalak_hpp
#define surf_skalak_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/verify.hpp>

//------------------------------------------------------------------------------
namespace surf{
    
    class Skalak;
}

//------------------------------------------------------------------------------
/** Hyperelastic membrane material due to Skalak.
 *
 *  Reference: Skalak et al., Biophysical Journal, 1973.
 */
class surf::Skalak
{
public:
    
    //! Construct with Lame parameters
    Skalak( const double A, const double B )
        : A_( A ), B_( B ) 
    {
        VERIFY_MSG( (A_ >  0.0 ), "First paramter has to be positive" );
        VERIFY_MSG( (B_ > -0.5 ), "Second parameter has to be larger than 1/2" );
    }


public:
    //--------------------------------------------------------------------------
    /** The hyperelastic energy due to Skalak reads
     *  \f[
     *      W = \frac{A}{4} \left( I_1^2 + 2 I_1 - 2 I_2 + B I_2^2 \right)
     *  \f]
     *  \param[in] F Deformation gradient
     *  \return      Value of elastic energy
     */
    double energy( const double I1, const double I2 ) const
    {
        // give value of energy
        return A_/4. * ( I1*I1 + 2. * I1 - 2. * I2 + B_ * I2*I2 );
    }

public:
    
    //! @name Partial derivatives of energy w.r.t invariants
    //@{
    double dWdI1( const double I1, const double I2 ) const
    {
        return A_/2. * (I1 + 1.);
    }

    double dWdI2( const double I1, const double I2 ) const
    {
        return A_/2. * (B_ * I2 - 2.);
    }

    double d2WdI1dI1( const double I1, const double I2 ) const
    {
        return A_/2.;
    }
    
    double d2WdI2dI2( const double I1, const double I2 ) const
    {
        return B_ * A_/2.;
    }

    double d2WdI1dI2( const double I1, const double I2 ) const
    {
        return 0.;
    }
    //@}

private:
    //! @name Material parameters
    //@{
    const double A_;
    const double B_;
    //@}
};

#endif
