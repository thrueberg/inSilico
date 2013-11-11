#ifndef terzaghi_h
#define terzaghi_h

#include <cmath>


#include <mat/Lame.hpp>

//------------------------------------------------------------------------------
/** Terzaghi's one-dimensional consolidation solution.
 *  The rotated situation looks as follows
 *  \code{.txt}
 *           
 *              F    p=0, t_x = -F H(t)               u=0        x
 *             ----> +================================+       ---> 
 *   (surface)      x=0                              x=H     (depth)
 *
 *  \endcode
 *  For this problem analytic series solutions of the pressure field
 *  \f$ p(x,t) \f$ and the vertical displacement \f$ u(x,t) \f$ exist and are
 *  implemented in this class.
 */
class Terzaghi
{
private:
    //--------------------------------------------------------------------------
    /** Compute the undrained Poisson ratio.
     *  \f[
     *       \nu_u = \frac{3 K_u - 2 \mu}{2 (3 K_u + \mu)}
     *  \f]
     *  which is the same formula as computing the Poisson ratio \f$ \nu \f$
     *  from the bulk modulus \f$ K \f$ and the shear modulus \f$ \mu \f$, but
     *  with the bulk modulus replaced by the undrained bulk modulus
     *  \f[
     *         K_u = K + \frac{ \alpha^2 }{ c_0 }
     *  \f]
     */
    static double coeffNuU_( const double E, const double nu,
                             const double alpha, const double c0 )
    {
        // shear modulus mu
        const double mu = mat::Lame::mu( E, nu );
        // buld modulus K
        const double K  = mat::Lame::bulk( E, nu );
        // undrained bulk modulus Ku
        const double KU = K + (alpha*alpha) / c0;
        // undrained Poisson ratio
        const double nuU = (3. * KU - 2. * mu) / (2. * (3. * KU + mu) );

        return nuU;
    }

    //--------------------------------------------------------------------------
    /** Compute some exotic storage coefficient.
     *  \f[
     *       S = c_0 \frac{ (1 - \nu_u)(1 - 2\nu) }{ (1 - \nu)(1 - 2\nu_u)}
     *  \f]
     *
     */
    static double coeffS_( const double E, const double nu,
                           const double alpha, const double c0 )
    {
        // undrained Poisson ratio
        const double nuU = coeffNuU_( E, nu, alpha, c0 );
        
        // coefficient S
        const double S =
            c0 * (1. - nuU) * (1. - 2.*nu) / (1. - nu) / (1. - 2.*nuU);

        return S;
    }

    //--------------------------------------------------------------------------
    /** Compute yet another coefficient.
     *  \f[
     *       \eta = \frac{\alpha}{2} \frac{1- 2\nu}{1 - \nu}
     *  \f]
     */
    static double coeffEta_( const double alpha, const double nu )
    {
        const double eta = (alpha / 2.) * (1 - 2. * nu) / (1. - nu);
        return eta;
    }

    
public:
    //--------------------------------------------------------------------------
    /** Constructor with basic poro-elastic material parameters.
     *  \param[in]  E       Young's modulus
     *  \param[in]  nu      Poisson ratio
     *  \param[in]  alpha   Biot coefficient
     *  \param[in]  c0      Storage coefficient
     *  \param[in]  k       Mobility coefficient (permeability / viscosity)
     *  \param[in]  H       Height of soil layer
     *  \param[in]  F       Applied surface traction
     */
    Terzaghi( const double E,  const double nu, const double alpha,
              const double c0, const double  k, const double H,
              const double F,  const unsigned numTerms = 25 )
        : mu_(  mat::Lame::mu( E, nu ) ),
          nu_(  nu ),
          nuU_( coeffNuU_( E, nu, alpha, c0 ) ), 
          S_(   coeffS_(E, nu, alpha, c0) ),
          eta_( coeffEta_(alpha, nu) ),
          c_(   k / S_ ),
          H_( H ), F_( F ),
          numTerms_( numTerms )
    {
        // -- Range checks --
        // - E and nu are already checked in mat::Lame
        VERIFY_MSG( (alpha > 0.) and (alpha <= 1.0),
                    "Value of alpha violates bounds" );

        VERIFY_MSG( (c0 > 0. ),
                    "Terzaghi's solution not valid for zero storage coeff. " );

        VERIFY_MSG( (k > 0.) and (H > 0.) and (F > 0.),
                    "Values have to be positive" );
    };

    //--------------------------------------------------------------------------
    /** Initial pore pressure field.
     *  At time \f$ t = 0^+ \f$, the instantaneous reaction is the initial pore
     *  pressure given as
     *  \f[
     *      p_0 = \frac{\eta F}{\mu S}
     *  \f]
     */
    double initialPressure( ) const
    {
        return (eta_ * F_)/(mu_ * S_);
    }

    //--------------------------------------------------------------------------
    /** Pressure response at given time and depth.
     *  \f[
     *        p(x,t) = p_0 [1 - F_1( (x/H), (ct/4/L^2) )]
     *  \f]
     */
    double pressure( const double x, const double t )
    {
        const double chi = x / H_;
        const double tau = c_ * t / 4. / H_ / H_;
        return initialPressure() * (1. - F1_( chi, tau ) );
    }

    //--------------------------------------------------------------------------
    /** Instantaneous displacement at given depth.
     *  
     */
    double initialDisplacement( const double x )
    {
        const double chi = x / H_;
        return F_ * ( H_ * (1. - 2.*nuU_)) / (2. * mu_ * (1. - nuU_)) * (1. - chi);
    }

    //--------------------------------------------------------------------------
    /** Additional, time-dependent displacement at given depth.
     *
     */
    double incrementalDisplacement( const double x, const double t )
    {
        const double chi = x / H_;
        const double tau = c_ * t / 4. / H_ / H_;

        return
            F_ * ( H_ * (nuU_ - nu_) )/( 2.*mu_ * (1. - nu_)*(1. - nuU_)) *
            F2_( chi, tau );
    }

    //--------------------------------------------------------------------------
    /** Total vertical displacement at given depth and time.
     */
    double displacement( const double x, const double t )
    {
        return initialDisplacement( x ) + incrementalDisplacement( x, t );
    }

    
private:
    //--------------------------------------------------------------------------
    //! @name Series expansions
    //@{
    double F1_( const double chi, const double tau )
    {
        double f1 = 1.;

        for ( unsigned i = 0; i < numTerms_; i++ ) {
            const double m = static_cast<double>(2*i+1);
            f1 -= (4. / m / M_PI) *
                std::sin( m * M_PI * chi / 2. ) *
                std::exp( -m*m * M_PI*M_PI * tau );
        }
        
        return f1;
    }

    double F2_( const double chi, const double tau )
    {
        double f2 = 0.;

        for ( unsigned i = 0; i < numTerms_; i++ ) {
            const double m = static_cast<double>(2*i+1);

            f2 += (8. / (m*m * M_PI*M_PI)) *
                std::cos( m * M_PI * chi / 2.) *
                (1. - std::exp( -m*m * M_PI*M_PI * tau));
        }

        return f2;
    }
    //@}

private:
    const double mu_;  //!< Shear modulus
    const double nu_;  //!< Poisson ratio
    const double nuU_; //!< Undrained Poisson ratio
    const double S_;   //!< Storage coefficient
    const double eta_; //!< Poroelastic stress coefficient
    const double c_;   //!< Diffusivity coefficient (gen. consolidation coeff.)
    const double H_;   //!< Height of soil layer
    const double F_;   //!< Applied surface traction at t = 0+

    const unsigned numTerms_; //!< Number of terms in series expansion
};





#endif
