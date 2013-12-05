//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   HyperElasticMembrane.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef surf_hyperelasticmembrane_hpp
#define surf_hyperelasticmembrane_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/verify.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>

//------------------------------------------------------------------------------
namespace surf{
    
    template<typename ENERGY>
    class HyperElasticMembrane;
}

//------------------------------------------------------------------------------
/** Hyperelastic membrane material.
 */
template<typename ENERGY>
class surf::HyperElasticMembrane
{
public:
    //! Type of energy providing derivatives w.r.t to invariants
    typedef ENERGY Energy;
    
    //! @name Types of tensors used
    //@{
    typedef mat::Tensor                     Tensor;
    typedef mat::ElastTensor                ElastTensor;
    //@}

    //! Constructor with access to energy
    HyperElasticMembrane( const Energy& energy )
        : energy_( energy )
    { }

public:

    //--------------------------------------------------------------------------
    //! @name First and second strain invariant
    //@{

    /** First invariant is
     *  \f[
     *       I_1 = tr( F^T F ) - 2 = tr(C) - 2
     *  \f]
     *  Note that in derivatives, it acts like the standard first invariant
     *  \f$ tr(C) \f$.
     */
    double firstInvariant( const Tensor& C ) const
    {
        return mat::trace( C) - 2.;
    }

    /** Second invariant is
     *  \f[
     *       I_2 = det( F^T F ) - 1 = \det(C) - 2
     *  \f]
     *  Note that in derivatives, it acts like the standard third invariant
     *  \f$ \det(C) \f$.
     */
    double secondInvariant( const Tensor& C ) const
    {
        return mat::determinant( C ) - 1.;
    }
    //@}
    

public:

    //--------------------------------------------------------------------------
    /** The second Piola-Kirchhoff stress tensor S.
     *  The definition \f$ S = \partial \psi / \partial E \f$ yields
     *  \f[
     *    S = 2 \left( \frac{\partial W}{\partial I_1} I +
     *                 \frac{\partial W}{\partial I_2} \det{C} C^{-1} \right)
     *  \f]
     *  \param[in]  F  Deformation gradient
     *  \param[out] S  2nd Piola-Kirchhoff stress tensor
     */
    void secondPiolaKirchhoff( const Tensor& F, Tensor& S ) const
    {
        // compute Cauchy-Green stretch tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        Tensor Cinv;
        const double detC = mat::inverseAndDet( C, Cinv );
        
        // invariants
        const double I1 = firstInvariant( C );
        const double I2 = secondInvariant( C );

        // energy derivatives
        const double dWdI1 = energy_.dWdI1( I1, I2 );
        const double dWdI2 = energy_.dWdI2( I1, I2 );

        // evaluate the second Piola-Kirchhoff stress tensor
        S.noalias() =
            2. * dWdI1 * Tensor::Identity() + 2. * dWdI2 * detC * Cinv;
    }

    //--------------------------------------------------------------------------
    /** The elasticity tensor in material description.
     *  Using expressions (6.193/4) from Holzapfel (2000), identifying the
     *  invariant \f$ I_2 \f$ from here with his \f$ I_3 \f$, and discarding all
     *  terms in that book with derivatives w.r.t \f$ I_2 \f$ gives finally
     *  \f[
     *       C_{ABCD} =
     *          \delta_1 \delta_{AB} \delta_{CD} +
     *          \delta_3 (\delta_{AB} C^{-1}_{CD} + C^{-1}_{AB} \delta_{CD})
     *          \delta_6 C^{-1}_{AB} C^{-1}_{CD} +
     *    \frac{\delta_7}{2} (C^{-1}_{AC} C^{-1}_{BD} + C^{-1}_{AD} C^{-1}_{BC} )
     *  \f]
     *  where the factors are
     *  \f[
     *     \delta_1 = 4 \frac{\partial^2 W}{\partial I_1^2} \qquad
     *     \delta_3 = 4 \det{C} \frac{\partial^2 W}{\partial I_1 \partial I_2}
     *     \qquad
     *     \delta_6 = 4 \left( \det(C) \frac{\partial   W}{\partial I_2} +
     *                       \det(C)^2 \frac{\partial^2 W}{\partial I_2^2} \right)
     *     \qquad
     *     \delta_7 = -4 \det{C} \frac{\partial W}{\partial I_2}
     *  \f]
     *
     *  Note: The storage is based on the minor symmetries.
     *  \param[in]  F  Deformation gradient
     *  \param[out] C  Elasticity tensor
     */
    void materialElasticityTensor( const Tensor& F, ElastTensor& CE ) const
    {
        // compute Cauchy-Green stretch tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        Tensor Cinv;
        const double detC = mat::inverseAndDet( C, Cinv );

        // invariants
        const double I1 = firstInvariant( C );
        const double I2 = secondInvariant( C );

        // energy derivatives
        const double dWdI2 = energy_.dWdI2( I1, I2 );

        const double d2WdI1dI1 = energy_.d2WdI1dI1( I1, I2 );
        const double d2WdI2dI2 = energy_.d2WdI2dI2( I1, I2 );
        const double d2WdI1dI2 = energy_.d2WdI1dI2( I1, I2 );

        // constant factors (non-zero ones)
        const double d1 =  4. * d2WdI1dI1;
        const double d3 =  4. * detC * d2WdI1dI2;
        const double d6 =  4. * ( detC * dWdI2 + detC*detC * d2WdI2dI2 );
        const double d7 = -4. * detC * dWdI2;
        
        // clear to zero 
        CE = ElastTensor::Constant( 0. );

        // go through all 81 indices
        for ( unsigned I = 0; I < 3; I++ ) {
            for ( unsigned J = 0; J <=I ; J++ ) {
                const unsigned v1 = mat::Voigt::apply( I, J );
                
                for ( unsigned K = 0; K < 3; K++ ) {
                    for ( unsigned L = 0; L <=K ; L++ ) {

                        const unsigned v2 = mat::Voigt::apply( K, L );
                        
                        CE( v1, v2 ) +=
                            d1    * (I==J and K==L ? 1. : 0.) +
                            d3    * (I==J ? Cinv(K,L) : 0.) +
                            d3    * (K==L ? Cinv(I,J) : 0.) + 
                            d6    *  Cinv(I,J) * Cinv(K,L) +
                            d7/2. * (Cinv(I,K) * Cinv(J,L) + Cinv(I,L) * Cinv(J,K));
                    }
                }
            }
        }

        return;
    }
    
    
private:
    //! Energy
    const Energy energy_;
};

#endif
