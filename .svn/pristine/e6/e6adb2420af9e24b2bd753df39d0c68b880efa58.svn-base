//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   NearlyIncompNeoHookean.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef mat_hypel_nearlyincompneohookean_hpp
#define mat_hypel_nearlyincompneohookean_hpp

#include <mat/TensorAlgebra.hpp>


//------------------------------------------------------------------------------
namespace mat{
    namespace hypel{

        class NearlyIncompNeoHookean;
    }
}

//------------------------------------------------------------------------------
/** A nearly incompressible Neo-Hookean material.
 *  This hyperelastic material serves as test case in which the transition from
 *  nearly incompressible to incompressible is taken. Therefore, a conventional
 *  Neo-Hookean energy function for an incompressible material
 *  \f[
 *      \psi_{iso} ( F ) = \frac{\mu}{2} (F:F - 3)
 *  \f]
 *  and augmentented for the nearly incompressible case with the volume energy
 *  \f[
 *      \psi_{vol} ( F ) = \frac{\kappa}{2} ( \det F - 1 )^2
 *  \f]
 *  The parameters are the shear modulus \f$ \mu \f$ and the bulk modulus
 *  \f$ \kappa \f$.
 *  This object returns the iso-choric second Piola-Kirchhoff stress tensor
 *  \f$ S_{iso} \f$, the iso-choric elasticity tensor \f$ C_{iso} \f$, the
 *  inverse bulk modulus (in fact the 2nd derivative of the volume energy) and
 *  the bulk ratio, which is the quotient of first and second derivatives of
 *  the volume energy.
 *  A flag isIncompressible_ is stored in order to turn on and off the full
 *  incompressibility, i.e. \f$ \kappa \to \infty \f$.
 */
class mat::hypel::NearlyIncompNeoHookean
{
public:
    
    //! @name Types of tensors used
    //@{
    typedef mat::Tensor                     Tensor;
    typedef mat::ElastTensor                ElastTensor;
    //@}
    
    //! Construct with Lame parameters
    NearlyIncompNeoHookean( const double bulk, const double mu,
                            const bool isIncompressible = false )
        : bulk_( bulk ), mu_( mu ),
          isIncompressible_( isIncompressible )
    { }

    //--------------------------------------------------------------------------
    /** Strain energy density function.
     *  \f[
     *    \psi(F) = \psi_{vol} + \psi_{iso}
     *            = \frac{\kappa}{2} (J - 1)^2
     *            + \frac{\mu}{2}( tr C - 3 ), 
     *            J  = \det F, C = F^T F
     *  \f]
     *  In case of full incompressibilty, only the iso-choric part is returned.
     *  \param[in]  F  Deformation gradient
     *  \return        Elastic energy
     */
    double energy( const Tensor& F ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        // compute trace of which
        const double trC = mat::trace( C );

        // compute determinant of F and its logarithm
        const double J = mat::determinant( F );

        // complete energy
        const double psi = mu_/2. * (trC - 3.) +
            ( isIncompressible_ ? 0. : (bulk_/2.)*(J-1.)*(J-1.) );

        // return energy
        return psi;
    }

    //--------------------------------------------------------------------------
    /** The iso-choric part of the second Piola-Kirchhoff stress tensor.
     *  The definition \f$ S^{iso} = \partial \psi_{iso} / \partial E \f$ yields
     *  \f[
     *         S^{iso} = \mu J^{-2/3} ( I - \frac{1}{3} tr(C) C^{-1} )
     *  \f]
     *  \param[in]  F  Deformation gradient
     *  \param[out] S  2nd Piola-Kirchhoff stress tensor (iso-choric part)
     */
    void secondPiolaKirchhoff( const Tensor& F, Tensor& S ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        // compute trace of which
        const double trC = mat::trace( C );

        // compute inverse of Cauchy-Green tensor
        Tensor Cinv;
        mat::inverse( C, Cinv );

        // compute determinant of F and its logarithm
        const double J = mat::determinant( F );

        const double factor = mu_ * (std::pow( J, -2./3. ) );

        // evaluate the second Piola-Kirchhoff stress tensor
        S.noalias() = factor* (Tensor::Identity() - trC/3. * Cinv);
    }

    //--------------------------------------------------------------------------
    /** The iso-choric elasticity tensor in material description.
     *  \f[
     *     C^{iso}_{ABCD} = 2 \mu J^{-2/3} \left[
     *        \frac{tr(C)}{3} \tilde{I}_{ABCD}
     *      - \frac{1}{3} \delta_{AB} C^{-1}_{CD}
     *      - \frac{1}{3} C^{-1}_{AB} \delta_{CD}
     *      + \frac{tr(C)}{9} C^{-1}_{AB} C^{-1}_{CD} \right]
     *  \f]
     *  where the tensor \f$ \tilde{I} \f$ is defined as
     *  \f[
     *     \tilde{I}_{ABCD} = \frac{1}{2} \left[
     *           C^{-1}_{AC} C^{-1}_{BD} + C^{-1}_{AD} C^{-1}_{BC}
     *                         \right]
     *  \f]
     *  Note: The storage is based on the minor symmetries.
     *  \param[in]  F    Deformation gradient
     *  \param[out] elC  Elasticity tensor
     */
    void materialElasticityTensor( const Tensor& F, ElastTensor& elC ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor CG;
        mat::rightCauchyGreen( F, CG );

        // compute trace of which
        const double trC = mat::trace( CG );

        // compute inverse of Cauchy-Green tensor
        Tensor Cinv;
        mat::inverse( CG, Cinv );

        // compute determinant of F 
        const double J = mat::determinant( F );

        // factor 
        const double fac = 2 * mu_ * std::pow( J, -2./3. );

        // go through all indices
        for ( unsigned A = 0; A < 3; A++ ) {
            for ( unsigned B = A; B < 3; B++ ) {

                // Voigt-Index I
                const unsigned voigtI = mat::Voigt::apply( A, B );
                
                for ( unsigned C = 0; C < 3; C++ ) {
                    for ( unsigned D = C; D < 3; D++ ) {

                        // Voigt-Index J 
                        const unsigned voigtJ = mat::Voigt::apply( C, D );

                        // I_{ABCD} = - d C^{-1}_{AB}/d C_{CD}
                        const double iABCD =
                            0.5 * ( Cinv(A,C) * Cinv(B,D) + Cinv(A,D) * Cinv(B,C) );

                        // delta_{AB}, delta_{CD}
                        const double deltaAB = (A==B? 1. : 0.);
                        const double deltaCD = (C==D? 1. : 0.);

                        // compute C_{ABCD}
                        const double cEntry = fac *
                            ( trC/3. * iABCD -
                              1.0/3. * deltaAB   * Cinv(C,D) -
                              1.0/3. * Cinv(A,B) * deltaCD   +
                              trC/9. * Cinv(A,B) * Cinv(C,D) );
                        
                        // insert to elasticity tensor
                        elC( voigtI, voigtJ ) = cEntry;
                        
                    } // D
                } // C
                
            } // B
        } // A
        
    }

    //--------------------------------------------------------------------------
    /** Ratio of first to second derivative of the volume energy.
     *  Returns
     *  \f[
     *             \frac{ \partial   \psi_{vol} }{\partial J  }
     *      \left( \frac{ \partial^2 \psi_{vol} }{\partial J^2} \right)^{-1}
     *  \f]
     *  which for the considered energy function \f$ \psi_{vol} \f$ becomes
     *  \f[
     *       \frac{ \kappa (J - 1) }{ \kappa } = J - 1
     *  \f]
     *  \param[in] F  Deformation gradient
     *  \return       Ratio of derivatives
     */
    double giveBulkRatio( const Tensor& F ) const
    {
        // compute determinant of F 
        const double J = mat::determinant( F );
        return ( J - 1. );
    }

    //--------------------------------------------------------------------------
    /** Give the inverse of the second derivative of the volume energy.
     *  Returns
     *  \f[
     *      \left( \frac{\partial^2 \psi_{vol} }{\partial J^2} \right)^{-1} =
     *         \frac{1}{ \kappa }
     *  \f]
     *  In case of incompressibility, the value zero is returned.
     *  \param[in] F   Deformation gradient  (unused in this case)
     *  \return        Inverse of bulk modulus
     */
    double giveInverseBulk( const Tensor& F) const
    {
        return ( isIncompressible_ ? 0. : 1./bulk_ );
    }
    
private:
    //! @name Material parameters
    //@{
    const double bulk_;             //!< Bulk modulus \f$ \kappa \f$
    const double mu_;               //!< Shear modulus
    const bool   isIncompressible_; //!< Flag to turn on full incompressibility
    //@}
};

#endif
