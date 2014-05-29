//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   NeoHookeanCompressible.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef mat_hypel_neohookeancompressible_hpp
#define mat_hypel_neohookeancompressible_hpp

#include <mat/TensorAlgebra.hpp>

//------------------------------------------------------------------------------
namespace mat{
    namespace hypel{

        class NeoHookeanCompressible;
    }
}

//------------------------------------------------------------------------------
/** A compressible Neo-Hookean material.
 */
class mat::hypel::NeoHookeanCompressible
{
public:
    
    //! @name Types of tensors used
    //@{
    typedef mat::Tensor                     Tensor;
    typedef mat::ElastTensor                ElastTensor;
    //@}
    
    //! Construct with Lame parameters
    NeoHookeanCompressible( const double lambda, const double mu )
        : lambda_( lambda ), mu_(     mu )
    { }

    //--------------------------------------------------------------------------
    /** Strain energy density function.
     *  \f[
     *    \psi(F) = \frac{\lambda}{2} (\log J)^2 - \mu \log J
     *            + \frac{\mu}{2}( tr C - 3 ), 
     *            J  = \det F, C = F^T F
     *  \f]
     *  \param[in]  F  Deformation gradient
     *  \returns       Elastic energy
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
        const double logJ = std::log( J );

        // return energy
        return lambda_/2. * (logJ * logJ) - mu_ * logJ + mu_/2. * (trC - 3.);
    }

    //--------------------------------------------------------------------------
    /** The second Piola-Kirchhoff stress tensor S.
     *  The definition \f$ S = \partial \psi / \partial E \f$ yields
     *  \f[
     *         S = (\lambda \log J - \mu) C^{-1} + \mu I
     *  \f]
     *  \param[in]  F  Deformation gradient
     *  \param[out] S  2nd Piola-Kirchhoff stress tensor
     */
    void secondPiolaKirchhoff( const Tensor& F, Tensor& S ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        // compute inverse of Cauchy-Green tensor
        Tensor Cinv;
        mat::inverse( C, Cinv );

        // compute determinant of F and its logarithm
        const double J = mat::determinant( F );
        const double logJ = std::log( J );

        // evaluate the second Piola-Kirchhoff stress tensor
        S.noalias() = (lambda_ * logJ - mu_) * Cinv + mu_ * Tensor::Identity();
    }

    //--------------------------------------------------------------------------
    /** The Cauch stress tensor sigma.
     *  The definition \f$ \sigma = \partial \psi / \partial e \f$ yields
     *  \f[
     *         \sigma = \frac{\mu}{J} (b - I) + \frac{\lambda}{J}(\log J) I
     *  \f]
     *  \param[in]  F      Deformation gradient
     *  \param[out] sigma  Cauchy stress tensor
     */
    void cauchyStress( const Tensor& F, Tensor& sigma ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor b;
        mat::leftCauchyGreen( F, b );

        // compute determinant of F and its logarithm
        const double J = mat::determinant( F );
        const double logJ = std::log( J );

        // evaluate the second Piola-Kirchhoff stress tensor
        sigma.noalias() = mu_/J * b + (lambda_ * logJ - mu_)/J * Tensor::Identity();
    }

    //--------------------------------------------------------------------------
    /** The elasticity tensor in material description.
     *  \f[
     *       C_{ABCD} = \lambda C^{-1}_{AB} C^{-1}_{CD} +
     *          (\mu-\log J) (C^{-1}_{AC} C^{-1}_{BD} + C^{-1}_{AD} C^{-1}_{BC})
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

        // compute inverse of Cauchy-Green tensor
        Tensor Cinv;
        mat::inverse( CG, Cinv );

        // compute determinant of F 
        const double J = mat::determinant( F );

        // factor of second term
        const double fac2 = mu_ - lambda_ * std::log( J );

        // go through all indices
        for ( unsigned A = 0; A < 3; A++ ) {
            for ( unsigned B = A; B < 3; B++ ) {

                // Voigt-Index 1
                const unsigned AB = mat::Voigt::apply( A, B );
                
                for ( unsigned C = 0; C < 3; C++ ) {
                    for ( unsigned D = C; D < 3; D++ ) {

                        // compute C_{ABCD}
                        const double cEntry =
                            lambda_ * Cinv(A,B) * Cinv(C,D) +
                            fac2 * ( Cinv(A,C) * Cinv(B,D) + Cinv(A,D) * Cinv(B,C) );
                        
                        // Voigt-Index 2
                        const unsigned CD = mat::Voigt::apply( C, D );

                        // insert to elasticity tensor
                        elC( AB, CD ) = cEntry;
                        
                    } // D
                } // C
                
            } // B
        } // A

        return;
    }

    //--------------------------------------------------------------------------
    /** The elasticity tensor in spatial description.
     *  \f[
     *       C_{ijkl} =
     *         \frac{\lambda}{J}             \delta_{ij} \delta_{kl} +
     *         \frac{\mu-\lambda \log J}{J} (\delta_{ik} \delta_{jl} +
     *                                       \delta_{il} \delta_{jk} )
     *  \f]
     *  Note: The storage is based on the minor symmetries.
     *  \param[in]  F    Deformation gradient
     *  \param[out] elC  Elasticity tensor
     */
    void spatialElasticityTensor( const Tensor& F, ElastTensor& elC ) const
    {
        // compute determinant of F 
        const double J = mat::determinant( F );

        // factors
        const double fac1 =  lambda_ / J;
        const double fac2 = ( mu_ - lambda_ * std::log( J ) ) / J;

        // go through all indices
        for ( unsigned i = 0; i < 3; i++ ) {
            for ( unsigned j = i; j < 3; j++ ) {

                // Voigt-Index 1
                const unsigned ij = mat::Voigt::apply( i, j );
                
                for ( unsigned k = 0; k < 3; k++ ) {
                    for ( unsigned el = k; el < 3; el++ ) {

                        // compute C_{ijkl}
                        const double cEntry =
                            fac1 *   (i == j ? 1. : 0.) * (k ==el ? 1. : 0.) +
                            fac2 * ( (i == k ? 1. : 0.) * (j ==el ? 1. : 0.) +
                                     (i ==el ? 1. : 0.) * (j == k ? 1. : 0.) );
                        
                        // Voigt-Index 2 
                        const unsigned kl = mat::Voigt::apply( k, el );

                        // insert to elasticity tensor
                        elC( ij, kl ) = cEntry;
                        
                    } // el
                } // k
                
            } // j
        } // i

        return;
    }

    
private:
    //! @name Lame parameters
    //@{
    const double lambda_;
    const double mu_;
    //@}
};

#endif
