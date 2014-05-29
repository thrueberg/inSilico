//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Ogden.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef mat_hypel_ogden_hpp
#define mat_hypel_ogden_hpp

// boost includes
#include <boost/array.hpp>
#include <boost/math/special_functions/cbrt.hpp>
// material includes
#include <mat/TensorAlgebra.hpp>

//------------------------------------------------------------------------------
namespace mat{
    namespace hypel{

        template<unsigned N>
        class Ogden;
    }
}

//------------------------------------------------------------------------------
/** A (nearly) incompressible Ogden material.
 *  \tparam N Number of parameter pairs in iso-choric energy.
 */
template<unsigned N>
class mat::hypel::Ogden
{
public:
    //! Template paramter
    static const unsigned nParam = N;
    
    //! @name Types of tensors used
    //@{
    typedef mat::Tensor                     Tensor;
    typedef mat::ElastTensor                ElastTensor;
    //@}

    //! Array of material parameters
    typedef boost::array<double,nParam>     ParamArray;
    
    //! Construct incompressible material
    Ogden( const ParamArray& mu,
           const ParamArray& alpha )
        : mu_( mu), alpha_( alpha ),
          kappa_( base::invalidReal() ),
          beta_(  base::invalidReal() ),
          isIncompressible_( true ),
          tol_( 1.e-8 )
    { }

    //! Construct slightly compressible material
    Ogden( const ParamArray& mu,
           const ParamArray& alpha,
           const double kappa,
           const double beta )
        : mu_( mu), alpha_( alpha ),
          kappa_( kappa ), beta_( beta ),
          isIncompressible_( false ),
          tol_( 1.e-8 )
    { }
    
   //--------------------------------------------------------------------------
    /** Strain energy density function.
     *  \f{eqnarray*}{
     *   \psi(F)       &=& \psi_{vol} + \psi_{iso} \cr
     *   \psi_{vol}(F) &=& \kappa \beta^{-2}(\beta \log J + J^{-\beta} - 1) \cr
     *   \psi_{iso}(F) &=& \sum_{p=1}^N \frac{\mu_p}{\alpha_p}
     *            \left( \bar{\lambda}_1^{\alpha_p} +
     *                   \bar{\lambda}_2^{\alpha_p} +
     *                   \bar{\lambda}_3^{\alpha_p} -3 \right) \cr
     *     J  &=& \det F, \quad C = F^T F, \quad \lambda_i^2 = \lambda( C ),
     *        \quad \bar{\lambda}_i = J^{-1/3} \lambda_i
     *  \f}
     *  In case of full incompressibilty, only the iso-choric part is returned.
     *  \param[in]  F  Deformation gradient
     *  \return        Elastic energy
     */
    double energy( const Tensor& F ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        // get eigenvalues of C
        const Vector evalC = base::eigenValues( C );

        // compute determinant of F and its logarithm
        const double J = mat::determinant( F );

        // rescale the stretches
        double scale = ( isIncompressible_ ? 1.0 : boost::math::cbrt( 1./ J ) );
        scale = scale * scale;

        // iso-choric energy
        double psiIso = 0.;
        for ( unsigned p = 0; p < nParam; p++ ) {
            psiIso += (mu_[p] / alpha_[p]) *
                ( std::pow( scale * evalC[0], 0.5 * alpha_[p] ) +
                  std::pow( scale * evalC[1], 0.5 * alpha_[p] ) +
                  std::pow( scale * evalC[2], 0.5 * alpha_[p] ) - 3. );
        }

        // volume energy
        double psiVol = 0.;
        if ( not isIncompressible_ ) {

            psiVol = (kappa_/beta_/beta_) *
                (beta_ * std::log( J ) + std::pow( J, -beta_) - 1.);
        }

        // return energy
        return psiIso + psiVol;
    }

    //--------------------------------------------------------------------------
    /** The iso-choric part of the second Piola-Kirchhoff stress tensor.
     *  The definition \f$ S^{iso} = \partial \psi^{iso} / \partial E \f$ yields
     *  \f[
     *      S^{iso} = \sum_{a=1}^3 S^{iso}_a N_a \otimes N_a
     *  \f]
     *  with the principal stresses
     *  \f[
     *    S^{iso}_a = \frac{1}{\lambda_a}
     *                \frac{\partial \psi^{iso}}{\partial \lambda_a}
     *              = \frac{1}{\lambda_a^2} \left(
     *                       \bar{\lambda_a} \psi^{iso}_{,\bar{a}}
     *     - \frac{1}{3} \sum_b \bar{\lambda_b} \psi^{iso}_{,\bar{b}}
     *               \right)
     *  \f]
     *  and the principal directions \f$ N_a \f$. The partial derivatives are
     *  computed as
     *  \f[
     *       \psi^{iso}_{,\bar{a}} =
     *             \frac{\partial \psi^{iso}}{\partial \bar{\lambda}_a}
     *        = \sum_{p=1}^N \mu_p \bar{\lambda}_a^{\alpha_p -1}
     *  \f]
     *  \param[in]  F  Deformation gradient
     *  \param[out] S  2nd Piola-Kirchhoff stress tensor (iso-choric part)
     */
    void secondPiolaKirchhoff( const Tensor& F, Tensor& S ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        // get eigen-pairs of C
        Tensor eVec;
        const Vector eVal = base::eigenPairs( C, eVec );
        

        // compute determinant of F and its logarithm
        const double J = mat::determinant( F );

        // scale of the stretches
        const double scale = boost::math::cbrt( 1./ J );

        // compute principal stretches
        boost::array<double,3> princStretch;
        for ( unsigned a = 0; a < 3; a++ ) {
            princStretch[a] = scale * std::sqrt( eVal[a] );
        }
        
        // get all partial derivatives
        boost::array<double,3> partialEnergyDerivatives;
        double sumPartial = 0.;
        for ( unsigned a = 0; a < 3; a++ ) {
            partialEnergyDerivatives[a] =
                this -> partialPsiIsoPartialLambdaBar_( princStretch[a] );
            sumPartial += partialEnergyDerivatives[a];
        }

        // computation via principal components and directions
        S = Tensor::Zero();
        for ( unsigned a = 0; a < 3; a ++ ) {
            const double SisoA =
                (1./eVal[a]) * (partialEnergyDerivatives[a] - sumPartial/3.);
            S += SisoA * (eVec.col(a) * (eVec.col(a)).transpose());
        }
        

    }

    
public:
    //--------------------------------------------------------------------------
    /** The iso-choric elasticity tensor in material description.
     *  \f[
     *       C^{iso} = A_{ab}   N_a \otimes N_a \otimes N_b \otimes N_b +
     *                 B_{ab} ( N_a \otimes N_b \otimes N_b \otimes N_a +
     *                          N_a \otimes N_b \otimes N_a \otimes N_b )
     *  \f]
     *  with the components
     *  \f{eqnarray*}{
     *    A_{ab} &=& \frac{1}{\lambda_b}{\partial S^{iso}_a}{\partial \lambda_b}
     *           \cr
     *    B_{ab} &=& (1-\delta_{ab})
     *                 \frac{S^{iso}_b - S^{iso}_a}{\lambda_b^2 - \lambda_a^2}
     *  \f}
     *  Note that \f$ B_{aa} = 0 \f$ and the other components evaluate to
     *  \f[
     *    A_{ab} = \frac{1}{\lambda_a^2 \lambda_b^2} \sum_{p=1}^N
     *             \mu_p \alpha_p  \left(
     *          C_{ab}^p + \frac{1}{9} \sum_{c=1}^3 \bar{\lambda}_c^{\alpha_p}
     *                            \right)
     *  \f]
     *  with
     *  \f{eqnarray*}{
     *       C_{ab,a=b}^p &=&   \frac{1}{3} \bar{\lambda}_a^{\alpha_p}  \cr
     *       C_{ab,a\neq b}^p &=& - \frac{1}{3} ( \bar{\lambda}_a^{\alpha_p} +
     *                                   \bar{\lambda}_b^{\alpha_p} ) 
     *  \f}
     *  Moreover, one has to be careful with the divided difference in the
     *  components \f$ B_{ab} \f$. Using the rule of l'H&ocirc;pital, we get
     *  \f[
     *     \lim_{\lambda_a \to \lambda_b}
     *             \frac{S^{iso}_b - S^{iso}_a}{\lambda_b^2 - \lambda_a^2} =
     *       \frac{\partial S^{iso}_b}{\partial \lambda_b^2} -
     *       \frac{\partial S^{iso}_a}{\partial \lambda_b^2}
     *  \f]
     *  Note that \f$ \partial / \partial \lambda_b^2 =
     *                (1/2\lambda_b) \partial / \partial \lambda_b \f$.
     *  The storage is based on the minor symmetries.
     *  \param[in]  F    Deformation gradient
     *  \param[out] elC  Elasticity tensor
     */
    void materialElasticityTensor( const Tensor& F, ElastTensor& elC ) const
    {
        // compute Right Cauchy-Green deformation tensor
        Tensor C;
        mat::rightCauchyGreen( F, C );

        // get eigen-pairs of C
        Tensor eVec;
        const Vector eVal = base::eigenPairs( C, eVec );

        // compute determinant of F and its logarithm
        const double J = mat::determinant( F );

        // scale of the stretches
        const double scale = boost::math::cbrt( 1./ J );

        // compute principal stretches
        boost::array<double,3> princStretch;
        for ( unsigned a = 0; a < 3; a++ ) {
            princStretch[a] = scale * std::sqrt( eVal[a] );
        }
        
        // get all partial derivatives
        boost::array<double,3> partialEnergyDerivatives;
        double sumPartial = 0.;
        for ( unsigned a = 0; a < 3; a++ ) {
            partialEnergyDerivatives[a] =
                this -> partialPsiIsoPartialLambdaBar_( princStretch[a] );
            sumPartial += partialEnergyDerivatives[a];
        }

        // compute components of S in principal directions
        Vector SisoPrinc;
        // go through principal directions
        for ( unsigned a = 0; a < 3; a++ ) {
            SisoPrinc[a] = (1./eVal[a]) * ( partialEnergyDerivatives[a] - sumPartial/3.);
        }

        // compute partial derivatives of S in principal directions
        Tensor partialSisoPrinc; // components A_{ab} above
        
        // go through principal directions
        for ( unsigned a = 0; a < 3; a++ ) {
            for ( unsigned b = 0; b < 3; b++ ) {
                const double factor = 1./eVal[a]/eVal[b];

                double outerSum = 0.;
                for ( unsigned p = 0; p < nParam; p++ ) {
                    double innerSum = 0.;
                    for ( unsigned c = 0; c < 3; c ++ ) {
                        innerSum += std::pow( princStretch[c], alpha_[p] );
                    }

                    const double aux =
                        ( a==b ?
                          std::pow(  princStretch[a], alpha_[p] ) / 3. :
                          -std::pow( princStretch[a], alpha_[p] ) / 3.
                          -std::pow( princStretch[b], alpha_[p] ) / 3. );
                    
                    outerSum += mu_[p] * alpha_[p] * (aux + innerSum/9.);
                }

                partialSisoPrinc(a,b) = factor * outerSum;
            }
        }

        // compute divided difference terms
        Tensor dividedDifference; // compute B_{ab} above
        for ( unsigned a = 0; a < 3; a++ ) {
            for ( unsigned b = 0; b < 3; b++ ) {

                if ( a != b ) {
                    
                    // squares of stretches
                    const double lamA2 = eVal[a];
                    const double lamB2 = eVal[b];

                    // tolerance check
                    if ( std::abs( lamB2 - lamA2 ) > tol_ ) {
                        // regular expression
                        dividedDifference(a,b) =
                            (SisoPrinc[b] - SisoPrinc[a]) /
                            (lamB2        - lamB2);
                    }
                    else{
                        // limit expression a l'Hopital
                        dividedDifference(a,b) = 
                            0.5 * ( partialSisoPrinc(b,b) -
                                    partialSisoPrinc(a,b) );
                    }

                }
                dividedDifference(a,b) = 0.;
            }
        }
        

        // Insert into storage matrix
        // go through all indices
        for ( unsigned A = 0; A < 3; A++ ) {
            for ( unsigned B = A; B < 3; B++ ) {

                // Voigt-Index I
                const unsigned voigtI = mat::Voigt::apply( A, B );
                
                for ( unsigned C = 0; C < 3; C++ ) {
                    for ( unsigned D = C; D < 3; D++ ) {

                        // Voigt-Index J 
                        const unsigned voigtJ = mat::Voigt::apply( C, D );

                        // go through principal directions
                        double cEntry = 0.;
                        for ( unsigned a = 0; a < 3; a++ ) {
                            for ( unsigned b = 0; b < 3; b++ ) {

                                // A_{ab} term
                                cEntry += partialSisoPrinc(a,b) *
                                    (eVec(a,A) * eVec(a,B) * eVec(b,C) * eVec(b,D));

                                // B_{ab} term
                                if ( a != b ) {
                                    cEntry += dividedDifference(a,b) * 
                                        ( (eVec(a,A) * eVec(b,B) * eVec(a,C) * eVec(b,D)) +
                                          (eVec(a,A) * eVec(b,B) * eVec(b,C) * eVec(a,D)) );
                                }
                            }
                        }
                                
                        // insert to elasticity tensor
                        elC( voigtI, voigtJ ) = cEntry;
                        
                    } // D
                } // C
                
            } // B
        } // A
        
    }
    
private:
    //--------------------------------------------------------------------------
    /** Partial derivative of the iso-choric energy with respect to principal
     *  stretch. In detail, this function returns
     *  \f[
     *   \bar{\lambda}_a \frac{\partial \psi_{iso}}{\partial \bar{\lambda}_a}
     *   = \bar{\lambda}_a \sum_{p=1}^N \mu_p \bar{\lambda}_a^{\alpha_p -1}
     *  \f]
     */
    double partialPsiIsoPartialLambdaBar_( double lambdaBar ) const
    {
        double sum = 0.;
        for ( unsigned p = 0; p < nParam; p++ ) {
            sum += mu_[p] * std::pow( lambdaBar, (alpha_[p] - 1.) );
        }
        return sum * lambdaBar;
    }

public:
    //--------------------------------------------------------------------------
    /** Ratio of first to second derivative of the volume energy.
     *  Returns
     *  \f[
     *             \frac{ \partial   \psi_{vol} }{\partial J  }
     *      \left( \frac{ \partial^2 \psi_{vol} }{\partial J^2} \right)^{-1}
     *  \f]
     *  which for the considered energy function \f$ \psi_{vol} \f$ becomes
     *  \f[
     *        J \frac{J^\beta -1}{\beta+1 - J^\beta}
     *  \f]
     *  \param[in] F  Deformation gradient
     *  \return       Ratio of derivatives
     */
    double giveBulkRatio( const Tensor& F ) const
    {
        // compute determinant of F 
        const double J = mat::determinant( F );

        const double JBeta = std::pow( J, beta_ );

        VERIFY_MSG( JBeta < beta_+1., "Too much strain input" );

        return  J * (JBeta - 1) /( beta_+1.0 - JBeta);
    }

    //--------------------------------------------------------------------------
    /** Give the inverse of the second derivative of the volume energy.
     *  Returns
     *  \f[
     *      \left( \frac{\partial^2 \psi_{vol} }{\partial J^2} \right)^{-1} =
     *        \frac{\beta J^{\beta+2}}{\kappa (\beta+1 -J^\beta)}
     *  \f]
     *  In case of incompressibility, the value zero is returned.
     *  \param[in] F   Deformation gradient  (unused in this case)
     *  \return        Inverse of bulk modulus
     */
    double giveInverseBulk( const Tensor& F) const
    {
        if ( isIncompressible_ ) return 0.;

        // compute determinant of F 
        const double J = mat::determinant( F );
        
        const double JBeta = std::pow( J, beta_ );

        return (beta_ * JBeta * J * J)/(beta_ + 1. - JBeta) / kappa_;
    }
    
private:
    //! @name Material parameters
    //@{
    const ParamArray mu_;     //!< Shear parameters
    const ParamArray alpha_;  //!< Shear exponents
    const double kappa_;      //!< Bulk modulus, unused in incompressible case
    const double beta_;       //!< Volume param, unused in incompressible case
    const bool   isIncompressible_; //!< Flag to turn on full incompressibility
    //@}

    const double tol_; //!< Evil tolerance for equality of principal stretches
};

#endif
