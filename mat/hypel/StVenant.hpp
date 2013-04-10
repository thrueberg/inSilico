//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   StVenant.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef mat_hypel_stvenant_hpp
#define mat_hypel_stvenant_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/verify.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>

//------------------------------------------------------------------------------
namespace mat{
    namespace hypel{

        class StVenant;
    }
}

//------------------------------------------------------------------------------
/** St.\ Venant-Kirchhoff hyperelastic material.
 *
 *
 */
class mat::hypel::StVenant
{
public:
    
    //! @name Types of tensors used
    //@{
    typedef mat::Tensor                     Tensor;
    typedef mat::ElastTensor                ElastTensor;
    //@}
    
    //! Construct with Lame parameters
    StVenant( const double lambda, const double mu )
        : lambda_( lambda ), mu_(     mu )
    { }

    //--------------------------------------------------------------------------
    /** Strain energy density function.
     *  \f[
     *       \psi(F) = \frac{\lambda}{2} ( tr E )^2 + \mu E:E ,
     *            E  = \frac{1}{2}( F^T F - I )
     *  \f]
     *  \param[in]  F  Deformation gradient
     */
    double energy( const Tensor& F ) const
    {
        // compute Green-Lagrange strain tensor
        Tensor E;
        mat::greenLagrange( F, E );

        // compute trace of which
        const double trE = mat::trace( E );

        // return energy
        return lambda_/2. * (trE * trE) + mu_ * mat::doubleContraction( E, E );
    }

    //--------------------------------------------------------------------------
    /** The second Piola-Kirchhoff stress tensor S.
     *  The definition \f$ S = \partial \psi / \partial E \f$ yields
     *  \f[
     *         S = \lambda tr( E ) I + 2 \mu E
     *  \f]
     *  \param[in]  F  Deformation gradient
     *  \param[out] S  2nd Piola-Kirchhoff stress tensor
     */
    void secondPiolaKirchhoff( const Tensor& F, Tensor& S ) const
    {
        // compute Green-Lagrange strain tensor
        Tensor E;
        mat::greenLagrange( F, E );

        // compute trace of which
        const double trE = mat::trace( E );

        // evaluate the second Piola-Kirchhoff stress tensor
        S.noalias() = lambda_ *  trE * Tensor::Identity() + 2. * mu_ * E;
    }

    //--------------------------------------------------------------------------
    /** The elasticity tensor in material description.
     *  \f[
     *       C_{ABCD} = \lambda \delta_{AB} \delta_{CD} +
     *                  \mu (\delta_{AC} \delta_{BD} + \delta_{AD} \delta_{BC} )
     *  \f]
     *  Note: The storage is based on the minor symmetries.
     *  \param[in]  F  Deformation gradient
     *  \param[out] C  Elasticity tensor
     */
    void materialElasticityTensor( const Tensor& F, ElastTensor& C ) const
    {
        // The deformation gradient is ignored, the material behaviour in
        // reference configuration is independent of F

        // clear to zero 
        C = ElastTensor::Constant( 0. );

        // 
        C.block( 0, 0, 3, 3 ) = Tensor::Constant( lambda_ );

        //
        for ( unsigned i = 0; i < 6; i ++ ) C( i, i ) +=  2. * mu_;
    }
    
    
private:
    //! @name Lame parameters
    //@{
    const double lambda_;
    const double mu_;
    //@}
};

#endif
