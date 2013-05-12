//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Stress.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef solid_stress_hpp
#define solid_stress_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>
// solid includes
#include <solid/Deformation.hpp>

//------------------------------------------------------------------------------
namespace solid{

    /** Computation of the Cauchy stress tensor.
     *  Given the state of deformation and the material law, the Caucy
     *  stress tensor will be computed by means of the Piola-transform
     *  \f[
     *       \sigma = \frac{1}{\det (F)} F S F^T
     *  \f]
     *  using the deformation gradient \f$ F \f$ and the 2nd Piola-Kirchhoff
     *  stress tensor \f$ S \f$.
     *  \tparam GEOMELEMENT  Geometry element
     *  \tparam FIELDELEMENT Field element
     *  \tparam MATERIAL     Material type
     */
    template<typename GEOMELEMENT, typename TRIALELEMENT, typename MATERIAL>
    mat::Tensor cauchy( const GEOMELEMENT*  geomEp,
                        const TRIALELEMENT* trialEp,
                        const MATERIAL& material,
                        const typename TRIALELEMENT::FEFun::VecDim& xi )
    {
        // get deformation gradient
        mat::Tensor F;
        deformationGradient( geomEp, trialEp, xi, F );

        // Material evaluations
        typename mat::Tensor S;      // 2nd PK
        material.secondPiolaKirchhoff( F, S );

        // Push-forward and divide by det(F)
        const double J = mat::determinant( F );
        mat::Tensor sigma;
        sigma.noalias() = (1./J) * F * S * F.transpose();

        return sigma;
    }

    //! Overload call of solid::cauchy for element centroid
    template<typename GEOMELEMENT, typename TRIALELEMENT, typename MATERIAL>
    mat::Tensor cauchy( const GEOMELEMENT*  geomEp,
                        const TRIALELEMENT* trialEp,
                        const MATERIAL& material )
    {
        return cauchy<GEOMELEMENT,TRIALELEMENT,MATERIAL>(
            geomEp, trialEp, material,
            base::ShapeCentroid<GEOMELEMENT::shape>::apply() );
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename PRESSELEMENT>
    double pressureHistory( const GEOMELEMENT*  geomEp,
                            const PRESSELEMENT* pressEp,
                            const typename PRESSELEMENT::FEFun::VecDim& xi )
    {
        // Evaluate the displacement gradient
        const typename base::Vector<PRESSELEMENT::DegreeOfFreedom::size,
                                        double>::Type
            p = base::post::evaluateFieldHistory<HIST>( geomEp, pressEp, xi );
        
        // return the entry
        return p[0];
    }

    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename PRESSELEMENT>
    double pressure( const GEOMELEMENT*  geomEp,
                     const PRESSELEMENT* pressEp,
                     const typename PRESSELEMENT::FEFun::VecDim& xi )
    {
        return pressureHistory<0>( geomEp, pressEp, xi );
    }
    
}
#endif
