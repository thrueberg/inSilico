//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Deformation.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef solid_deformation_hpp
#define solid_deformation_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>

//------------------------------------------------------------------------------
namespace solid{

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    void deformationGradientHistory( const GEOMELEMENT*  geomEp,
                                     const FIELDELEMENT* fieldEp,
                                     const typename FIELDELEMENT::FEFun::VecDim& xi,
                                     typename mat::Tensor& F )
    {
        // Compute the deformation gradient
        F = mat::Tensor::Identity();

        // Evaluate the displacement gradient
        const typename base::MatrixType<GEOMELEMENT::Node::dim,
                                        FIELDELEMENT::DegreeOfFreedom::size,
                                        double>::Type
            GradU = base::post::evaluateFieldGradientHistory<HIST>( geomEp,
                                                                    fieldEp, xi );
        
        // Add displacement gradient to tensor
        F.block( 0, 0, GradU.rows(), GradU.cols() ) += GradU;
    }

    //--------------------------------------------------------------------------
    /** Computation of the deformation gradient based on the current field.
     *  The deformation gradient is defined as
     *  \f[
     *        F = \frac{ \partial x }{ \partial X }
     *  \f]
     *  relating the coordinates of current configuration \f$ x \f$ with
     *  the reference one in \f$ X \f$.
     *  Using the displacement field \f$ u \f$, this tensor becomes
     *  \f[
     *       F_{iJ} = \delta_{iJ} + \frac{ \partial u_i }{ \partial X_J}
     *  \f]
     *  Note that here the displacement might have less then 3 components.
     *  It is extended by zero in the missing components.
     *  \param[in]  geomEp   Pointer to geometry element
     *  \param[in]  fieldEp  Pointer to field element
     *  \param[in]  xi       Local evaluation coordinate
     *  \param[out] F        Deformation gradient tensor
     */
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    void deformationGradient( const GEOMELEMENT*  geomEp,
                              const FIELDELEMENT* fieldEp,
                              const typename FIELDELEMENT::FEFun::VecDim& xi,
                              typename mat::Tensor& F )
    {
        deformationGradientHistory<0>( geomEp, fieldEp, xi, F );
    }

    //--------------------------------------------------------------------------
    /**
     */
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    double jacobian( const GEOMELEMENT* geomEp,
                     const FIELDELEMENT* fieldEp,
                     const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        mat::Tensor F;
        deformationGradientHistory<0>( geomEp, fieldEp, xi, F );

        return mat::determinant( F );
    }

    template<typename GEOMELEMENT, typename FIELDELEMENT>
    double jacobian( const GEOMELEMENT* geomEp,
                     const FIELDELEMENT* fieldEp )
    {
        return jacobian( geomEp, fieldEp,
                         base::ShapeCentroid<GEOMELEMENT::shape>::apply() );
    }
                   

}
#endif
