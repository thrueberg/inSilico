//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   BSplineShapeFun.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_bsplineshapefun_hpp
#define base_bsplineshapefun_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
// base  includes
#include <base/geometry.hpp>
#include <base/funSpace.hpp>
// base/sfun includes
#include <base/sfun/BSpline.hpp>
#include <base/sfun/TensorProduct.hpp>

//------------------------------------------------------------------------------
namespace base{

    template<unsigned DIM, unsigned DEGREE, int CONTINUITY = DEGREE-1>
    class BSplineShapeFun;
}


//------------------------------------------------------------------------------
/** Shape function based on tensor-product BSplines.
 */
template<unsigned DIM, unsigned DEGREE, int CONTINUITY>
class base::BSplineShapeFun
    : public base::sfun::TensorProduct<base::sfun::BSpline<DEGREE,CONTINUITY>,
                                       DIM, base::sfun::LEXICOGRAPHIC>
{
public:
    //! @name Template parameter: dim, degree and continuity
    //@{
    static const unsigned    dim        = DIM;
    static const unsigned    degree     = DEGREE;
    static const int         continuity = CONTINUITY;
    //@}

    //! Identifier for introspection
    static const base::FunSpace funSpace = base::BSPLINE;

    
    typedef base::sfun::TensorProduct<base::sfun::BSpline<degree,continuity>,
                                      dim,base::sfun::LEXICOGRAPHIC>  Base;

    //! Evaluate function in physical space
    template<typename GEOMELEMENT>
    void evaluate( const GEOMELEMENT* geomElemPtr,
                   const typename Base::VecDim& xi,
                   typename Base::FunArray &result ) const
    {
        Base::fun( xi, result );
    }

    //! Evaluate gradient in physical space
    template<typename GEOMELEM>
    double evaluateGradient( const GEOMELEM* geomElemPtr,
                             const typename Base::VecDim& xi,
                             std::vector<typename GEOMELEM::Node::VecDim>& result ) const
    {
        typename Base::GradArray gradXiPhi;
        Base::gradient( xi, gradXiPhi );

        // Get contra-variant basis and jacobian
        typename base::ContraVariantBasis<GEOMELEM>::MatDimLDim contraBasis;
        const double detJ =
            base::ContraVariantBasis<GEOMELEM>()( geomElemPtr, xi, contraBasis );

        // Transform every gradient
        result.resize( Base::numFun );
        for ( unsigned i = 0; i < Base::numFun; i++ )
            result[i] = contraBasis * gradXiPhi[i];

        // return det J
        return detJ;
    }
        
};

#endif
