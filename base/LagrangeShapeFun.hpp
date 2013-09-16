//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   LagrangeShapeFun.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_lagrangeshapefun_hpp
#define base_lagrangeshapefun_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
// base  includes
#include <base/shape.hpp>
#include <base/geometry.hpp>
#include <base/funSpace.hpp>
// basse/sfun includes
#include <base/sfun/ZeroDimensionalDummy.hpp>
#include <base/sfun/Lagrange1D.hpp>
#include <base/sfun/LagrangeTriangle.hpp>
#include <base/sfun/LagrangeTetrahedron.hpp>
#include <base/sfun/TensorProduct.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace detail_{

        template<unsigned DEGREE, base::Shape SHAPE>
        struct LagrangeShapeFunImpl;
        
        //----------------------------------------------------------------------
        //! Lagrangian shape function on a line
        template<unsigned DEGREE>
        struct LagrangeShapeFunImpl<DEGREE,base::POINT>
        {
            typedef base::sfun::ZeroDimensionalDummy<DEGREE> Type;
        };

        //----------------------------------------------------------------------
        //! Lagrangian shape function on a line
        template<unsigned DEGREE>
        struct LagrangeShapeFunImpl<DEGREE,base::LINE>
        {
            typedef base::sfun::TensorProduct<base::sfun::Lagrange1D<DEGREE>,
                                              1, base::sfun::HIERARCHIC>
            Type;
        };

        //----------------------------------------------------------------------
        //! Lagrangian shape function on a triangle
        template<unsigned DEGREE>
        struct LagrangeShapeFunImpl<DEGREE,base::TRI>
        {
            typedef base::sfun::LagrangeTriangle<DEGREE> Type;
        };

        //----------------------------------------------------------------------
        //! Tensor-product Lagrangian shape function on a quadrilateral
        template<unsigned DEGREE>
        struct LagrangeShapeFunImpl<DEGREE,base::QUAD>
        {
            typedef base::sfun::TensorProduct< base::sfun::Lagrange1D<DEGREE>,
                                               2, base::sfun::HIERARCHIC>
            Type;
        };

        //----------------------------------------------------------------------
        //! Lagrangian shape function on a tetrahedron
        template<unsigned DEGREE>
        struct LagrangeShapeFunImpl<DEGREE,base::TET>
        {
            typedef base::sfun::LagrangeTetrahedron<DEGREE> Type;
        };
        
        //----------------------------------------------------------------------
        //! Tensor-product Lagrangian shape function on a hexahedron
        template<unsigned DEGREE>
        struct LagrangeShapeFunImpl<DEGREE,base::HEX>
        {
            typedef base::sfun::TensorProduct< base::sfun::Lagrange1D<DEGREE>,
                                               3, base::sfun::HIERARCHIC>
            Type;
        };
    }

    template<unsigned DEGREE, base::Shape SHAPE>
    class LagrangeShapeFun;
}

//------------------------------------------------------------------------------
/** Shape function of the Lagrangian type.
 *  Inherits from specific implementations as given in the sfun namesapce.
 *  Moreover, it provides the interfaces evaluate and evaluateGradient which,
 *  given access to a geometry element, allow for the proper mapping from
 *  parameter to physical space. Note that in the function evaluation case
 *  this is trivial, but the interface remains for conformity with other
 *  types of shape functions.
 *  \tparam DEGREE  Polynomial degree of the shape function
 *  \tparam SHAPE   Geometric chape of the reference element
 */
template<unsigned DEGREE, base::Shape SHAPE>
class base::LagrangeShapeFun
    : public base::detail_::LagrangeShapeFunImpl<DEGREE,SHAPE>::Type
{
public:
    //! @name Template parameter: degree and shape
    //@{
    static const unsigned    degree = DEGREE;
    static const base::Shape shape  = SHAPE;
    //@}

    //! Identifier for introspection
    static const base::FunSpace funSpace = base::LAGRANGE;

    //! Base class contains the parametric shape function implementation
    typedef typename detail_::LagrangeShapeFunImpl<DEGREE,SHAPE>::Type Base;

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
