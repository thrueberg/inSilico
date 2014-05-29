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

    namespace detail_{
        template<typename SELEMENT,bool ISSURF>
        struct EnableSurfaceNormal
        {
            static double apply( const SELEMENT* sep,
                                 const typename base::GeomTraits<SELEMENT>::LocalVecDim& xi,
                                 typename base::GeomTraits<SELEMENT>::GlobalVecDim& normal )
            {
                return base::SurfaceNormal<SELEMENT>()( sep, xi, normal );
            }
        };

        template<typename SELEMENT>
        struct EnableSurfaceNormal<SELEMENT,false>
        {
            static double apply( const SELEMENT* sep,
                                 const typename base::GeomTraits<SELEMENT>::LocalVecDim& xi,
                                 typename base::GeomTraits<SELEMENT>::GlobalVecDim& normal )
            {
                normal = base::constantVector<base::GeomTraits<SELEMENT>::globalDim>( 0. );
                return 0.;
            }
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

    //--------------------------------------------------------------------------
    /** Evaluate function in physical space.
     *  The function values in the physical coordinate system are defined by
     *  \f[
     *        \phi(x) := \phi( \xi(x) ) = \phi \circ x^{-1}
     *  \f]
     *  where \f$ x = x(\xi) \f$ is the coordinate map.
     *  \tparam GEOMELEMENT Type of geometry description element
     *  \param  geomElemPtr Pointer to geometry element (here unused)
     *  \param  xi          Local evaluation coordinate
     *  \param  result      Function values of all shape functions
     */
    template<typename GEOMELEMENT>
    void evaluate( const GEOMELEMENT* geomElemPtr,
                   const typename Base::VecDim& xi,
                   typename Base::FunArray &result ) const
    {
        Base::fun( xi, result );
    }

    //--------------------------------------------------------------------------
    /** Evaluate gradient in physical space.
     *  Compute the (tangential) gradient of the shape functions
     *  \f[
     *       \nabla \phi = g^\alpha \phi_{,\alpha}
     *  \f]
     *  where \f$ g^\alpha \f$ are the contra-variante basis vectors and
     *  \f$ \phi_{,\alpha} \f$ the partial derivatives in parametric coordinates.
     *  \tparam GEOMELEMENT Type of geometry description element
     *  \param  geomElemPtr Pointer to geometry element 
     *  \param  xi          Local evaluation coordinate
     *  \param  result      Gradient values of all shape functions
     */
    template<typename GEOMELEM>
    double evaluateGradient( const GEOMELEM* geomElemPtr,
                             const typename Base::VecDim& xi,
                             std::vector<typename GEOMELEM::Node::VecDim>& result ) const
    {
        // get local gradient
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

    //--------------------------------------------------------------------------
    /** Evaluate Hessian in physical space.
     *  Compute the (tangential) Hessian of the shape functions
     *  \f[
     *       \nabla \otimes \nabla \phi =
     *       (\phi_{,\alpha\beta} - \Gamma^\delta_{\alpha \beta} \phi_{,\delta})
     *       g^\alpha \otimes g^\beta +
     *       \phi_{,alpha} L^\alpha_\beta g^\beta \otimes n
     *  \f]
     *  where \f$ g^\alpha \f$ are the contra-variante basis vectors, 
     *  \f$ \phi_{,\alpha} \f$ the partial derivatives in parametric coordinates,
     *  and we have the Christoffel symbol
     *  \f[
     *       \Gamma^\delta_{\alpha \beta} = g_{\alpha,\beta} \cdot g^\delta
     *  \f]
     *  based on the partial derivatives of the co-variante basis vectors
     *  \f$ g_\alpha \f$. Moreover, the second fundamental form
     *  \f[
     *       L^\alpha_\beta
     *         = G^{\alpha \gamma} L_{\gamma \beta}
     *         = G^{\alpha \gamma} (g_{\gamma, \beta} \cdot n)
     *  \f]
     *  is used for the case of surface derivatives.
     *
     *  \tparam GEOMELEMENT Type of geometry description element
     *  \param  geomElemPtr Pointer to geometry element 
     *  \param  xi          Local evaluation coordinate
     *  \param  result      Gradient values of all shape functions
     */
    template<typename GEOMELEM>
    double evaluateHessian( const GEOMELEM* geomElemPtr,
                            const typename Base::VecDim& xi,
                            std::vector<typename base::Matrix<
                            GEOMELEM::Node::dim,GEOMELEM::Node::dim>::Type
                            >& result ) const
    {
        // parameter space gradient
        typename Base::GradArray gradXi;
        Base::gradient(xi, gradXi );

        // parameter space Hessian
        typename Base::HessianArray hessXi;
        Base::hessian( xi, hessXi );

        // geometry items: contra-variant basis and its derivatives
        typename base::ContraVariantBasis<GEOMELEM>::MatDimLDim contraBasis;
        const double detJ = 
            base::ContraVariantBasis<GEOMELEM>()( geomElemPtr, xi, contraBasis );

        typename base::DerivativesOfTangents<GEOMELEM>::result_type tangDeriv =
            base::DerivativesOfTangents<GEOMELEM>()( geomElemPtr, xi );

        // compute normal vector in case of a surface derivative
        static const bool isSurf =
            (base::GeomTraits<GEOMELEM>::globalDim ==
             base::GeomTraits<GEOMELEM>::localDim+1);
        typename base::Vector<GEOMELEM::Node::dim,double>::Type normal;
        detail_::EnableSurfaceNormal<GEOMELEM,isSurf>::
            apply( geomElemPtr, xi, normal );

        // 2nd ff
        typename base::Matrix<Base::dim,Base::dim>::Type FF, Ginv, FF2;
        for ( unsigned a = 0; a < Base::dim; a++ ) {
            for ( unsigned b = 0; b < Base::dim; b++ ) {
                FF(a,b)   = tangDeriv(a,b).dot( normal );
                Ginv(a,b) = contraBasis.col(a).dot( contraBasis.col(b) );
            }
        }
        FF2 = Ginv * FF;

        result.resize( Base::numFun );
        for ( unsigned i = 0; i < Base::numFun; i++ ) {

            // initalise with zeros
            typename base::Matrix<GEOMELEM::Node::dim,
                                  GEOMELEM::Node::dim>::Type hessX
                = base::constantMatrix<GEOMELEM::Node::dim,
                                       GEOMELEM::Node::dim>(0.);
            
            for ( unsigned a = 0; a < Base::dim; a++ ) {
                for ( unsigned b = 0; b < Base::dim; b++ ) {
                    
                    // partial second derivative
                    const double aux1 = hessXi[i](a,b);

                    // Christoffel term
                    double aux2 = 0.;
                    for ( unsigned d = 0; d < Base::dim; d++ )
                        aux2 += gradXi[i][d] *
                            ( tangDeriv(a,b).dot( contraBasis.col(d) ) );

                    // normal term
                    const double aux3 = gradXi[i][a] * FF2(a,b);

                    //
                    hessX += (aux1 - aux2) *
                        (contraBasis.col(a) * contraBasis.col(b).transpose() ) +
                        aux3 * (contraBasis.col(b) * normal.transpose());
                }
            }

            result[i] = hessX;
        }


        return detJ;
    }

};

#endif
