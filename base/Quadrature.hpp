//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   base/Quadrature.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_quadrature_hpp
#define base_quadrature_hpp

//------------------------------------------------------------------------------
// base  includes
#include <base/quad/PointEvaluation.hpp>
#include <base/quad/GaussLegendre.hpp>
#include <base/quad/GaussTriangle.hpp>
#include <base/quad/GaussTetrahedron.hpp>
#include <base/quad/TensorProduct.hpp>
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace detail_{

        //! Helper to define from which quadrature to inherit the implementation
        template<unsigned DEGREE, base::Shape SHAPE>
        struct QuadratureBase;

        //----------------------------------------------------------------------
        //! Quadrature on a point, point evaluation
        template<unsigned DEGREE>
        struct QuadratureBase<DEGREE,base::POINT>
        {
            typedef base::quad::PointEvaluation Type;
        };

        //----------------------------------------------------------------------
        //! Quadrature on a line, by default Gauss-Legendre
        template<unsigned DEGREE>
        struct QuadratureBase<DEGREE,base::LINE>
        {
            typedef base::quad::GaussLegendre<DEGREE> Type;
        };

        //----------------------------------------------------------------------
        //! Quadrature on a triangle
        template<unsigned DEGREE>
        struct QuadratureBase<DEGREE,base::TRI>
        {
            typedef base::quad::GaussTriangle<DEGREE> Type;
        };

        //---------------------------------------------------------------------
        //! Quadrature on a quadrilateral, be default tensor-product GaussLeg
        template<unsigned DEGREE>
        struct QuadratureBase<DEGREE,base::QUAD>
        {
            typedef
            base::quad::TensorProduct<base::quad::GaussLegendre<DEGREE>,2> Type;
        };

        //----------------------------------------------------------------------
        //! Quadrature on a tetrahedron
        template<unsigned DEGREE>
        struct QuadratureBase<DEGREE,base::TET>
        {
            typedef base::quad::GaussTetrahedron<DEGREE> Type;
        };
        
        //---------------------------------------------------------------------
        //! Quadrature on a quadrilateral, be default tensor-product GaussLeg
        template<unsigned DEGREE>
        struct QuadratureBase<DEGREE,base::HEX>
        {
            typedef
            base::quad::TensorProduct<base::quad::GaussLegendre<DEGREE>,3> Type;
        };
    }

    //--------------------------------------------------------------------------
    /** Main object for numerical integration.
     *  The integral of the form
     *  \f[
     *       I[f] = \int_{\hat \tau} f(\xi) d\xi
     *  \f]
     *  with the integrand kernel \f$ f \f$ and the shape \f$ \hat{\tau} \f$ in
     *  reference space is approximated by the weighted sum
     *  \f[
     *       Q[f] = \sum_{i=1}^N f(\xi_i) w_i
     *  \f]
     *  with the quadrature points \f$ \xi_i \in \hat{\tau} \f$ and the
     *  corresponding weights \f$ w_i \f$.
     *
     *  This object inherits its specific implementation based on the values of
     *  DEGREE and SHAPE:
     *    - DEGREE denotes the degree of any coordiante monomial which is
     *      integrated exactly by the rule
     *    - SHAPE is the base::Shape which represents the integration region
     *      in parametric coordinates
     *
     *  Moreover, the interface of begin and end iterators to the storage of
     *  weighted quadrature points is inherited. In addition, this object
     *  provides an apply function which calls a kernel function of type
     *  \code{.cpp}
     *  void( Arg1&, const VecDim&, const double, Arg4& )
     *  \endcode
     *  for every quadrature point and weight.
     *  \tparam DEGREE Polynomial degree to be integrated exactly
     *  \tparam SHAPE  Convex integration domain in parameter space
     */
    template<unsigned DEGREE, base::Shape SHAPE>
    class Quadrature
        : public base::detail_::QuadratureBase<DEGREE,SHAPE>::Type
    {
        typedef typename base::detail_::QuadratureBase<DEGREE,SHAPE>::Type Base;
        
    public:
        //----------------------------------------------------------------------
        /** Generic apply function.
         *  Evaluate provided kernel function object at every pair of weight and
         *  point of the quadrature rule.
         *  \tparam KERNEL Type of function object representing the
         *                 integral kernel
         *  \param[in] kernel   Reference to kernel object
         *  \param[in] arg1     First argument of the kernel function
         *                      (typically, a tuple of element pointers)
         *  \param[in,out] arg4 Fourth argument of the kernel function
         *                      (typically, a result container)
         */
        template<typename KERNEL>
        void apply( KERNEL& kernel,
                    typename KERNEL::arg1_type& arg1,
                    typename KERNEL::arg4_type& arg4 ) const
        {
            typename Base::Iter qIter = this -> begin();
            typename Base::Iter qEnd  = this -> end();
            for ( ; qIter != qEnd; ++qIter )
                kernel( arg1, qIter -> second, qIter -> first, arg4 );
        }

    };


    //--------------------------------------------------------------------------
    //! Quadrature over the surface of a shape
    template<unsigned DEGREE, base::Shape SHAPE>
    class SurfaceQuadrature
        : public base::Quadrature<DEGREE,base::FaceShape<SHAPE,1>::value>
    { };
}

#endif
