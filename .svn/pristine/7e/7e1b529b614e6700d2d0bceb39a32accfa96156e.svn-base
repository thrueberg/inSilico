//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Quadrature.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_quadrature_hpp
#define base_quadrature_hpp

//------------------------------------------------------------------------------
// base  includes
#include <base/quad/GaussLegendre.hpp>
#include <base/quad/GaussTriangle.hpp>
#include <base/quad/GaussTetrahedron.hpp>
#include <base/quad/TensorProduct.hpp>
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{

    template<unsigned DEGREE, base::Shape SHAPE>
    class Quadrature;

    //! \cond SKIPDOX
    //--------------------------------------------------------------------------
    //! Quadrature on a line, by default Gauss-Legendre
    template<unsigned DEGREE>
    class Quadrature<DEGREE,base::LINE>
        : public base::quad::GaussLegendre<DEGREE>
    { };

    //--------------------------------------------------------------------------
    //! Quadrature on a triangle
    template<unsigned DEGREE>
    class Quadrature<DEGREE,base::TRI>
        : public base::quad::GaussTriangle<DEGREE>
    { };

    //-------------------------------------------------------------------------
    //! Quadrature on a quadrilateral, be default tensor-product GaussLeg
    template<unsigned DEGREE>
    class Quadrature<DEGREE,base::QUAD>
        : public base::quad::TensorProduct< base::quad::GaussLegendre<DEGREE>,
                                            2 >
    { };

    //--------------------------------------------------------------------------
    //! Quadrature on a tetrahedron
    template<unsigned DEGREE>
    class Quadrature<DEGREE,base::TET>
        : public base::quad::GaussTetrahedron<DEGREE>
    { };

    //-------------------------------------------------------------------------
    //! Quadrature on a quadrilateral, be default tensor-product GaussLeg
    template<unsigned DEGREE>
    class Quadrature<DEGREE,base::HEX>
        : public base::quad::TensorProduct< base::quad::GaussLegendre<DEGREE>,
                                            3 >
    { };

    //! \endcond

    //--------------------------------------------------------------------------
    //! Quadrature over the surface of a shape
    template<unsigned DEGREE, base::Shape SHAPE>
    class SurfaceQuadrature
        : public base::Quadrature<DEGREE,
                                  base::FaceShape<SHAPE>::value >
    { };
}

#endif
