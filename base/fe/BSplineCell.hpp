//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   BSplineCell.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_fe_bsplinecell_hpp
#define base_fe_bsplinecell_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/shape.hpp>
#include <base/meta.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace fe{


        //----------------------------------------------------------------------
        //! Lagrangian DIM-HyperCube of polynomial DEGREE.
        template<unsigned DIM, unsigned DEGREE, int CONTINUITY>
        struct BSplineCell
        {
            //! 1D Shape function            
            typedef base::sfun::BSpline<DEGREE,CONTINUITY> BSpline;

            //! Tensor-product gives element shape function
            typedef base::sfun::TensorProduct<BSpline,DIM> ShapeFun;

            //! Total number of DoFs
            static const unsigned numTotalDoFs = ShapeFun::numFun;
        };

    }
}

#endif
