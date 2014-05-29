//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   cut/ParametricMeasure.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_cut_parametricmeasure_hpp
#define base_cut_parametricmeasure_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/linearAlgebra.hpp>
// base/kernel includes
#include <base/kernel/KernelFun.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename FIELDTUPLE>
        class ParametricMeasure;
    }
}

//------------------------------------------------------------------------------
/** Kernel for area/volume computations in the parameter space.
 *  Using the simple formula
 *  \f[
 *        V = \int_{\hat{\tau}} d\xi
 *  \f]
 *  for conversion from physical to parameter space, the parametric area of an
 *  element becomes the sum of the quadrature weights.
 *  Note that this seemingly trivial task is of interest in case of an immersed
 *  method with cut elements. In that case the quadrature is carried out only
 *  over the physically active portion of the elements.
 *  \tparam FIELDTUPLE  Type of tuple of geometry and field for the integration
 */
template<typename FIELDTUPLE>
class base::cut::ParametricMeasure
    : public base::kernel::KernelFun<FIELDTUPLE,double>::Type
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! Local coordinate
    typedef typename
    base::GeomTraits<typename FieldTuple::GeomElement>::LocalVecDim
    LocalVecDim;

    //--------------------------------------------------------------------------
    //! Add quadrature weight to get the parametric mesarure
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     double& sumOfMeasures ) const
    {
        sumOfMeasures += weight;
    }
};

#endif
