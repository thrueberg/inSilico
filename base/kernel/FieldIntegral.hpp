//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   kernel/FieldIntegral.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_kernel_fieldintegral_hpp
#define base_kernel_fieldintegral_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/kernel includes
#include <base/kernel/KernelFun.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace kernel{

        template<typename FIELDTUPLE>
        class FieldIntegral;
    }
}

//------------------------------------------------------------------------------
/** Kernel for the integral over a field.
 *  The integral over the field \$ u \f$ is computed as
 *  \f[
 *        M = \int_\Omega u d x
 *  \f]
 *
 *  \tparam FIELDTUPLE  Type of tuple of geometry and field for the integration
 */
template<typename FIELDTUPLE>
class base::kernel::FieldIntegral
    : public base::kernel::KernelFun<FIELDTUPLE,
                                     typename base::Vector<FIELDTUPLE::TrialElement::
                                                           DegreeOfFreedom::size>::Type>::Type
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Result type of field evaluation
    typedef typename base::Vector<TrialElement::DegreeOfFreedom::size>::Type
    VecDoF;

    //--------------------------------------------------------------------------
    //! Added weighted field result to the storage variable
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     VecDoF& fieldIntegral ) const
    {
        // Extract element pointer from tuple
        const GeomElement*   geomEp  = fieldTuple.geomElementPtr();
        const TrialElement* trialEp  = fieldTuple.trialElementPtr();

        // Evaluate the field in the given local coordinate
        const VecDoF fieldValue =
            base::post::evaluateField( geomEp, trialEp, xi );

        // element jacobian
        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );

        // summation of result
        fieldIntegral += fieldValue * detJ * weight;

        return;
    }
};

#endif
