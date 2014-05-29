//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   kernel/Measure.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_kernel_measure_hpp
#define base_kernel_measure_hpp

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
        class Measure;

        template<typename FIELDTUPLE>
        class Moment;
        
        template<typename FIELDTUPLE>
        class SecondMoment;
    }
}

//------------------------------------------------------------------------------
/** Kernel for area/volume computations.
 *  Using the simple formula
 *  \f[
 *        V = \int_\tau dx = \int_{\hat{\tau}} det J d\xi
 *  \f]
 *  for conversion from physical to parameter space, the area of an element
 *  becomes the weighted sum of its metrics, i.e. the \f$ det J \f$ values.
 *  \tparam FIELDTUPLE  Type of tuple of geometry and field for the integration
 */
template<typename FIELDTUPLE>
class base::kernel::Measure
    : public base::kernel::KernelFun<FIELDTUPLE,double>::Type
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    //@}
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //--------------------------------------------------------------------------
    //! Added weighted metric to the result storatge
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     double& sumOfMeasures ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();

        // element jacobian
        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );

        sumOfMeasures += detJ * weight;

    }
};

//------------------------------------------------------------------------------
/** Kernel for volume (area) moment computation.
 *  The volume moment \f$ S \f$ is defined by the integral
 *  \f[
 *       S = \int_\Omega x d x
 *  \f]
 *  and computed by this function using the transformation to the reference
 *  element. This quantity is typically used to compute the centroid of
 *  a volume by means of
 *  \f[
 *       C = \frac{1}{V} S
 *  \f]
 *  where \f$ V \f$ is the volume of the domain \f$ \Omega \f$.
 *  \tparam FIELDTUPLE  Type of tuple of geometry and field for the integration
 */
template<typename FIELDTUPLE>
class base::kernel::Moment
    : public base::kernel::KernelFun<
    FIELDTUPLE,
    typename base::Vector<FIELDTUPLE::GeomElement::Node::dim>::Type>::Type
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    //@}
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim   LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim  GlobalVecDim;
    
    //--------------------------------------------------------------------------
    //! Added weighted metric to the result storatge
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     GlobalVecDim&      sumOfMoments ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();

        // evalute geometry
        const GlobalVecDim x = base::Geometry<GeomElement>()( geomEp, xi );

        // element jacobian
        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );

        sumOfMoments += x * detJ * weight;
    }
};

//------------------------------------------------------------------------------
/** Kernel for 2nd moment of volume (area) computation.
 *  The 2nd moment of volume is defined as
 *  \f[
 *      J = \int_\Omega x \otimes x d x
 *  \f]
 *  and is a square matrix of size of the spatial dimension. By means of this
 *  quantity the moment of inertia (assuming a homogeneous material) becomes
 *  \f[
 *      I = J - V (C \otimes C)
 *  \f]
 *  based on the volume \f$ V \f$ and the centroid \f$ C \f$ of the domain
 *  \f$ \Omega \f$.
 *  \tparam FIELDTUPLE  Type of tuple of geometry and field for the integration
 */
template<typename FIELDTUPLE>
class base::kernel::SecondMoment
    : public base::kernel::KernelFun<
    FIELDTUPLE,
    typename base::Matrix<FIELDTUPLE::GeomElement::Node::dim,
                          FIELDTUPLE::GeomElement::Node::dim>::Type>::Type
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    //@}
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim   LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim  GlobalVecDim;

    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;
    typedef typename base::Matrix<globalDim,globalDim,double>::Type  Result;
    
    //--------------------------------------------------------------------------
    //! Added weighted metric to the result storatge
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     Result&            sumOfMoments ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();

        // evalute geometry
        const GlobalVecDim x = base::Geometry<GeomElement>()( geomEp, xi );

        // element jacobian
        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );

        for ( unsigned d1 = 0; d1 < globalDim; d1++ ) {
            for ( unsigned d2 = 0; d2 < globalDim; d2++ ) {
        
                sumOfMoments(d1,d2) += x[d1] * x[d2] * detJ * weight;
            }
        }
    }
};

#endif
