#ifndef fieldintegralpart_h
#define fieldintegralpart_h

#include <utility>
#include <base/kernel/KernelFun.hpp>
#include <base/linearAlgebra.hpp>
#include <base/geometry.hpp>
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
// Integrate over the area and the field value and assign the values to
// different regions according to the physical coordinate of the integration
// point
template<typename FIELDTUPLE>
class FieldIntegralPart
    : public base::kernel::KernelFun<FIELDTUPLE,
                                     std::pair<base::Vector<3>::Type,
                                               base::Vector<3>::Type> >::Type
{
public:
    typedef FIELDTUPLE FieldTuple;
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TrialElement TrialElement;

    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim  GlobalVecDim;
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim   LocalVecDim;

    typedef base::Vector<1>::Type Vec1;
    typedef base::Vector<3>::Type Vec3;
    typedef std::pair<Vec3,Vec3> ResultPair;

    FieldIntegralPart( const double x1, const double x2 )
        : x1_( x1 ), x2_( x2 ) {  }

    //--------------------------------------------------------------------------
    //! Added weighted field result to the storage variable
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     ResultPair& fieldIntegral ) const
    {
        // Extract element pointer from tuple
        const GeomElement*   geomEp  = fieldTuple.geomElementPtr();
        const TrialElement* trialEp  = fieldTuple.trialElementPtr();

        const GlobalVecDim x = base::Geometry<GeomElement>()( geomEp, xi );

        // decide the region in which the integration point lies
        unsigned domainNum = 0;
        if      ( x[0] <= x1_ ) domainNum = 1;
        else if ( x[0] <  x2_ ) domainNum = 2;
        else                    domainNum = 3;
        
        // Evaluate the field in the given local coordinate
        const Vec1 fieldValue =
            base::post::evaluateField( geomEp, trialEp, xi );

        // element jacobian
        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );

        // summation of result
        //fieldIntegral += fieldValue * detJ * weight;
        fieldIntegral.first[  domainNum-1 ] += detJ * weight;
        fieldIntegral.second[ domainNum-1 ] += fieldValue[0] * detJ * weight;

        return;
    }

private:
    const double x1_;
    const double x2_;
};


#endif
