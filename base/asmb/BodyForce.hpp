//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   BodyForce.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_bodyforce_hpp
#define base_asmb_bodyforce_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/asmb includes
#include <base/asmb/ForceIntegrator.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace asmb{

        template<typename FIELDTUPLE>
        class BodyForce;

        //----------------------------------------------------------------------
        /** Convenience function for the computation of the body force term.
         *
         */
        template<typename QUADRATURE, typename SOLVER, typename FIELDBINDER>
        void bodyForceComputation(
            const QUADRATURE& quadrature,
            SOLVER& solver,
            const FIELDBINDER& fieldBinder,
            const typename BodyForce<typename FIELDBINDER::ElementPtrTuple>::ForceFun& ff )
        {
            // body force wrapper
            BodyForce<typename FIELDBINDER::ElementPtrTuple> bodyForce( ff );
            
            // force integrator
            typedef ForceIntegrator<QUADRATURE,SOLVER,
                                    typename FIELDBINDER::ElementPtrTuple> ForceInt;
            
            typename ForceInt::ForceKernel forceKernel =
                boost::bind( bodyForce, _1, _2, _3, _4 );
            ForceInt forceInt( forceKernel, quadrature, solver );
            std::for_each( fieldBinder.elementsBegin(), fieldBinder.elementsEnd(),
                           forceInt );
        }


    }
}
//------------------------------------------------------------------------------
/** Wrapper around a body force function.
 *  Given any forcing function \f$ f(x) \f$, the corresponding linear form reads
 *  \f[
 *      L(\phi) = \int_\Omega f(x) \phi(x) d x
 *  \f]
 *  Here, the integral kernel function \f$ f(x) \phi(x) \f$ is represented. The
 *  computation of the integral and the assembly is done outside in
 *  base::ForceIntegrator.
 *  \tparam FIELDTUPLE   Tuple of field element pointers
 */
template<typename FIELDTUPLE>
class base::asmb::BodyForce
    : public boost::function<void( const FIELDTUPLE&,
                                   const typename
                                   base::GeomTraits<typename FIELDTUPLE::GeomElement>::
                                   LocalVecDim&,
                                   const double, base::VectorD& ) >
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types
    //@{
    typedef typename FieldTuple::GeomElement GeomElement;
    typedef typename FieldTuple::TestElement TestElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Size of a DoF
    static const unsigned doFSize = TestElement::DegreeOfFreedom::size;

    //! Type of result vector
    typedef typename base::Vector<doFSize,number>::Type VecDof;

    //! Type of force function
    typedef boost::function<VecDof( const GlobalVecDim& )> ForceFun;

    //--------------------------------------------------------------------------
    //! Constructor setting all references
    BodyForce( const ForceFun& forceFun  )
        : forceFun_( forceFun )
    { }

    //--------------------------------------------------------------------------
    /** Evaluation of the kernel function due to a body force term.
     *  \param[in]   fieldTuple Tuple of field elements
     *  \param[in]   xi         Local evaluation coordinate
     *  \param[in]   weight     Corresponding quadrature weight
     *  \param[out]  vector     Result container (pre-sized and zero-initialised)
     */
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     base::VectorD&     vector ) const
    {
        // extract test and trial elements from tuple
        const GeomElement* geomEp = fieldTuple.geomElementPtr();
        const TestElement* testEp = fieldTuple.testElementPtr();
        
        // Evaluate element geometry
        const GlobalVecDim x = base::Geometry<GeomElement>()( geomEp, xi );

        // Get Jacobian of element
        const double detJ = base::Jacobian<GeomElement>()(    geomEp, xi );

        // Evaluate force function
        const typename ForceFun::result_type f = forceFun_( x );

        // Evaluate the shape function
        typename TestElement::FEFun::FunArray funValues;
        (testEp -> fEFun()).evaluate( geomEp, xi, funValues );

        // deduce the size of every contribution
        const unsigned numFun = static_cast<unsigned>( funValues.size() );

        // Loop over shape functions
        for ( unsigned s = 0; s < numFun; s++ ) {

            vector.segment( s * doFSize,
                            doFSize )  += f * funValues[s] * weight * detJ;
        }
                
        return;
        
    }
    
    //--------------------------------------------------------------------------
private:
    const ForceFun&      forceFun_;  //!< Function producing body force
};


#endif

