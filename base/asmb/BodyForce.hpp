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
#include <base/asmb/FunEvaluationPolicy.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace asmb{

        template<typename FIELDTUPLE, typename EP>
        class BodyForce;

        namespace detail_{

            //----------------------------------------------------------------------
            template<typename FIELDTUPLEBINDER, typename QUADRATURE,
                     typename SOLVER, typename FIELDBINDER, typename BODYFORCE>
            void computeBodyForce(
                const QUADRATURE& quadrature,
                SOLVER& solver,
                const FIELDBINDER& fieldBinder,
                BODYFORCE& bodyForce )
            {
                // force integrator
                typedef ForceIntegrator<QUADRATURE,SOLVER,
                                        typename FIELDTUPLEBINDER::Tuple> ForceInt;
            
                typename ForceInt::ForceKernel forceKernel =
                    boost::bind( bodyForce, _1, _2, _3, _4 );
                ForceInt forceInt( forceKernel, quadrature, solver );

                // Apply to all elements
                typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
                typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
                for ( ; iter != end; ++iter ) {
                    forceInt( FIELDTUPLEBINDER::makeTuple( *iter ) );
                }
                return;
            }
        }

        //----------------------------------------------------------------------
        /** Convenience function for the computation of the body force term.
         *  Called with a function of type \f$ f(x) \f$, i.e. evaluated at
         *  physical coordinates.
         */
        template<typename FIELDTUPLEBINDER, typename QUADRATURE,
                 typename SOLVER, typename FIELDBINDER, typename FUN>
        void bodyForceComputation(
            const QUADRATURE& quadrature,
            SOLVER& solver,
            const FIELDBINDER& fieldBinder,
            const FUN& forceFun )
        {
            typedef base::asmb::EvaluateDirectly<
                typename FIELDTUPLEBINDER::Tuple::GeomElement,FUN> Evaluate;
            
            // body force wrapper
            BodyForce<typename FIELDTUPLEBINDER::Tuple,Evaluate> bodyForce( forceFun );


            detail_::computeBodyForce<FIELDTUPLEBINDER>( quadrature, solver,
                                                         fieldBinder, bodyForce );

            return;
        }

        //----------------------------------------------------------------------
        /** Convenience function for the computation of the body force term.
         *  Called with a function of type \f$ f( Ep*, xi ) \f$, i.e. evaluated
         *  for an element pointer and a local coordinate.
         */
        template<typename FIELDTUPLEBINDER, typename QUADRATURE,
                 typename SOLVER, typename FIELDBINDER, typename FUN>
        void bodyForceComputation2(
            const QUADRATURE& quadrature,
            SOLVER& solver,
            const FIELDBINDER& fieldBinder,
            const FUN& forceFun )
        {
            typedef base::asmb::EvaluateViaElement<
                typename FIELDTUPLEBINDER::Tuple::GeomElement,FUN> Evaluate;
            
            // body force wrapper
            BodyForce<typename FIELDTUPLEBINDER::Tuple,Evaluate> bodyForce( forceFun );


            detail_::computeBodyForce<FIELDTUPLEBINDER>( quadrature, solver,
                                                         fieldBinder, bodyForce );

            return;
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
template<typename FIELDTUPLE, typename EVALUTEPOLICY>
class base::asmb::BodyForce
    : public boost::function<void( const FIELDTUPLE&,
                                   const typename
                                   base::GeomTraits<typename FIELDTUPLE::GeomElement>::
                                   LocalVecDim&,
                                   const double, base::VectorD& ) >
{
public:
    //! @name Template parameter
    //@{
    typedef FIELDTUPLE    FieldTuple;
    typedef EVALUTEPOLICY EvaluatePolicy;
    //@}

    //! @name Extract element types
    //@{
    typedef typename FieldTuple::GeomElement GeomElement;
    typedef typename FieldTuple::TestElement TestElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Size of a DoF
    static const unsigned doFSize = TestElement::DegreeOfFreedom::size;

    //! Type of result vector
    typedef typename base::Vector<doFSize,number>::Type VecDof;

    //! How to evaluate the force function
    typedef typename EvaluatePolicy::Fun ForceFun;
    
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

        // Get Jacobian of element
        const double detJ = base::Jacobian<GeomElement>()(    geomEp, xi );

        // Evaluate force function
        const typename EvaluatePolicy::result_type f
            = EvaluatePolicy::apply( geomEp, xi, forceFun_ );

        // Evaluate the shape function
        typename TestElement::FEFun::FunArray funValues;
        (testEp -> fEFun()).evaluate( geomEp, xi, funValues );

        // deduce the size of every contribution
        const unsigned numFun = static_cast<unsigned>( funValues.size() );

        // Loop over shape functions
        for ( unsigned s = 0; s < numFun; s++ ) {

            for ( unsigned d = 0; d < doFSize; d++ )
                vector( s*doFSize+d ) += f[d] * funValues[s] * weight * detJ;
            
        }
                
        return;
        
    }
    
    //--------------------------------------------------------------------------
private:
    const ForceFun&      forceFun_;  //!< Function producing body force
};


#endif

