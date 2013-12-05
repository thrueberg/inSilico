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
#include <boost/typeof/typeof.hpp>
#include <boost/function_types/function_arity.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/asmb includes
#include <base/asmb/ForceIntegrator.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace asmb{

        template<typename FIELDTUPLE, typename EP>
        class BodyForce;

        namespace detail_{

            //----------------------------------------------------------------------
            template<typename FIELDTUPLEBINDER>
            struct ElementTypeBinder
            {
                typedef typename FIELDTUPLEBINDER::Tuple      Tuple;
                typedef typename Tuple::GeomElement           GeomElement;
                typedef typename GeomElement::GeomFun::VecDim LocalVecDim;
                typedef typename GeomElement::Node::VecDim    GlobalVecDim;

                typedef typename Tuple::TrialElement           TrialElement;
                typedef typename TrialElement::DegreeOfFreedom DegreeOfFreedom;
                typedef typename
                base::Vector<DegreeOfFreedom::size>::Type VecDof;
            };


            //----------------------------------------------------------------------
            template<typename FIELDTUPLEBINDER>
            struct EvaluateViaElement
            {
                typedef ElementTypeBinder<FIELDTUPLEBINDER> ETB;
                
                typedef typename ETB::VecDof                   result_type;
                typedef typename ETB::GeomElement              GeomElement;
                typedef typename ETB::LocalVecDim              LocalVecDim;

                typedef boost::function<result_type( const GeomElement*,
                                                     const LocalVecDim& ) > Fun;

                static result_type apply( const GeomElement* gep,
                                          const LocalVecDim& xi, 
                                          Fun& fFun )
                {
                    return fFun( gep, xi );
                }
            };

            //----------------------------------------------------------------------
            template<typename FIELDTUPLEBINDER>
            struct EvaluateDirectly
            {
                typedef ElementTypeBinder<FIELDTUPLEBINDER> ETB;
                
                typedef typename ETB::VecDof                   result_type;
                typedef typename ETB::GeomElement              GeomElement;
                typedef typename ETB::GlobalVecDim             GlobalVecDim;

                typedef boost::function<result_type( const GlobalVecDim& )> Fun;
                
                static result_type apply( const GeomElement* gep,
                                          const typename GeomElement::GeomFun::VecDim& xi,
                                          Fun& fFun )
                {
                    return fFun( base::Geometry<GeomElement>()( gep, xi ) );
                }
            };

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
                 typename SOLVER, typename FIELDBINDER>
        void bodyForceComputation(
            const QUADRATURE& quadrature,
            SOLVER& solver,
            const FIELDBINDER& fieldBinder,
            boost::function<
                typename detail_::ElementTypeBinder<FIELDTUPLEBINDER>::VecDof(
                    const typename detail_::ElementTypeBinder<FIELDTUPLEBINDER>::GlobalVecDim& )>
            forceFun )
        {
            typedef detail_::EvaluateDirectly<FIELDTUPLEBINDER> Evaluate;
            
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
                 typename SOLVER, typename FIELDBINDER>
        void bodyForceComputation2(
            const QUADRATURE& quadrature,
            SOLVER& solver,
            const FIELDBINDER& fieldBinder,
            boost::function<
                typename detail_::ElementTypeBinder<FIELDTUPLEBINDER>::VecDof(
                    const typename detail_::ElementTypeBinder<FIELDTUPLEBINDER>::GeomElement* , 
                    const typename detail_::ElementTypeBinder<FIELDTUPLEBINDER>::LocalVecDim& )>
            forceFun )
        {
            typedef detail_::EvaluateViaElement<FIELDTUPLEBINDER> Evaluate;
            
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
    typedef FIELDTUPLE FieldTuple;
    typedef EVALUTEPOLICY EvaluatePolicy;
    //@}

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

    //! How to evaluate the force function
    typedef typename EvaluatePolicy::Fun ForceFun;
    
    //--------------------------------------------------------------------------
    //! Constructor setting all references
    BodyForce( ForceFun& forceFun  )
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

            vector.segment( s * doFSize,
                            doFSize )  += f * funValues[s] * weight * detJ;
        }
                
        return;
        
    }
    
    //--------------------------------------------------------------------------
private:
    ForceFun&      forceFun_;  //!< Function producing body force
};


#endif

