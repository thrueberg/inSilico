//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   MultiStep.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_time_multistep_hpp
#define base_time_multistep_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
// boost includes
#include <boost/integer/static_min_max.hpp>


//------------------------------------------------------------------------------
namespace base{
    namespace time{

        template<unsigned ORDER, template<unsigned> class METHOD>
        class MultiStep;

        namespace detail_{

            //------------------------------------------------------------------
            // Recursive call for LHS weights
            template<unsigned ORDER, template<unsigned> class METHOD>
            struct LHSWeights
            {
                static std::vector<double> apply( const unsigned step )
                {
                    std::vector<double> result;

                    // refer to lower order method
                    if ( step+2 < METHOD<ORDER>::numLHS ) { 
                        result = LHSWeights<ORDER-1,METHOD>::apply( step );
                    }
                    else {
                        const typename METHOD<ORDER>::LHSWeightsArray aux = 
                            METHOD<ORDER>::getLHSWeights();
                        
                        std::copy( aux.begin(), aux.end(), // copy to dynamic
                                   std::back_inserter( result ) );
                    }
                    
                    return result;
                }
            };

            // Specialisation for ORDER=1 (recursion stop)
            template<template<unsigned> class METHOD>
            struct LHSWeights<1,METHOD>
            {
                static std::vector<double> apply( const unsigned step )
                {
                    std::vector<double> result;
                    const typename METHOD<1>::LHSWeightsArray aux = 
                            METHOD<1>::getLHSWeights();
                    std::copy( aux.begin(), aux.end(), // copy to dynamic
                               std::back_inserter( result ) );
                    return result;
                }
            };
            
            //------------------------------------------------------------------
            // Recursive call for RHS weights
            template<unsigned ORDER, template<unsigned> class METHOD>
            struct RHSWeights
            {
                static std::vector<double> apply( const unsigned step )
                {
                    std::vector<double> result;

                    // refer to lower-order method
                    if ( step+2 < METHOD<ORDER>::numRHS ) { 
                        result = RHSWeights<ORDER-1,METHOD>::apply( step );
                    }
                    else {
                        const typename METHOD<ORDER>::RHSWeightsArray aux = 
                            METHOD<ORDER>::getRHSWeights();

                        std::copy( aux.begin(), aux.end(), // copy to dynamic
                                   std::back_inserter( result ) );
                    }
                    return result;
                }
            };

            // Specialisation for ORDER=1 (recursion stop)
            template<template<unsigned> class METHOD>
            struct RHSWeights<1,METHOD>
            {
                static std::vector<double> apply( const unsigned step )
                {
                    std::vector<double> result;
                    const typename METHOD<1>::RHSWeightsArray aux = 
                            METHOD<1>::getRHSWeights();
                    std::copy( aux.begin(), aux.end(), // copy to dynamic
                               std::back_inserter( result ) );
                    return result;
                }
            };

        } // namespace detail_
    }
}
//------------------------------------------------------------------------------
/** Interface for a linear multi-step method.
 *  The application of a linear multi-step method to the generic ODE
 *  \f[
 *        \dot{y} = f( t, y )
 *  \f]
 *  leads to the algebraic expression
 *  \f[
 *        \sum_{s=0}^A a_s y^{n+1-s} =
 *                 \Delta t \sum_{s=0}^B b_s f(t^{n+1-s}, y^{n+1-s})
 *  \f]
 *  where the coefficients \f$ a_s \f$ and \f$ b_s \f$ are specific to the
 *  chosen method, e.g., base::time::BDF or base::time::AdamsMoulton.
 *  In this context, a finite discretisation usually leads to a system of
 *  equations of the form
 *  \f[
 *       M \dot{x} + F(x) = 0
 *  \f]
 *  with the mass matrix \f$ M \f$, the vector of unknowns \f$ x \f$ and the
 *  vector of internal forces \f$ F \f$. An application of above multi-step
 *  method to this system leads to a Newton iteration of the form
 *  \f[
 *      \left[ \frac{a_0}{b_0 \Delta t} M + K_i \right]
 *      (x_{i+1}^{n+1} - x_i^{n+1}) =
 *      - \left[ \frac{a_0}{b_0 \Delta t} M x^{n+1}_i + F(x^{n+1}_i) \right]
 *      - M \left( \sum_{s=1}^A \frac{a_s}{b_0 \Delta t} x^{n+1-s} \right)
 *      - \sum_{s=1}^B \frac{b_s}{b_0} F(x^{n+1-s})
 *  \f]
 *  Here the sub-script \f$ i \f$ refers to the Newton iteration. Note that
 *  in the static case, the iteration would have the form
 *  \f[
 *    K_i (x_{i+1} - x_i) = - F(x_i)
 *  \f]
 *  This interface provides the weights needed for the extra terms which
 *  occur in the time-stepping but not in the static case. These are the
 *  weight of the system mass matrix \f$ M \f$, the weights of the right
 *  hand side reaction terms and the weights for the right hand side force
 *  history:
 *   - systemMassWeight() gives
 *                       \f$ a_0 /b_0 \f$
 *   - reactionWeights()  gives
 *                       \f$ \left\{ a_s/ b_0 \right\}_{s=0}^A \f$
 *   - forceWeights()     gives
 *                       \f$ \left\{ b_s / b_0 \right\}_{s=1}^B \f$
 *
 *  \par Startup
 *  Since higher-order methods require a solution history beyond the
 *  latest results, they need a startup phase. That means for the first
 *  steps, lower-order methods are used in order to build up the solution
 *  history. For example, a BDF-3 method begins with one step of BDF-1,
 *  a second step of BDF-2 and then every following step has the necessary
 *  history available for the complete BDF-3. In order to achieve this
 *  mechanism, recursive calls to lower-order methods are used.
 *
 *  \par Note on the implementation
 *  This class provides an interface to a specific implementation of
 *  a linear multi-step method. Therefore, these implementations inherit
 *  from this class. This is achieved by a variant of the curiously
 *  recurring template pattern
 *  (http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern).
 *  Other than
 *  \code
 *  template<typename METHOD> class MultiStep{ };
 *  \endcode
 *  it is preferred to have the ORDER and the METHOD as seperate
 *  template parameters.
 *  \code
 *  template<unsigned ORDER,template<unsigned> class METHOD> class MultiStep{ };
 *  \endcode
 *  This allows the interface to construct lower
 *  order methods for the startup phase as explained above.
 *
 *  \tparam ORDER  The order of the method
 *  \tparam METHOD The specific method (inherits from this interface)
 */
template<unsigned ORDER, template<unsigned> class METHOD>
class base::time::MultiStep
{
public:

    //! @name Template parameters
    //@{
    //! Order of the method
    static const unsigned order = ORDER;
    //! The actual implementation
    typedef METHOD<order>         Method;
    //@}

    //! Number of steps required deduced from number of LHS and RHS terms
    static const unsigned numSteps =
        boost::static_unsigned_max<Method::numLHS,Method::numRHS>::value-1;
    
    //! Weight for the system mass matrix
    static double systemMassWeight( const unsigned step )
    {
        VERIFY_MSG( Method::isImplicit,
                    "Illegal call for explicit time stepping" );
        
        std::vector<double> lhs = detail_::LHSWeights<order,METHOD>::apply( step );
        std::vector<double> rhs = detail_::RHSWeights<order,METHOD>::apply( step );

        return lhs[0] / rhs[0];
    }

    //! Weights for the reaction terms
    static void reactionWeights( const unsigned step,
                                 std::vector<double>& weights )
    {
        std::vector<double> lhs = detail_::LHSWeights<order,METHOD>::apply( step );
        std::vector<double> rhs = detail_::RHSWeights<order,METHOD>::apply( step );

        weights.resize( lhs.size() );
        for ( unsigned s = 0; s < weights.size(); s++ )
            weights[s] = lhs[s] / rhs[0];

        return;
    }

    //! Weights for the force terms
    static void forceWeights( const unsigned step,
                              std::vector<double>& weights )
    {
        std::vector<double> rhs = detail_::RHSWeights<order,METHOD>::apply( step );

        weights.resize( rhs.size() - 1 );
        for ( unsigned s = 0; s < weights.size(); s++ )
            weights[s] = rhs[s+1] / rhs[0];

        return;
    }

    //! Weights for the derivative computation
    static void derivativeWeights( const unsigned step,
                                   std::vector<double>& weights )
    {
        weights = detail_::LHSWeights<order,METHOD>::apply( step );
    }

};

#endif
