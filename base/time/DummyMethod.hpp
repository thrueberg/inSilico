//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   DummyMethod.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_time_dummymethod_hpp
#define base_time_dummymethod_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/verify.hpp>
// base/time includes
#include <base/time/MultiStep.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace time{

        template<unsigned DUMMY>
        class DummyMethod;

        namespace detail_{

            //! Compiler workaround
            template<unsigned Q>
            struct NumSteps<Q,DummyMethod>
            {
                static const unsigned numLHS   = 0;
                static const unsigned numRHS   = 1;
                static const unsigned numSteps = 0;
            };
        }

    }
}

//------------------------------------------------------------------------------
/** Time integration method which does not integrate in time.
 *  This method is useful if dynamic and static calculations are carried out
 *  in the same environment. Calling the functions
 *  base::time::computeInertiaTerms() and computeResidualForceHistory() with
 *  this method does not activate any time-integration terms.
 *  \tparam DUMMY Used for conformity, but not relevant
 */
template<unsigned DUMMY>
class base::time::DummyMethod
    : public base::time::MultiStep<1, base::time::DummyMethod >
{
public:
    //! Order of the method
    static const unsigned order = 1;

    //! For introspection
    static const bool isImplicit = true;

    //! Number of LHS terms needed
    static const unsigned numLHS = 1;

    //! Number of RHS terms needed
    static const unsigned numRHS = 1;

    typedef boost::array<double,1> LHSWeightsArray;
    typedef boost::array<double,1> RHSWeightsArray;


    //! @name Accessors
    //@{

    //! Return zero for \f$ a_0 \f$
    static LHSWeightsArray getLHSWeights()
    {
        RHSWeightsArray dummy = {{ 0.0 }};
        return dummy;
    }

    //! Return one for \f$ b_0 \f$
    static RHSWeightsArray getRHSWeights()
    {
        RHSWeightsArray dummy = {{ 1.0 }};
        return dummy;
    }
    //@}
};

//------------------------------------------------------------------------------
namespace base{
    namespace time{

        //! Use this, if you don't like an unused template parameter
        class NoTimeIntegration
            : public base::time::DummyMethod<1>
        { };

    }
}
#endif


