//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   BDF.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_time_bdf_hpp
#define base_time_bdf_hpp

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
        template<unsigned Q> class BDF;

        namespace detail_{

            //! Workaround for compilers 
            template<unsigned Q>
            struct NumSteps<Q,BDF>
            {
                static const unsigned numLHS = 2;
                static const unsigned numRHS = Q;
                static const unsigned numSteps =
                    boost::static_unsigned_max<numLHS,numRHS>::value-1;
            };
            
        }
    }
}

//------------------------------------------------------------------------------
/** Backward Differentiation Formula (%BDF) weights.
 *  BDF rules operating on an ODE of the type
 *  \f[
 *         \dot{y} = f(t,y)
 *  \f]
 *  are derived by using polynomial approximations of the derivative term and
 *  evaluating the force term always at the new time step. Therefore they
 *  assume the general form for \f$ q \f$-th order
 *  \f[
 *      \sum_{s=0}^q a_s y^{n+1-s} = \Delta t f( t^{n+1}, y^{n+1} )
 *  \f]
 *  In the framework of a general linear multi-step method, the RHS weights
 *  assume the values \f$ b_0 = 1 \f$ and \f$ b_s \f$ else.
 *
 *  Table of weights:
 *
 *  |  q  |  a0   |   a1  |  a2  |  a3  |  a4  |
 *  |-----|-------|-------|------|------|------|
 *  |  1  |   1   | -1    |      |      |      |
 *  |  2  |  3/2  | -2    |  1/2 |      |      |
 *  |  3  | 11/6  | -3    |  3/2 | -1/3 |      |
 *  |  4  | 25/12 | -4    |   3  | -4/3 |  1/4 |
 *
 *  For any order \f$ q \f$, \f$ q \f$ LHS terms \f$ y^{n+1-s} \f$ are needed
 *  and therefore the method is implemented as a \f$ q \f$-step method.
 *
 *  \sa http://en.wikipedia.org/wiki/Backward_differentiation_formula
 *  \tparam Q  The methods order \f$ q \f$
 */
template<unsigned Q>
class base::time::BDF
    : public base::time::MultiStep< Q, base::time::BDF >
{
public:
    //! Template parameter: order of the method
    static const unsigned order = Q;

    //! For introspection
    static const bool isImplicit = true;

    //! Number of LHS terms needed
    static const unsigned numLHS = order+1;

    //! Number of RHS terms needed
    static const unsigned numRHS = 1;

    //! Number of total steps required
    static const unsigned numSteps =
        boost::static_unsigned_max<numLHS,numRHS>::value-1;

    //! Array of RHS weights
    typedef boost::array<double,numLHS> LHSWeightsArray;

    //! Array of RHS weights
    typedef boost::array<double,numRHS> RHSWeightsArray;

    //! @name Accessors
    //@{

    //! LHS weights, 
    static LHSWeightsArray getLHSWeights( ) { return a_; }

    //! RHS weights are always \f$ b_0 = 1\f$ and \f$ b_s = 0\f$ else
    static RHSWeightsArray getRHSWeights( )
    {
        RHSWeightsArray b = {{ 1.0 }};
        return b;
    }
    //@}
    
private:
    static const LHSWeightsArray a_;
};

//------------------------------------------------------------------------------
// \cond SKIPDOX
//------------------------------------------------------------------------------
// Specialisation via initialisation of static members
namespace base{
    namespace time{

        //----------------------------------------------------------------------
        // Q=1
        template<> const BDF<1>::LHSWeightsArray BDF<1>::a_ = {{ 1.0, -1.0 }};

        //----------------------------------------------------------------------
        // Q=2
        template<> const BDF<2>::LHSWeightsArray BDF<2>::a_ =
        {{ 1.5, -2.0, 0.5 }};
        
        //----------------------------------------------------------------------
        // Q=3
        template<> const BDF<3>::LHSWeightsArray BDF<3>::a_ =
        {{ 11./6., -3.0, 1.5, -1./3. }};

        //----------------------------------------------------------------------
        // Q=4
        template<> const BDF<4>::LHSWeightsArray BDF<4>::a_ =
        {{ 25./12., -4.0, 3.0, 4./3., 1./4.}};

    }
}

//! \endcond

#endif


