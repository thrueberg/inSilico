//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   AdamsMoulton.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_time_adamsmoulton_hpp
#define base_time_adamsmoulton_hpp

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
        template<unsigned Q> class AdamsMoulton;
    }
}

//------------------------------------------------------------------------------
/** Adams-Moulton method weights.
 *  AM methods operating on an ODE of the type
 *  \f[
 *         \dot{y} = f(t,y)
 *  \f]
 *  are derived by approximating the integral form of this ODE.
 *  Therefore they assume the general form for \f$ q \f$-th order
 *  \f[
 *      \frac{1}{\Delta t} (y^{n+1} - y^n) =
 *      \sum_{s=0}^q b_s f(t^{n+1-s}, y^{n+1-s})
 *  \f]
 *  In the framework of a general linear multi-step method, they LHS weights
 *  assume the values \f$ a_0 = 1 \f$, \f$ a_1 = -1 \f$ and \f$ a_s \f$ else.
 *
 *  Table of weights:
 *
 *  |  q  |  b0   |   b1  |  b2  |  b3  |
 *  |-----|-------|-------|------|------|
 *  |  1  |   1   |       |      |      |
 *  |  2  |  1/2  |  1/2  |      |      |
 *  |  3  | 5/12  |  2/3  |-1/12 |      |
 *  |  4  |  3/8  | 19/24 |-5/24 | 1/24 |
 *
 *  Note that \f$ b_q = 0 \f$ and therefore the method can be implemented
 *  as a \f$ q-1 \f$-step method, exception \f$ q = 1 \f$ yields still a
 *  1-step method because of the LHS terms.
 *
 *
 *  \sa http://en.wikipedia.org/wiki/Adams-Moulton_method
 *  \tparam Q  The methods order \f$ q \f$
 */
template<unsigned Q>
class base::time::AdamsMoulton
    : public base::time::MultiStep<Q, base::time::AdamsMoulton >
{
public:
    //! Template parameter: order of the method
    static const unsigned order = Q;

    //! For introspection
    static const bool isImplicit = true;

    //! Number of LHS terms needed
    static const unsigned numLHS = 2;

    //! Number of RHS terms needed
    static const unsigned numRHS = order;

    //! Array of RHS weights
    typedef boost::array<double,numLHS> LHSWeightsArray;

    //! Array of RHS weights
    typedef boost::array<double,numRHS> RHSWeightsArray;

    //! @name Accessors
    //@{

    /** LHS weights are always \f$ a_0 = 1\f$, \f$ a_1 = -1 \f$
     *  and \f$ a_s = 0 \f$ else
     */
    static LHSWeightsArray getLHSWeights()
    {
        LHSWeightsArray a = {{ 1.0, -1.0 }};
        return a;
    }

    //! RHS weights 
    static RHSWeightsArray getRHSWeights() { return b_; }
    //@}
    
private:
    static const RHSWeightsArray b_;
};

//------------------------------------------------------------------------------
//! \cond SKIPDOX
//------------------------------------------------------------------------------
// Specialisation via initialisation of static members
namespace base{
    namespace time{

        //----------------------------------------------------------------------
        // Q=1
        template<>
        const typename AdamsMoulton<1>::RHSWeightsArray AdamsMoulton<1>::b_ =
        {{ 1.0 }};

        //----------------------------------------------------------------------
        // Q=2
        template<>
        const typename AdamsMoulton<2>::RHSWeightsArray AdamsMoulton<2>::b_ =
        {{ 0.5, 0.5 }};
        
        //----------------------------------------------------------------------
        // Q=3
        template<>
        const typename AdamsMoulton<3>::RHSWeightsArray AdamsMoulton<3>::b_ =
        {{ 5./12., 2./3., -1./12 }};

        //----------------------------------------------------------------------
        // Q=4
        template<>
        const typename AdamsMoulton<4>::RHSWeightsArray AdamsMoulton<4>::b_ =
        {{ 3./8., 19./24., -5./24., 1./24. }};

    }
}
//! \endcond

#endif


