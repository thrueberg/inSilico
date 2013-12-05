//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ResidualForceHistory.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_time_residualforcehistory_hpp
#define base_time_residualforcehistory_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// base includes
#include <base/linearAlgebra.hpp>
#include <base/types.hpp>
// base/asmb includes
#include <base/asmb/ForceIntegrator.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace time{

        //----------------------------------------------------------------------
        template<typename KERNEL, typename QUAD, typename SOLVER, typename TSMETHOD>
        class ResidualForceHistory;

        //----------------------------------------------------------------------
        namespace detail_ {

            template<typename KERNEL, int CTR>
            struct CollectResidualForces
            {
                typedef typename KERNEL::FieldTuple       FieldTuple;
                typedef typename FieldTuple::TrialElement TrialElement;

                // number of abvailable history terms
                static const unsigned past =
                    TrialElement::DegreeOfFreedom::nHist - CTR;

                // size of individual dof
                static const unsigned doFSize =
                    TrialElement::DegreeOfFreedom::size;


                static void apply( const KERNEL&                       kernel,
                                   const FieldTuple&                   fieldTuple,
                                   const typename KERNEL::LocalVecDim& xi,
                                   const double                        quadWeight, 
                                   std::vector<double>&                weights,
                                   base::VectorD&                      result )
                {
                    // quit in case for no additional reaction terms
                    const unsigned numWeights = static_cast<unsigned>( weights.size() );
                    if ( numWeights == 0 ) return;

                    base::VectorD tmp = base::VectorD::Zero( result.size() );

                    // residual force caller
                    kernel.template residualForceHistory<past>( fieldTuple, xi,
                                                                quadWeight, tmp );

                    // weight by weights
                    result += weights[0] * tmp;

                    if ( numWeights > 1 ) {
                        // truncate weights
                        std::vector<double>::iterator first = weights.begin();
                        ++first;
                        weights = std::vector<double>( first, weights.end() );

                        // recursive call
                        CollectResidualForces<KERNEL,
                                              CTR-1>::apply( kernel, fieldTuple,
                                                             xi, quadWeight,
                                                             weights, result );
                    }
                }
            };

            template<typename KERNEL>
            struct CollectResidualForces<KERNEL,-1>
            {
                static void apply( const KERNEL&                       kernel,
                                   const typename KERNEL::FieldTuple&  fieldTuple,
                                   const typename KERNEL::LocalVecDim& xi,
                                   const double                        quadWeight, 
                                   std::vector<double>&                weights,
                                   base::VectorD&                      result )
                {
                    // empty in order to stop recursion
                }
            };

        }// namespace detail_
        
    }
}

//------------------------------------------------------------------------------
/** Computation of the reaction terms due to a time-integration process.
 *  Application of a linear multi-step method to the coupled system of ODEs
 *  \f[
 *        M \dot{y} + F(t,y) = 0
 *  \f]
 *  leads to the system of equations (in the \f$ i \f$-th Newton step)
 *  \f[
 *      \left[ \frac{a_0}{b_0 \Delta t} M + K_i \right]
 *      (x_{i+1}^{n+1} - x_i^{n+1}) =
 *      - \left[ \frac{a_0}{b_0 \Delta t} M x^{n+1}_i + F(x^{n+1}_i) \right]
 *      - M \left( \sum_{s=1}^A \frac{a_s}{b_0 \Delta t} x^{n+1-s} \right)
 *      - \sum_{s=1}^B \frac{b_s}{b_0} F(x^{n+1-s})
 *  \f]
 *
 */
template<typename KERNEL, typename QUAD, typename SOLVER, typename TSMETHOD>
class base::time::ResidualForceHistory
{
public:
    //! @name Template parameter
    //@{
    typedef KERNEL        Kernel;
    typedef QUAD          Quadrature;
    typedef SOLVER        Solver;
    typedef TSMETHOD      TimeSteppingMethod;
    //@}

    //! FieldTuple
    typedef typename Kernel::FieldTuple FieldTuple;

    typedef typename FieldTuple::TestElementPtr   TestElementPtr;
    typedef typename FieldTuple::TrialElement     TrialElement;

    // number of abvailable history terms
    static const unsigned nHist = TrialElement::DegreeOfFreedom::nHist;

    //! Constructor
    ResidualForceHistory( const Kernel&       kernel,
                          const Quadrature&    quadrature,
                          Solver&              solver,
                          const unsigned       step )
        : kernel_(     kernel ),
          quadrature_( quadrature ),
          solver_(     solver ),
          step_(       step )
    { }

    //--------------------------------------------------------------------------
    //! General case: possibly different test and trial spaces
    void operator()( const FieldTuple& fieldTuple )
    {
        // extract test and trial elements from tuple
        TestElementPtr  testEp  = fieldTuple.testElementPtr();

        // dof activities
        std::vector<bool> doFActivity;

        // dof IDs
        std::vector<std::size_t> doFIDs;

        // Collect dof entities from element
        base::asmb::detail_::collectFromDoFs( testEp, doFActivity, doFIDs );

        // RHS weights for reaction terms
        std::vector<double> forceWeights;
        TimeSteppingMethod::forceWeights( step_, forceWeights );

        // negative weights for moving the forces to the right hand side
        for ( unsigned w = 0; w < forceWeights.size(); w++ )
            forceWeights[w] *= -1.0;

        // Compute the element force contribution
        {
            base::VectorD elemForce = base::VectorD::Zero( doFIDs.size() );

            // do the quadrature loop
            typename Quadrature::Iter qIter = quadrature_.begin();
            typename Quadrature::Iter qEnd  = quadrature_.end();
            for ( ; qIter != qEnd; ++qIter ) {
                
                detail_::CollectResidualForces<Kernel,
                                               nHist-1>::apply( kernel_,
                                                                fieldTuple, 
                                                                qIter -> second,
                                                                qIter -> first,
                                                                forceWeights,
                                                                elemForce );
            }
            
            // Assemble to solver
            base::asmb::detail_::assembleForces( elemForce, doFActivity,
                                                 doFIDs, solver_ );
        }
        
        return;
    }

private:
    const Kernel&        kernel_;      //!< Kernel object
    const Quadrature&    quadrature_;  //!< Quadrature object
    Solver&              solver_;      //!< Solver object
    const unsigned       step_;        //!< Number of time step
};

#endif
