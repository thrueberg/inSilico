//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ReactionTerms.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_time_reactionterms_hpp
#define base_time_reactionterms_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// base includes
#include <base/linearAlgebra.hpp>
// base/asmb includes
#include <base/asmb/collectFromDoFs.hpp>
#include <base/asmb/assembleMatrix.hpp>
#include <base/asmb/assembleForces.hpp>

// base/kernel includes
#include <base/kernel/Mass.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace time{

        template<typename QUAD, typename SOLVER, typename TSMETHOD,
                 typename FIELDTUPLE>
        class ReactionTerms;

        //----------------------------------------------------------------------
        template<typename FIELDTUPLEBINDER, typename MSM, typename KERNEL,
                 typename QUADRATURE, typename SOLVER, typename FIELDBINDER>
        void computeReactionTerms( const KERNEL& kernel, 
                                   const QUADRATURE& quadrature,
                                   SOLVER&           solver,
                                   const FIELDBINDER& fieldBinder,
                                   const double stepSize,
                                   const unsigned step, 
                                   const bool zeroConstraints = false )
        {
            typedef typename FIELDTUPLEBINDER::Tuple ElementPtrTuple;
            typedef base::time::ReactionTerms<QUADRATURE,SOLVER,MSM,ElementPtrTuple> RT;
            typename RT::Kernel kernelFun = boost::bind( &KERNEL::tangentStiffness,
                                                         &kernel, _1, _2, _3, _4 );

            RT rt( kernelFun, quadrature, solver, stepSize, step, zeroConstraints );

            // Apply to all elements
            typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
            typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
            for ( ; iter != end; ++iter ) {
                rt( FIELDTUPLEBINDER::makeTuple( *iter ) );
            }

        }

        //----------------------------------------------------------------------
        template<typename FIELDTUPLEBINDER, typename MSM,
                 typename QUADRATURE, typename SOLVER, typename FIELDBINDER>
        void computeInertiaTerms( const QUADRATURE& quadrature,
                                  SOLVER&           solver,
                                  const FIELDBINDER& fieldBinder,
                                  const double stepSize,
                                  const unsigned step, 
                                  const double density,
                                  const bool zeroConstraints = false )
        {
            typedef typename FIELDTUPLEBINDER::Tuple ElementPtrTuple;
            
            base::kernel::Mass<ElementPtrTuple> mass( density );
            typedef base::time::ReactionTerms<QUADRATURE,SOLVER,MSM,ElementPtrTuple> RT;
            typename RT::Kernel kernelFun = boost::bind( mass, _1, _2, _3, _4 );

            RT rt( kernelFun, quadrature, solver, stepSize, step, zeroConstraints );

            // Apply to all elements
            typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
            typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
            for ( ; iter != end; ++iter ) {
                rt( FIELDTUPLEBINDER::makeTuple( *iter ) );
            }

        }

                                   

        //----------------------------------------------------------------------
        namespace detail_{

            template<typename TRIALELEMENT, int CTR>
            struct CollectWeightedHistory
            {
                // which history to access
                static const unsigned past =
                    TRIALELEMENT::DegreeOfFreedom::nHist - CTR;// + 1;

                // size of individual dof
                static const unsigned doFSize =
                    TRIALELEMENT::DegreeOfFreedom::size;

                static void apply( const TRIALELEMENT* trialEp,
                                   std::vector<double>& weights,
                                   base::VectorD& result )
                {
                    
                    // quit in case for no additional reaction terms
                    const unsigned numWeights = static_cast<unsigned>( weights.size() );
                    if ( numWeights == 0 ) return;

                    base::VectorD tmp = base::VectorD::Zero( result.size() );
                    
                    // collect history result from element
                    typename TRIALELEMENT::DoFPtrConstIter doFIter =
                        trialEp -> doFsBegin();
                    typename TRIALELEMENT::DoFPtrConstIter doFEnd  =
                        trialEp -> doFsEnd();
                    for ( unsigned d = 0; doFIter != doFEnd; d++, ++doFIter ) {
                        for ( unsigned s = 0; s < doFSize; s++)
                            tmp[d * doFSize + s]
                                = (*doFIter) ->template getHistoryValue<past>( s );
                    }

                    // weight by weights
                    result += weights[0] * tmp;

                    if ( numWeights > 1 ) {
                        // truncate weights
                        std::vector<double>::iterator first = weights.begin();
                        ++first;
                        weights = std::vector<double>( first, weights.end() );

                        // recursive call
                        CollectWeightedHistory<TRIALELEMENT,CTR-1>::apply( trialEp,
                                                                           weights,
                                                                           result );
                    }
                }
            };

            // specialisation for 
            template<typename TRIALELEMENT>
            struct CollectWeightedHistory<TRIALELEMENT,-1>
            {
                static void apply( const TRIALELEMENT* trialEp,
                                   std::vector<double>& weights,
                                   base::VectorD& result )
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
 *  leads to the system of equations (in the \f$i \f$-th Newton step)
 *  \f[
 *      \left[ \frac{a_0}{b_0 \Delta t} M + K_i \right]
 *      (x_{i+1}^{n+1} - x_i^{n+1}) =
 *      - \left[ \frac{a_0}{b_0 \Delta t} M x^{n+1}_i + F(x^{n+1}_i) \right]
 *      - M \left( \sum_{s=1}^A \frac{a_s}{b_0 \Delta t} x^{n+1-s} \right)
 *      - \sum_{s=1}^B \frac{b_s}{b_0} F(x^{n+1-s})
 *  \f]
 *  This class computes and assembles all terms containing the matrix
 *  \f$ M \f$. By passing a flag, it is decided whether the the latest iterate
 *  \f$ x^{n+1}_i \f$ contributes to the system or not. In most applications,
 *  \f$ M \f$ is the so-called mass matrix with defining bilinear form
 *  \f[
 *       m(u,v) = \int_\Omega \rho u v d x
 *  \f]
 *  The default behaviour of this object is to use  base::kernel::Mass and 
 *  for this reason a constructor with the mass density \f$ \rho \f$ is
 *  provided. Alternatively, the user can provide a custom kernel function
 *  which defines \f$ M \f$.
 */
template<typename QUAD, typename SOLVER, typename TSMETHOD,
         typename FIELDTUPLE>
class base::time::ReactionTerms
{
public:
    //! @name Template parameter
    //@{
    typedef QUAD          Quadrature;
    typedef SOLVER        Solver;
    typedef TSMETHOD      TimeSteppingMethod;
    typedef FIELDTUPLE    FieldTuple;
    //@}

    //! @name Access types of the tuple
    //@{
    typedef typename FieldTuple::TestElementPtr   TestElementPtr;
    typedef typename FieldTuple::TrialElementPtr  TrialElementPtr;
    typedef typename FieldTuple::TrialElement     TrialElement;
    //@}

    //! Bubnov-Galerkin method
    static const bool isBubnov =
        boost::is_same<TestElementPtr,TrialElementPtr>::value;

    //! General Kernel function
    typedef boost::function<void( const FieldTuple&,
                                  const typename Quadrature::VecDim&,
                                  const double,
                                  base::MatrixD& ) >  Kernel;

    //! Class for mass matrix computation
    typedef base::kernel::Mass<FieldTuple> Mass;
    
    //! Constructor with use-provided kernel function
    ReactionTerms( const Kernel&        kernel, 
                   const Quadrature&    quadrature,
                   Solver&              solver,
                   const double         stepSize, 
                   const unsigned       step,
                   const bool           prevIterate = false )
        : mass_(            base::invalidReal() ), // not used
          kernel_(          kernel ), 
          quadrature_(      quadrature ),
          solver_(          solver ),
          stepSize_(        stepSize ),
          step_(            step ),
          prevIterate_(     prevIterate )
    { }

    //--------------------------------------------------------------------------
    //! General case: possibly different test and trial spaces
    void operator()( const FieldTuple&   fieldTuple )
    {
        // extract test and trial elements from tuple
        TestElementPtr  testEp  = fieldTuple.testElementPtr();
        TrialElementPtr trialEp = fieldTuple.trialElementPtr();

        // dof activities
        std::vector<bool> rowDoFActivity, colDoFActivity;

        // dof IDs
        std::vector<std::size_t> rowDoFIDs, colDoFIDs;

        // dof values (for constraints)
        std::vector<number> rowDoFValues, colDoFValues;

        // dof constraints
        typedef std::pair<unsigned, std::vector< std::pair<base::number,std::size_t> > >
            WeightedDoFIDs;
        std::vector<WeightedDoFIDs> rowConstraints, colConstraints;

                // Collect dof entities from element
        base::asmb::collectFromDoFs( testEp, rowDoFActivity,
                                     rowDoFIDs, rowDoFValues,
                                     rowConstraints );

        if ( isBubnov ) {
            colDoFActivity = rowDoFActivity;
            colDoFIDs      = rowDoFIDs;
            colDoFValues   = rowDoFValues;
            colConstraints = rowConstraints;
        }
        else
            base::asmb::collectFromDoFs( trialEp, colDoFActivity,
                                         colDoFIDs, colDoFValues,
                                         colConstraints );

        // Compute the element matrix contribution
        base::MatrixD elemMat = base::MatrixD::Zero( rowDoFIDs.size(),
                                                     colDoFIDs.size() );

        // apply quadrature
        quadrature_.apply( kernel_, fieldTuple, elemMat );

        // LHS
        {
            // LHS weight
            const double alpha = TimeSteppingMethod::systemMassWeight( step_ );

            const base::MatrixD lhsMatrix = (alpha/stepSize_) * elemMat;
        
            // assemble element matrix to global system
            base::asmb::assembleMatrix( lhsMatrix,
                                        rowDoFActivity, colDoFActivity,
                                        rowDoFIDs, colDoFIDs,
                                        colDoFValues,
                                        rowConstraints, colConstraints, 
                                        solver_, isBubnov,
                                        prevIterate_ );
        }

        // RHS
        {
            // RHS weights for reaction terms
            std::vector<double> reactionWeights;
            TimeSteppingMethod::reactionWeights( step_, reactionWeights );

            // divide weights by step size
            for ( unsigned s = 0; s < reactionWeights.size(); s++ )
                reactionWeights[s] /= -stepSize_;

            // container for weighted history of solutions
            base::VectorD resultVec = base::VectorD::Zero( colDoFIDs.size() );

            // call solution collector (recursive)
            if ( prevIterate_ ) {
                // include terms of the previous iterate

                // call recursion beginning with the current iterate
                detail_::CollectWeightedHistory<
                    TrialElement,
                    TrialElement::DegreeOfFreedom::nHist>::apply( trialEp,
                                                                  reactionWeights,
                                                                  resultVec );
            }
            else {

                // exclude the leading term for the previous iterate

                // remove first element from the list of weights
                std::vector<double>::iterator beginWeights = reactionWeights.begin();
                beginWeights++;
                reactionWeights = std::vector<double>( beginWeights, reactionWeights.end() );

                // call recursion beginning with the previous solution
                detail_::CollectWeightedHistory<
                    TrialElement,
                    TrialElement::DegreeOfFreedom::nHist-1>::apply( trialEp,
                                                                    reactionWeights,
                                                                    resultVec );
            }

            // multiply by mass matrix
            base::VectorD forceVec;
            forceVec.noalias() = elemMat * resultVec;

            // Assemble to solver
            base::asmb::assembleForces( forceVec, rowDoFActivity,
                                        rowDoFIDs, rowConstraints, solver_ );
        }
        
        return;
    }


private:
    const Mass           mass_;        //!< Mass matrix provider
    const Kernel         kernel_;      //!< Kernel function
    const Quadrature&    quadrature_;  //!< Quadrature object
    Solver&              solver_;      //!< Solver object

    const double         stepSize_;    //!< Time step size
    const unsigned       step_;        //!< Number of time step
    const bool           prevIterate_; //!< There is a previous iterate
};

#endif
