//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ForceIntegrator.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_forceintegrator_hpp
#define base_asmb_forceintegrator_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
// boost includes
#include <boost/bind.hpp>
#include <boost/function.hpp>
// base includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{
        
        template<typename QUAD, typename SOLVER, typename FIELDTUPLE>
        class ForceIntegrator;

        //--------------------------------------------------------------------------
        // Convenience function to compute the residual forces
        template<typename QUADRATURE, typename SOLVER,
                 typename BOUNDFIELD, typename KERNEL>
        void computeResidualForces( const QUADRATURE& quadrature,
                                    SOLVER&           solver,
                                    const BOUNDFIELD& boundField,
                                    const KERNEL& kernelObj )
        {
            typedef ForceIntegrator<QUADRATURE,SOLVER,
                                    typename BOUNDFIELD::ElementPtrTuple>
                ForceIntegrator;

            typename ForceIntegrator::ForceKernel
                residualForce = boost::bind( &KERNEL::residualForce, 
                                             &kernelObj,
                                             _1, _2, _3, _4 );

            // Note the -1 for moving the forces to the RHS of the system
            ForceIntegrator forceInt( residualForce, quadrature, solver, -1.0 );

            // apply to all elements
            std::for_each( boundField.elementsBegin(), boundField.elementsEnd(),
                           forceInt );
        
            return;
        }

        
        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            // Helper function in order to reduce redundancies below
            template<typename TESTELEMENT>
            void collectFromDoFs( const TESTELEMENT*        testEp,
                                  std::vector<bool>&        doFActivity,
                                  std::vector<std::size_t>& doFIDs )
            {
                // Get the element's DoF objects
                std::vector<typename TESTELEMENT::DegreeOfFreedom*> doFs;
                std::copy( testEp -> doFsBegin(), testEp -> doFsEnd(),
                           std::back_inserter( doFs ) );

                // ask for activity and ID
                for ( unsigned d = 0; d < doFs.size(); d ++ ) {
                    doFs[d] -> getActivity( std::back_inserter( doFActivity ) ); 
                    doFs[d] -> getIndices(  std::back_inserter( doFIDs ) );
                }
        
            }

            //------------------------------------------------------------------
            // Helper function in order to reduce redundancies below
            template<typename SOLVER>
            void assembleForces( const base::VectorD&            forceVec,
                                 const std::vector<bool>&        doFActivity,
                                 const std::vector<std::size_t>& doFIDs,
                                 SOLVER& solver )
            {
                const unsigned numActiveDoFs =
                    std::count_if( doFActivity.begin(), doFActivity.end(),
                                   boost::bind( std::equal_to<bool>(), _1, true ) );

                // Result container
                base::VectorD sysVector = VectorD::Zero( numActiveDoFs );

                // Active ID numbers
                std::vector<std::size_t> activeDoFIDs( numActiveDoFs );

                // Collect for active DoFs
                unsigned ctr = 0;
                for ( unsigned d = 0; d < doFIDs.size(); d ++ ) {

                    if ( doFActivity[d] ) {

                        sysVector[    ctr ] = forceVec[ d ];
                        activeDoFIDs[ ctr ] = doFIDs[   d ];
                        ctr++;
                    }
                }

                // Pass on to system solver
                solver.insertToRHS( sysVector, activeDoFIDs );
            }


        } // namespace detail_
    } // namespace asmb    
} // namespace base

//------------------------------------------------------------------------------
/** Computation of the 'nodal' forces due to a given element  force function.
 *  The contribution of a force term to the system RHS comes from
 *  the linear form 
 *  \f[
 *        L( \phi ) =  \int_{\Omega} f(\phi) d x, 
 *  \f]
 *  where \f$ f \f$ is the force function. In the implementation, the evaluation
 *  of \f$ f \f$ will be represented by means of an element pointer and local
 *  quadrature point with weight. Finally, this term will be approximated with a
 *  given quadrature rule.
 *
 *  \tparam QUAD        Type of quadrature 
 *  \tparam SOLVER      Type of system solver with RHS storage
 *  \tparam FIELDTUPLE  Type of tuple holding field element pointers
 */
template<typename QUAD, typename SOLVER, typename FIELDTUPLE>
class base::asmb::ForceIntegrator
    : public boost::function<void( const FIELDTUPLE& )>
{
public:
    //! @name Template parameter
    //@{
    typedef QUAD        Quadrature;
    typedef SOLVER      Solver;
    typedef FIELDTUPLE  FieldTuple;
    //@}

    typedef typename FieldTuple::TestElementPtr   TestElementPtr;
    
    //! Type of forcing function
    typedef boost::function< void( const FieldTuple&,
                                   const typename Quadrature::VecDim&,
                                   const double,
                                   base::VectorD& )>  ForceKernel;

    //--------------------------------------------------------------------------
    //! Constructor setting all references
    ForceIntegrator( ForceKernel&         forceKernel,
                     const Quadrature&    quadrature,
                     Solver&              solver,
                     const double factor = 1.0 )
        : forceKernel_( forceKernel ),
          quadrature_(  quadrature ),
          solver_(      solver ),
          factor_(      factor )
    { }

    //--------------------------------------------------------------------------
    void operator()( const FieldTuple& fieldTuple )
    {
        // extract test and trial elements from tuple
        TestElementPtr  testEp  = fieldTuple.testElementPtr();
        
        // dof activities and IDs
        std::vector<bool> doFActivity;
        std::vector<std::size_t> doFIDs;
        detail_::collectFromDoFs( testEp, doFActivity, doFIDs );
        
        // Compute the element contribution to all its dofs
        base::VectorD forceVec = base::VectorD::Zero( doFIDs.size() );
        {
            typename Quadrature::Iter qIter = quadrature_.begin();
            typename Quadrature::Iter qEnd  = quadrature_.end();
            for ( ; qIter != qEnd; ++qIter ) {

                // Quadrature weight
                const double weight = qIter -> first;

                // Quadrature point
                const typename Quadrature::VecDim xi = qIter -> second;

                // add to force vector
                forceKernel_( fieldTuple, xi, weight, forceVec );
                
            } // end of quadrature
        }

        // apply factor
        forceVec *= factor_;
        
        // Assemble to solver
        detail_::assembleForces( forceVec, doFActivity, doFIDs, solver_ );
        
        return;
    }

private:
    ForceKernel&         forceKernel_; //!< Function producing the force term
    const Quadrature&    quadrature_;  //!< Quadrature
    Solver&              solver_;      //!< Solver with RHS storage
    const double         factor_;      //!< Scalar multiplier
};


#endif

