//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Distribute.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_distribute_hpp
#define base_distribute_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// base includes
#include <base/numbers.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        enum DoFOperation
        {
            CLEAR = 0, //!< Set to zero
            SET,       //!< Set to new value 
            ADD        //!< Add new value
        };

        template<typename DOF, typename SRC, DoFOperation DOFOP>
        class Distribute;

        //----------------------------------------------------------------------
        //! Convenience function to set dof values
        template<typename SOLVER, typename FIELD>
        void setDoFsFromSolver( const SOLVER& solver, FIELD& field )
        {
            base::dof::Distribute<typename FIELD::DegreeOfFreedom,
                                  SOLVER,SET> distributeDoF( solver );
            distributeDoF.apply( field.doFsBegin(), field.doFsEnd() );
        }

        //----------------------------------------------------------------------
        //! Convenience function to add to dof values
        template<typename SOLVER, typename FIELD>
        void addToDoFsFromSolver( const SOLVER& solver, FIELD& field )
        {
            base::dof::Distribute<typename FIELD::DegreeOfFreedom,SOLVER,ADD>
                distributeDoF( solver );
            distributeDoF.apply( field.doFsBegin(), field.doFsEnd() );
        }

        //----------------------------------------------------------------------
        //! Convenience function to clear the current dof values
        template<typename FIELD>
        void clearDoFs( FIELD& field )
        {
            typename FIELD::DoFPtrIter dIter = field.doFsBegin();
            typename FIELD::DoFPtrIter dEnd  = field.doFsEnd();
            for ( ; dIter != dEnd; ++dIter ) (*dIter) -> clearValue();
        }

        //----------------------------------------------------------------------
        //! Convenience function to clear the current dof values
        template<typename FIELD>
        void pushHistory( FIELD& field )
        {
            typename FIELD::DoFPtrIter dIter = field.doFsBegin();
            typename FIELD::DoFPtrIter dEnd  = field.doFsEnd();
            for ( ; dIter != dEnd; ++dIter ) (*dIter) -> pushHistory();
        }
        
        //----------------------------------------------------------------------
        //! Convenience function to multiply existing constraints by a factor
        template<typename FIELD>
        void scaleConstraints( FIELD& field, const double factor )
        {
            typename FIELD::DoFPtrIter dIter = field.doFsBegin();
            typename FIELD::DoFPtrIter dEnd  = field.doFsEnd();
            for ( ; dIter != dEnd; ++dIter ) (*dIter) -> scaleConstraint( factor );
        }

        //----------------------------------------------------------------------
        template<typename FIELD>
        void clearConstraints( FIELD& field )
        {
            typename FIELD::DoFPtrIter dIter = field.doFsBegin();
            typename FIELD::DoFPtrIter dEnd  = field.doFsEnd();
            for ( ; dIter != dEnd; ++dIter ) (*dIter) -> clearConstraints();
        }
        
        //----------------------------------------------------------------------
        namespace detail_{

            //! Type of operation  c = op( a, b )
            template<DoFOperation DOFOP> struct DoFManip;

            template<>
            struct DoFManip<CLEAR>
            {
                template<typename T>
                static T apply( const T& a, const T& b ) { return T(); }
            };

            template<>
            struct DoFManip<SET>
            {
                template<typename T>
                static T apply( const T& a, const T& b ) { return b; }
            };

            template<>
            struct DoFManip<ADD>
            {
                template<typename T>
                static T apply( const T& a, const T& b ) { return a+b; }
            };
            
        }
    }
}

//------------------------------------------------------------------------------
/** Pass values to degree of freedom.
 *  The main activity of this object is carried out in 2 steps:
 *
 *    1. Ask every active dof-component for its ID, get the value from the
 *       source (usually the solver), and apply that value to the dof-component
 *       (according to the DOFOP flag by setting or adding)
 *    2. Evaluate the constraint of every in-active dof (either as a homogeneous
 *       dof, in case of e.g. non-linear iterations, or fully) and apply that
 *       value to the dof-component (set, not add!)
 *
 *
 *  \tparam DOF   Type of degree of freedom
 *  \tparam SRC   Source of values (normally the solver)
 *  \tparam DOFOP Type of operation on the dof values
 */
template<typename DOF, typename SRC,
         base::dof::DoFOperation DOFOP = base::dof::SET>
class base::dof::Distribute
{
public:

    //! @name Template parameters
    //@{
    typedef  DOF DegreeOfFreedom;
    typedef  SRC SourceOfValues; //!< Solver normally
    //@}

    //! Size of the individual dof
    static const unsigned dofSize = DegreeOfFreedom::size;

    //! operation on the dofs
    typedef detail_::DoFManip<DOFOP>  DoFOp;

    //! Constructor with source access
    Distribute( const SourceOfValues & sourceOfValues )
        : sourceOfValues_(  sourceOfValues )
    {  }

    //! Distribute the solution values back to the degrees of freedom
    template<typename DOFITER>
    void apply( DOFITER first, DOFITER last )
    {
        // go through all dofs in order to collect solution values
        for ( DOFITER iter = first; iter != last; ++iter ) {

            DegreeOfFreedom* dof = *iter;
            
            // Collect dof IDs
            std::vector<unsigned> dofIDs( dofSize );
            dof -> getIndices( dofIDs.begin() );

            // Go through all components
            for ( unsigned d = 0; d < dofSize; d ++ ) {

                // only query solution for active dofs
                if ( dof -> isActive(d) ) {
                
                    dof -> setValue( d,
                                     DoFOp::apply( dof -> getValue( d ), 
                                                   sourceOfValues_.getValue( dofIDs[d] ) ) );
                }
                
            }
            
        }

        // go through all dofs in order to evaluate constraints
        for ( DOFITER iter = first; iter != last; ++iter ) {

            DegreeOfFreedom* dof = *iter;

            for ( unsigned d = 0; d < dofSize; d++ ) {

                // only evaluate constraints at inactive dofs
                if ( dof -> isConstrained(d) ) {

                    // evaluate the constraint of this dof-component
                    const base::number constraintValue =
                        dof -> getConstraint( d ) -> evaluate( );

                    // apply to dof
                    dof -> setValue( d, constraintValue );
                }
            }
            
        }
    }
    
private:
    //! Access to source of values (usually the solver object)
    const SourceOfValues & sourceOfValues_;
};

#endif
