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
    Distribute( const SourceOfValues & sourceOfValues,
                const bool zeroConstraints = false )
        : sourceOfValues_(  sourceOfValues ),
          zeroConstraints_( zeroConstraints )
    {  }

    //! Function call operator acting on an individual dof
    void operator()( DOF * dof ) const
    {
        // Collect dof IDs
        std::vector<unsigned> dofIDs( dofSize );
        dof -> getIndices( dofIDs.begin() );

        // Collect activity flags
        std::vector<bool>     dofActivity( dofSize );
        dof -> getActivity( dofActivity.begin() );

        // Set values of active dofs
        for ( unsigned d = 0; d < dofSize; d ++ ) {
            
            if ( dofActivity[d] ) {
                
                dof -> setValue( d,
                                 DoFOp::apply( dof -> getValue( d ), 
                                               sourceOfValues_.getValue( dofIDs[d] ) ) );
            }
            else {
                if ( not zeroConstraints_ )
                    dof -> setValue( d,
                                     DoFOp::apply( dof -> getValue( d ),
                                                   dof -> getConstraint( d ) ) );
            }
                                               
        }

        return;
    }
    
private:
    //! Access to source of values (usually the solver object)
    const SourceOfValues & sourceOfValues_;

    const bool zeroConstraints_;
};

#endif
