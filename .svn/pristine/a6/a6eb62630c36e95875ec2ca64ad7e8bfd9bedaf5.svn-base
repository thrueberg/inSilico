//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   assembleForces.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_asmb_assembleforces_hpp
#define base_asmb_assembleforces_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
// boost includes
#include <boost/bind.hpp>
// base includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

            //------------------------------------------------------------------
            template<typename SOLVER>
            void assembleForces( const base::VectorD&            forceVec,
                                 const std::vector<base::dof::DoFStatus>& doFStatus,
                                 const std::vector<std::size_t>& doFIDs,
                                 const std::vector<
                                     std::pair<unsigned,
                                               std::vector<std::pair<base::number,
                                                                     std::size_t> > > > &
                                 constraints,
                                 SOLVER& solver );
    }
}

//------------------------------------------------------------------------------
/** Assemble a local force vector to the global system vector.
 *  Given a local force vector, it needs to be assembled according to global
 *  degree-of-freedom IDs. Possibly, some of these degrees-of-freedom are
 *  constrained and represented by a linear combination of other degrees of
 *  freedom. In that case the vector entries are multiplied by the weights of
 *  this linear combination and assembled to the contributing degrees of
 *  freedom of this constraint.
 *
 *  \tparam SOLVER  Type of system solver
 *  \param[in] forceVec       Local force vector to be assembled
 *  \param[in] doFStatus      Status array for the entries
 *  \param[in] doFIDs         Degree-of-freedom indices for the entries
 *  \param[in] constraints    Linear constraints on the entries
 *  \param[in] solver         Access to the system solver
 */
template<typename SOLVER>
void base::asmb::assembleForces( const base::VectorD&            forceVec,
                                 const std::vector<base::dof::DoFStatus>& doFStatus,
                                 const std::vector<std::size_t>& doFIDs,
                                 const std::vector<
                                     std::pair<unsigned,
                                               std::vector<std::pair<base::number,
                                                                     std::size_t> > > > &
                                 constraints,
                                 SOLVER& solver )
    
{
    // number of active (i.e. non-constrained) dofs 
    const std::size_t numActiveDoFs =
        std::count_if( doFStatus.begin(), doFStatus.end(),
                       boost::bind( std::equal_to<base::dof::DoFStatus>(), _1,
                                    base::dof::ACTIVE ) );

    
    // Collect all dof IDs: active dof IDs + constrained contributing dof IDs
    std::vector<std::size_t> effectiveDoFIDs;
    {
        // active dofs
        for ( unsigned d = 0; d < doFIDs.size(); d++ )
            if ( doFStatus[d] == base::dof::ACTIVE )
                effectiveDoFIDs.push_back( doFIDs[d] );

        // contributing dofs
        for ( unsigned d = 0; d < constraints.size(); d++ ) {
            for ( unsigned d2 = 0; d2 < constraints[d].second.size(); d2++ ) {
                effectiveDoFIDs.push_back( constraints[d].second[d2].second );
            }
        }
    }

    // Result container
    base::VectorD sysVector( static_cast<int>( effectiveDoFIDs.size() ) );

    // Counter of active dofs, constrained dofs and resulting contributing dofs
    unsigned activeCtr = 0;
    unsigned cstrCtr   = 0;
    unsigned extraCtr  = 0;

    // go through all entries of the input force vector
    for ( std::size_t d = 0; d < doFIDs.size(); d ++ ) {

        // check activity 
        if ( doFStatus[d] == base::dof::ACTIVE ) {

            sysVector[    activeCtr ] = forceVec[ d ];
            activeCtr++;
            
        }
        else if ( doFStatus[d] == base::dof::CONSTRAINED ) {

            // local dof ID of constrained dof (sanity check)
            const unsigned localDoFID = constraints[ cstrCtr ].first;
            assert( localDoFID == d );

            // array of linear constraints contributing to this dof
            const std::vector<std::pair<base::number,std::size_t> > constraint =
                constraints[ cstrCtr ].second;

            // go through all contributors
            for ( std::size_t d2 = 0; d2 < constraint.size(); d2++) {

                // weight as mulitplier
                const base::number weight = constraint[d2].first;

                sysVector[ numActiveDoFs + extraCtr ] = weight * forceVec[ d ];
                
                extraCtr++;
            }
            cstrCtr++;
        }
    }

    // Pass on to system solver
    solver.insertToRHS( sysVector, effectiveDoFIDs );

    return;
}

#endif
