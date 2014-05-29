//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   collectFromDoFs.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_asmb_collectfromdofs_hpp
#define base_asmb_collectfromdofs_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
// base/dof includes
#include <base/dof/DegreeOfFreedom.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        //------------------------------------------------------------------
        //! Helper to collect relevant data from the dof objects
        template<typename FELEMENT>
        bool collectFromDoFs( const FELEMENT*           fElement, 
                              std::vector<base::dof::DoFStatus>& status,
                              std::vector<std::size_t>& ids,
                              std::vector<number>&      values,
                              std::vector<
                              std::pair<unsigned, std::vector<
                              std::pair<number, std::size_t> >
                              > >& constraints,
                              const bool incremental );
        
    }
}


//------------------------------------------------------------------------------
/** With access to a field element, collect all the assembly-relevant DoF-data.
 *  Going through all DoF pointers of the given element, collect for each vector
 *  component of every DoF
 *
 *  1) The status (ACTIVE, CONSTRAINED, INACTIVE)
 *  2) The global ID (will refer to system matrix row/column)
 *  3) Prescribed values (only relevant for CONSTRAINED), or the difference
 *     between prescribed and current value in case of an incremental analysis
 *
 *  After these data have been retrieved, the linear constraints of every
 *  DoF-component with status CONSTRAINED are registered. If a DoF component
 *  is constrained, from its constraint the slave DoFs are queried. This means
 *  the IDs and coefficients of the linear combination
 *  \f[
 *       u_i = \sum_{j \in C(i)} c_{ij} u_j + g_i
 *  \f]
 *  Here, \f$ C(i) \f$ is the set of indices of slave DoFs for the
 *  constrained DoF component \f$ u_i \f$, \f$ c_{ij} \f$ are the
 *  corresponding weights and \f$ g_i \f$ is the already collected
 *  prescribed value. The set of pairs
 *  \f[
 *       \{ j, c_{ij} \}_{j \in C(i)}
 *  \f]
 *  is associated with a local DoF counter and provided to the caller.
 *  For immersed FE simulations it is frequently the case that there is huge
 *  number of inactive DoFs. Therefore this function returns a flag which
 *  indicates if any DoF is ACTIVE or CONSTRAINED. If this is not the case,
 *  the caller can quit doing whatever followed for this element. 
 *
 *  \tparam FELEMENT  Type of Field Element
 *  \param[in]  fElement Pointer to field element
 *  \param[out] status       Status of all DoF components
 *  \param[out] ids          IDs of all DoF components
 *  \param[out] values       Prescribed values (or difference) for constraints
 *  \param[out] constraints  Weights and IDs of slave DoFs for constrained DoF
 *  \param[in]  incremental  True for an incremental analysis
 *  \return     Flag that indicates if any DoF is ACTIVE or CONSTRAINED
 */
template<typename FELEMENT>
bool base::asmb::collectFromDoFs( const FELEMENT*           fElement, 
                                  std::vector<base::dof::DoFStatus>& status,
                                  std::vector<std::size_t>& ids,
                                  std::vector<number>&      values,
                                  std::vector< std::pair<unsigned,
                                                         std::vector<std::pair<number,
                                                                               std::size_t> >
                                  > >& constraints,
                                  const bool incremental )
{
    typedef typename FELEMENT::DegreeOfFreedom DegreeOfFreedom;
                
    // Get the element's DoF pointers
    std::vector<DegreeOfFreedom*> doFs;
    std::copy( fElement -> doFsBegin(), fElement -> doFsEnd(),
               std::back_inserter( doFs ) );

    // fill arrays of activity, ID values and possible prescribed values
    for ( unsigned d = 0; d < doFs.size(); d ++ ) {
        doFs[d] -> getStatus(           std::back_inserter( status )   ); 
        doFs[d] -> getIndices(          std::back_inserter( ids )      );
        doFs[d] -> getPrescribedValues( std::back_inserter( values ), incremental );
    }

    // get linear constraints
    unsigned localDoF = 0; // number of the local dof which might be constrained

    // avoid computations for inactive dofs
    bool allDoFsAreInactive = true;

    for ( unsigned d = 0; d < doFs.size(); d++ ) {

        for ( unsigned s = 0; s < DegreeOfFreedom::size; s++ ){

            // check if DoF is possibly constrained
            if ( status[localDoF] == base::dof::CONSTRAINED ) { 

                // get pairs of weight and global DoF ID
                std::vector<std::pair<number,std::size_t> > weightedDoFIDs;
                doFs[d] -> getConstraint( s ) -> getWeightedDoFIDs( weightedDoFIDs );

                // store pair of local ID and weighted dof IDs
                constraints.push_back( std::make_pair( localDoF, weightedDoFIDs ) );

                allDoFsAreInactive = false;
            }
            else if ( status[localDoF] == base::dof::ACTIVE ) allDoFsAreInactive = false;
            
            localDoF++; // increment local DoF counter
        }
    }
    
    return (not allDoFsAreInactive);
}

#endif
