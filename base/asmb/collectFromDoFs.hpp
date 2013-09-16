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

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        //------------------------------------------------------------------
        //! Helper to collect relevant data from the dof objects
        template<typename FELEMENT>
        void collectFromDoFs( const FELEMENT*           fElement, 
                              std::vector<base::dof::DoFStatus>& status,
                              std::vector<std::size_t>& ids,
                              std::vector<number>&      values,
                              std::vector< std::pair<unsigned,
                                                     std::vector<std::pair<number,
                                                                           std::size_t> >
                                                     > >& constraints );
        
    }
}


//------------------------------------------------------------------------------
/** With access to a field element, collect all the assembly-relevant DoF-data.
 *
 */
template<typename FELEMENT>
void base::asmb::collectFromDoFs( const FELEMENT*           fElement, 
                                  std::vector<base::dof::DoFStatus>& status,
                                  std::vector<std::size_t>& ids,
                                  std::vector<number>&      values,
                                  std::vector< std::pair<unsigned,
                                                         std::vector<std::pair<number,
                                                                               std::size_t> >
                                                         > >& constraints )
{
    typedef typename FELEMENT::DegreeOfFreedom DegreeOfFreedom;
                
    // Get the element's DoF objects
    std::vector<DegreeOfFreedom*> doFs;
    std::copy( fElement -> doFsBegin(),
               fElement -> doFsEnd(),
               std::back_inserter( doFs ) );

    // fill arrays of activity, ID values and possible prescribed values
    for ( unsigned d = 0; d < doFs.size(); d ++ ) {
        doFs[d] -> getStatus(           std::back_inserter( status )   ); 
        doFs[d] -> getIndices(          std::back_inserter( ids )      );
        doFs[d] -> getPrescribedValues( std::back_inserter( values )   );
    }

    // get linear constraints
    typedef typename DegreeOfFreedom::Constraint Constraint;
    unsigned localDoF = 0; // number of the local dof which might be constrained

    for ( unsigned d = 0; d < doFs.size(); d++ ) {

        for ( unsigned s = 0; s < DegreeOfFreedom::size; s++ ){

            // check if DoF is possibly constrained
            if ( status[localDoF] == base::dof::CONSTRAINED ) { 

                // get pairs of weight and global DoF ID
                std::vector<std::pair<number,std::size_t> > weightedDoFIDs;
                doFs[d] -> getConstraint( s ) -> getWeightedDoFIDs( weightedDoFIDs );

                // store pair of local ID and weighted dof IDs
                constraints.push_back( std::make_pair( localDoF, weightedDoFIDs ) );
            }
            
            localDoF++; // increment local DoF counter
        }
    }
    return;
}

#endif
