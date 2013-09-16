//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   assembleMatrix.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_asmb_assemblematrix_hpp
#define base_asmb_assemblematrix_hpp

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
// base/dof includes
#include <base/dof/DegreeOfFreedom.hpp> // for the DoFStatus enum

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        template<typename SOLVER>
        void assembleMatrix( const base::MatrixD&             elemMatrix,
                             const std::vector<base::dof::DoFStatus>& rowDoFStatus,
                             const std::vector<base::dof::DoFStatus>& colDoFStatus,
                             const std::vector<std::size_t>&  rowDoFIDs,
                             const std::vector<std::size_t>&  colDoFIDs,
                             const std::vector<base::number>& colDoFValues,
                             const std::vector<
                                 std::pair<unsigned,
                                           std::vector<std::pair<base::number,
                                                                 std::size_t> > > > &
                             rowConstraints,
                             const std::vector<
                                 std::pair<unsigned,
                                           std::vector<std::pair<base::number,
                                                                 std::size_t> > > > &
                             colConstraints,
                             SOLVER& solver,
                             const bool isBubnov, 
                             const bool zeroConstraints );


        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            /** Helper function to assemble a row to the matrix and vector.
             */
            void assembleRow( const std::size_t                numCols,
                              const std::size_t                numActiveColDoFs, 
                              const std::vector<base::dof::DoFStatus>& colDoFStatus,
                              const std::vector<base::number>& colDoFValues,
                              const std::vector<
                                  std::pair<unsigned,
                                            std::vector<std::pair<base::number,
                                                                  std::size_t> > > > &
                              colConstraints,
                              const unsigned r,
                              const unsigned rowCtr,
                              const base::number   rowWeight,
                              const base::MatrixD& elemMatrix,
                              base::MatrixD&       sysMatrix,
                              base::VectorD&       sysVector )
            {
                // counters of active column dofs and constraint target dofs
                unsigned activeColCtr = 0;
                unsigned colCstrCtr   = 0;
                unsigned extraColCtr  = 0;

                // go through all columns of the element matrix
                for ( std::size_t c = 0; c < numCols; c++ ) {

                    // check if the column dof is active
                    if ( colDoFStatus[c] == base::dof::ACTIVE ) {

                        // insert to matrix
                        sysMatrix( rowCtr, activeColCtr ) =
                            rowWeight * elemMatrix( r, c );
                        
                        // increment corresponding column counter
                        activeColCtr++;
                        
                    }
                    else if ( colDoFStatus[c] == base::dof::CONSTRAINED ) {

                        // pass matrix entry, weighted by constraint value, to RHS 
                        sysVector[ rowCtr ] -=
                            colDoFValues[ c ] * rowWeight * elemMatrix( r, c );

                        // local ID of the constrained DoF for sanity check only
                        const unsigned localDofID = colConstraints[colCstrCtr].first;

                        // assertion of correct ID
                        assert( localDofID == c );

                        // the linear constraint of this DoF
                        std::vector<std::pair<base::number,std::size_t> > constraints =
                            colConstraints[ colCstrCtr ].second;

                        // go through all linear constraints 
                        for ( std::size_t c2 = 0; c2 < constraints.size(); c2++ ) {


                            // multiplier of master DoF entry
                            const base::number colWeight = constraints[c2].first;

                            // insert to additional columns
                            sysMatrix( rowCtr,
                                       numActiveColDoFs + extraColCtr ) =
                                rowWeight * colWeight * elemMatrix( r, c );

                            // increment counter of additional columns
                            extraColCtr++;
                        }

                        colCstrCtr++;
                        
                    } // end of if-else

                }
                
                return;
            }

        } // namespace detail_
        //----------------------------------------------------------------------
        
    }
}

//------------------------------------------------------------------------------
/** Assemble an element matrix to the system matrix.
 *  Given an element matrix, the task is to insert its entries into the system's
 *  global matrix and global vector according to the degree-of-freedom numbers
 *  and possible linear constraints on these.
 *
 *  In order to illustrate the mechanism, consider an example. It is given the
 *  element matrix \f$ K \f$ of dimension \f$ 3 \times 4 \f$
 *  \f[
 *       K = \left( \matrix{ K_{11} & K_{12} & K_{13} & K_{14} \cr
 *                           K_{21} & K_{22} & K_{23} & K_{24} \cr
 *                           K_{31} & K_{32} & K_{33} & K_{34} } \right)
 *  \f]
 *  The 3rd row degree of freedom and the 3rd column degree of freedom are
 *  constrained, such that we have the activity arrays
 *  \f[
 *      a_r = \left( \matrix{ 1 & 1 & 0 } \right) \quad
 *      a_c = \left( \matrix{ 1 & 1 & 0 & 1 } \right)
 *  \f]
 *  and the ID arrays
 *  \f[
 *      d_r = \left( \matrix{ r_1 & r_2 & X } \right) \quad
 *      d_c = \left( \matrix{ c_1 & c_2 & X & c_4 } \right)
 *  \f]
 *  with \f$ X \f$ indicating an invalid entry which will not be used.
 *  Let \f$ u_i \f$ be the \f$ i \f$-th column degree-of-freedom and
 *  \f$ v_j \f$ the \f$ j \f$-th row degree-of-freedom. We shall have the
 *  constraints
 *  \f[
 *         u_3 = g \quad v_3 = \alpha v_4 + \beta v_5
 *  \f]
 *  Note that \f$ g \f$ is a given number and would correspond to the nodal
 *  value of the degree of freedom \f$ u_3 \f$ due to e.g. a Dirichlet boundary
 *  condition. On the other hand, the FE space for the rows is constrained such
 *  that entries corresponding to \f$ v_3 \f$ are represented by some linear
 *  combination of  \f$ v_4 \f$ and \f$ v_5 \f$. The latter contributions have
 *  degree of freedom numbers \f$ r_4 \f$ and \f$ r_5 \f$ which can be but need
 *  not be from the \f$ r_j \f$ that are already used for this element.
 *  
 *  Since degree of freedom indices are not checked for redundancies, we will
 *  have in this example 4-1=3 active columns and 3-1+2=4 rows. Moreover, there
 *  are 4 entries to the system's right hand side. For this example we get the
 *  effective dof index arrays
 *  \f[
 *      \tilde{d}_r = \left( \matrix{ r_1 & r_2 & r_4 & r_5 } \right) \quad
 *      \tilde{d}_c = \left( \matrix{ c_1 & c_2 & c_4 } \right)
 *  \f]
 *  The \f$ 4 \times 3 \f$ matrix to enter the system matrix is
 *  \f[
 *     \tilde{K} = \left( \matrix{
 *                   K_{11} &        K_{12} &        K_{14} \cr
 *                   K_{21} &        K_{22} &        K_{24} \cr
 *            \alpha K_{31} & \alpha K_{32} & \alpha K_{34} \cr
 *            \beta  K_{31} & \beta  K_{32} & \beta  K_{34} \cr } \right)
 *  \f]
 *  The contribution to the system's right hand side is
 *  \f[
 *     \tilde{V} = -g \left( \matrix{
 *                     K_{13} & K_{23} & \alpha K_{33} & \beta K_{43}
 *                                  } \right)^\top
 *  \f]
 *
 *  \tparam SOLVER Type of system solver to receive the entries
 *  \param[in] elemMatrix     Element stiffness matrix
 *  \param[in] rowDoFStatus   Status array for the rows
 *  \param[in] colDoFStatus   Status array for the columns
 *  \param[in] rowDoFIDs      Degree-of-freedom indices for the rows
 *  \param[in] colDoFIDs      Degree-of-freedom indices for the columns
 *  \param[in] colDoFValues   Prescribed values of the column degrees-of-freedom
 *  \param[in] rowConstraints Linear constraints on the rows
 *  \param[in] colConstraints Linear constraints on the columns
 *  \param[in] solver         Access to the system solver
 *  \param[in] isBubnov       True if column and row spaces are identical
 *  \param[in] zeroConstraints True if the prescribed values are set to zero
 *                      as in the sub-sequent iterations of a nonlinear solve
 */
template<typename SOLVER>
void base::asmb::assembleMatrix( const base::MatrixD&             elemMatrix,
                                 const std::vector<base::dof::DoFStatus>& rowDoFStatus,
                                 const std::vector<base::dof::DoFStatus>& colDoFStatus,
                                 const std::vector<std::size_t>&  rowDoFIDs,
                                 const std::vector<std::size_t>&  colDoFIDs,
                                 const std::vector<base::number>& colDoFValues,
                                 const std::vector<
                                     std::pair<unsigned,
                                               std::vector<std::pair<base::number,
                                                                     std::size_t> > > > &
                                 rowConstraints,
                                 const std::vector<
                                     std::pair<unsigned,
                                               std::vector<std::pair<base::number,
                                                                     std::size_t> > > > &
                                 colConstraints,
                                 SOLVER& solver,
                                 const bool isBubnov, 
                                 const bool zeroConstraints )
{
    // number of active dofs in the row space
    const std::size_t numActiveRowDoFs =
        std::count_if( rowDoFStatus.begin(), rowDoFStatus.end(),
                       boost::bind( std::equal_to<base::dof::DoFStatus>(), _1,
                                    base::dof::ACTIVE ) );
    
    // number of active dofs in the column space
    const std::size_t numActiveColDoFs =
        ( isBubnov ? numActiveRowDoFs :
          std::count_if( colDoFStatus.begin(), colDoFStatus.end(),
                         boost::bind( std::equal_to<base::dof::DoFStatus>(), _1,
                                      base::dof::ACTIVE ) ) );

    // Collect all ID numbers
    std::vector<std::size_t> effRowDoFIDs, effColDoFIDs;
    {
        // collect active row dof IDs from this element
        for ( unsigned r = 0; r < rowDoFIDs.size(); r++ )
            if ( rowDoFStatus[r] == base::dof::ACTIVE )
                effRowDoFIDs.push_back( rowDoFIDs[r] );

        // collect master dof IDs from linear constraints
        for ( unsigned r = 0; r < rowConstraints.size(); r++ ) {
            for ( unsigned r2 = 0; r2 < (rowConstraints[r]).second.size(); r2++ ) {
                effRowDoFIDs.push_back( (rowConstraints[r].second[r2].second) );
            }
        }

        // if test- and trial-spaces are equal, just copy the IDs
        if ( isBubnov ) effColDoFIDs = effRowDoFIDs;
        // otherwise, do the same for the column (i.e trial) space
        else {
            for ( unsigned c = 0; c < colDoFIDs.size(); c++ )
                if ( colDoFStatus[c] == base::dof::ACTIVE )
                    effColDoFIDs.push_back( colDoFIDs[c] );

            for ( unsigned c = 0; c < colConstraints.size(); c++ ) {
                for ( unsigned c2 = 0; c2 < (colConstraints[c]).second.size(); c2++ ) {
                    effColDoFIDs.push_back( (colConstraints[c].second[c2].second) );
                }
            }
        }
    }
    

    // Result container for LHS and RHS contributions
    MatrixD sysMatrix( static_cast<int>( effRowDoFIDs.size() ),
                       static_cast<int>( effColDoFIDs.size() ) );
    VectorD sysVector =
        VectorD::Zero( static_cast<int>( effRowDoFIDs.size() ) );

    // Counters for active dofs, constraints and additional dofs 
    unsigned activeRowCtr = 0;
    unsigned rowCstrCtr   = 0;
    unsigned extraRowCtr  = 0; 

    // go through all rows of the element matrix
    for ( std::size_t r = 0; r < rowDoFIDs.size(); r++ ) {

        // check if row belongs to an active DoF
        if ( rowDoFStatus[r] == base::dof::ACTIVE ) {

            // assemble the entire row
            detail_::assembleRow( colDoFIDs.size(), numActiveColDoFs, colDoFStatus,
                                  colDoFValues, colConstraints,
                                  static_cast<unsigned>(r), activeRowCtr,
                                  1.0, elemMatrix, sysMatrix, sysVector );

            activeRowCtr++;
        }
        else if ( rowDoFStatus[r] == base::dof::CONSTRAINED ) {

            // sanity check of the passed constrained array
            const unsigned localDoFID = rowConstraints[ rowCstrCtr ].first;
            assert( localDoFID == r );

            // get vector of constrains for this row
            const std::vector<std::pair<base::number,std::size_t> > constraints =
                rowConstraints[ rowCstrCtr ].second;

            // go through all rows which contribute to this constrained one
            for ( std::size_t r2 = 0; r2 < constraints.size(); r2++) {

                // weight multiplies the row
                const base::number rowWeight = constraints[r2].first;

                // assemble additional row due to linear constraint
                detail_::assembleRow( colDoFIDs.size(), numActiveColDoFs, colDoFStatus,
                                      colDoFValues, colConstraints,
                                      static_cast<unsigned>(r),
                                      static_cast<unsigned>(numActiveRowDoFs) + extraRowCtr,
                                      rowWeight, elemMatrix, sysMatrix, sysVector );
                extraRowCtr++;
            }

            rowCstrCtr++;
            
        }
    }
    
    // Pass on to system matrix
    solver.insertToLHS( sysMatrix, effRowDoFIDs, effColDoFIDs );
    if ( not zeroConstraints )
        solver.insertToRHS( sysVector, effRowDoFIDs );
            
    return;
}

#endif
