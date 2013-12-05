//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Eigen3Block2x2.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_solver_eigen3block2x2_hpp
#define base_solver_eigen3block2x2_hpp

#error FIRST_DRAFT_DOES_NOT_PERFORM

//------------------------------------------------------------------------------
// std   includes
#include <ostream>

// Eigen includes
#include <Eigen/Sparse>
#include <Eigen/Core>

// base includes
#include <base/linearAlgebra.hpp>
#include <base/verify.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace solver{

        class Eigen3Block2x2;

    }
}

//------------------------------------------------------------------------------
/** Solver for 2x2 block system based on the Eigen sparse solver interface.
 *  Consider the system
 *  \f[ 
 *       \left( \matrix{ A & B \cr B^T & C } \right)
 *       \left( \matrix{ x_1 \cr x_2 } \right) =
 *       \left( \matrix{ f_1 \cr f_2 } \right)
 *  \f]
 *  This system is solved by first constructing the Schur complement
 *  \f[
 *       S = C - B^T A^{-1} B
 *  \f]
 *  and solving the reduced system
 *  \f[
 *       S x_2 = g = f_2 - B^T A^{-1} f_1
 *  \f]
 *  and then computing the other unknown field via
 *  \f[
 *       x_2 = A^{-1} (f_1 - B x_2)
 *  \f]
 *
 */
class base::solver::Eigen3Block2x2
{
public:
    typedef Eigen::Triplet<number> Triplet;

    Eigen3Block2x2( const std::size_t size1, const std::size_t size2 )
        : size1_( size1 ), size2_( size2 )
    {
        rhs1_.resize( size1_ ); rhs1_.fill( 0. );
        rhs2_.resize( size2_ ); rhs2_.fill( 0. );
        
        A_.resize( size1_, size1_ );
        B_.resize( size1_, size2_ );
        C_.resize( size2_, size2_ );
    }

    //--------------------------------------------------------------------------
    template<typename MATRIX, typename RDOFS, typename CDOFS>
    void insertToLHS( const MATRIX & matrix, 
                      const RDOFS  & rowDofs,
                      const CDOFS  & colDofs )
    {
        const unsigned numRowDofs = rowDofs.size();
        const unsigned numColDofs = colDofs.size();

        const bool top  = ( rowDofs[0] < size1_ );
        const bool left = ( colDofs[0] < size1_ );

        for ( unsigned i = 0; i < numRowDofs; i ++ ) {
            for ( unsigned j = 0; j < numColDofs; j ++ ) {

                if ( top and left ) 
                    tripletsA_.push_back( Triplet( rowDofs[i],
                                                   colDofs[j],
                                                   matrix(i,j) ) );
                else if ( top and not left )
                    tripletsB_.push_back( Triplet( rowDofs[i],
                                                   colDofs[j] - size1_,
                                                   matrix(i,j) ) );
                else if ( not top and not left )
                    tripletsC_.push_back( Triplet( rowDofs[i] - size1_,
                                                   colDofs[j] - size1_,
                                                   matrix(i,j) ) );

                else VERIFY_MSG( false, "No other block available" );
                
            }
        }

    }


    //--------------------------------------------------------------------------
    template<typename VECTOR, typename DOFS>
    void insertToRHS( const VECTOR & vector,
                      const DOFS   & dofs )
    {
        const unsigned numDofs = dofs.size();

        const bool top = ( dofs[0] < size1_ );

        for ( unsigned i = 0; i < numDofs; i ++ ) {

            if ( top ) rhs1_[ dofs[i]          ] += vector[i];
            else       rhs2_[ dofs[i] - size1_ ] += vector[i];
        }

    }

    //--------------------------------------------------------------------------
    double norm1() const
    {
        return ( rhs1_.norm() / static_cast<double>( rhs1_.size() ) );
    }

    //--------------------------------------------------------------------------
    double norm2() const
    {
        return ( rhs2_.norm() / static_cast<double>( rhs2_.size() ) );
    }

    //--------------------------------------------------------------------------
    void finishAssembly()
    {
        A_.setFromTriplets( tripletsA_.begin(), tripletsA_.end() );
        B_.setFromTriplets( tripletsB_.begin(), tripletsB_.end() );
        C_.setFromTriplets( tripletsC_.begin(), tripletsC_.end() );
        // make sure that the triplet dies
        tripletsA_.clear(); std::vector<Triplet>().swap( tripletsA_ );
        tripletsB_.clear(); std::vector<Triplet>().swap( tripletsB_ );
        tripletsC_.clear(); std::vector<Triplet>().swap( tripletsC_ );
    }

    //--------------------------------------------------------------------------
    void choleskyBlockSolve()
    {
        // LL^T composition of A
        Eigen::SimplicialLLT< Eigen::SparseMatrix<number> > cholA( A_ );
        // A^{-1} * B
        Eigen::SparseMatrix<number> aInvB( B_.rows(), B_.cols() );
        cholA._solve_sparse( B_, aInvB );
        // Schur complement S = C - B^T A^{-1} B
        const Eigen::SparseMatrix<number> S =
            C_ - Eigen::SparseMatrix<number>( B_.transpose() ) * aInvB;
        // combined RHS for Schur complement eqn
        rhs2_ -= Eigen::SparseMatrix<number>( aInvB.transpose() ) * rhs1_;
        // Solve for second unknown
        Eigen::SimplicialLLT< Eigen::SparseMatrix<number> > cholS( S );
        const VectorD x2 = cholS.solve( rhs2_ );
        // Solver auxiliary eqn
        const VectorD y = cholA.solve( rhs1_ );
        // First unknown
        const VectorD x1 = y - (aInvB * x2);
        // pass back to member variable storage in RHS
        rhs1_ = x1;
        rhs2_ = x2;
    }


    //--------------------------------------------------------------------------
    number getValue( const std::size_t index ) const
    {
        if ( index < size1_ ) return rhs1_[ index ];
        return rhs2_[ index - size1_ ];
    }

    //--------------------------------------------------------------------------
    template<typename MATRIX>
    void debug( std::ostream & out, MATRIX& matrix) const
    {
        for ( int k=0; k < matrix.outerSize(); k++ ) {
            for ( Eigen::SparseMatrix<number>::InnerIterator it(matrix,k); it; ++it) {

                out << it.row() << " " << it.col() << " " << it.value() << "\n";
            }
        }
    }

    
private:
    const std::size_t size1_, size2_;
    
    std::vector<Triplet> tripletsA_;
    std::vector<Triplet> tripletsB_;
    std::vector<Triplet> tripletsC_;

    VectorD                     rhs1_, rhs2_;
    Eigen::SparseMatrix<number> A_, B_, C_;
};

#endif
