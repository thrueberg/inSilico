//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Eigen3.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_solver_eigen3_hpp
#define base_solver_eigen3_hpp

//------------------------------------------------------------------------------
// std   includes
#include <ostream>

// Eigen includes
#include <Eigen/Sparse>
#include <Eigen/Core>

#ifdef LOAD_SUPERLU
#include <Eigen/SuperLUSupport>
#endif

#ifdef LOAD_PARDISO
#define EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif

#ifdef LOAD_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

// base includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace solver{

        class Eigen3;

    }
}

//------------------------------------------------------------------------------
class base::solver::Eigen3
{
public:
    typedef Eigen::Triplet<number> Triplet;

    Eigen3( const std::size_t size )
    {
        rhs_.resize( size );
        rhs_.fill( 0. );
        lhs_.resize( static_cast<int>(size),
                     static_cast<int>(size) );
    }

    //--------------------------------------------------------------------------
    template<typename MATRIX, typename RDOFS, typename CDOFS>
    void insertToLHS( const MATRIX & matrix, 
                      const RDOFS  & rowDofs,
                      const CDOFS  & colDofs )
    {
        const std::size_t numRowDofs = rowDofs.size();
        const std::size_t numColDofs = colDofs.size();

        for ( std::size_t i = 0; i < numRowDofs; i ++ ) {
            for ( std::size_t j = 0; j < numColDofs; j ++ ) {
                triplets_.push_back( Triplet( static_cast<unsigned>(rowDofs[i]),
                                              static_cast<unsigned>(colDofs[j]),
                                              matrix(i,j) ) );
            }
        }

    }


    //--------------------------------------------------------------------------
    template<typename VECTOR, typename DOFS>
    void insertToRHS( const VECTOR & vector,
                      const DOFS   & dofs )
    {
        const std::size_t numDofs = dofs.size();

        assert( static_cast<std::ptrdiff_t>(numDofs) == vector.size() );

        for ( std::size_t i = 0; i < numDofs; i ++ ) {

            rhs_[ dofs[i] ] += vector[i];
        }

    }

    //--------------------------------------------------------------------------
    double norm() const
    {
        return ( rhs_.norm() / static_cast<double>( rhs_.size() ) );
    }

    //--------------------------------------------------------------------------
    void finishAssembly()
    {
        lhs_.setFromTriplets( triplets_.begin(), triplets_.end() );
        
        // make sure that the triplet dies
        triplets_.clear();
        std::vector< Triplet >().swap( triplets_ );
    }

    //--------------------------------------------------------------------------
    void choleskySolve()
    {
        // create a LLT-decomposition from the system matrix
        Eigen::SimplicialLDLT< Eigen::SparseMatrix<number> > chol( lhs_ );

        /* This is the version, which should always work, but does not give
         * access to the factorisation, as needed in the block solver:
         *
         * const VectorD x = chol.solve( rhs_ );
         * rhs_ = x;
         *
         */

        /*  Note that in this version, the individual steps in the
         *  solve process are made explicit.  It appears that the
         *  latest version of Eigen has an inverse notion of P and
         *  Pinv with respect to previous versions. Therefore when
         *  using old versions, in the following the first and forth
         *  lines shall be swapped in position (as indicated by the
         *  post-fixes).
         */


        /*  Solution of
         *        A x = b
         *  deferred to
         *        A' y = AP A P^{-1} y = P b = b'
         *  Steps are based on the factorisation
         *        A' = L D L^T = L D U
         *
         *
         *  1)    b' = P b        (permutate right hand side)
         *  2) L y'' = b'         (solve for y'')
         *  3)   y'  = D^{-1} y'' (diagonal part)
         *  4) U y   = y'         (solve for y)
         *  5)   x   = P^{-1} y   (undo permutation)
         */
        
        rhs_  = chol.permutationP() * rhs_;                   // 1)
        chol.matrixL().solveInPlace( rhs_ );                  // 2)
        rhs_ = chol.vectorD().asDiagonal().inverse() * rhs_;  // 3)
        chol.matrixU().solveInPlace( rhs_ );                  // 4)
        rhs_ = chol.permutationPinv() * rhs_;                 // 5)
    }

    //--------------------------------------------------------------------------
#ifdef LOAD_SUPERLU
    void superLUSolve()
    {
        Eigen::SuperLU< Eigen::SparseMatrix<number> >  superLU( lhs_ );
        VectorD x = superLU.solve( rhs_ );
        rhs_ = x;
    }
#endif

    //--------------------------------------------------------------------------
#ifdef LOAD_PARDISO
    void pardisoLUSolve()
    {
        Eigen::PardisoLU<Eigen::SparseMatrix<number> > pardisoLU( lhs_ );
        VectorD x = pardisoLU.solve( rhs_ );
        rhs_ = x;
    }

    void pardisoCholeskySolve()
    {
        Eigen::PardisoLDLT<Eigen::SparseMatrix<number> > pardisoLDLT( lhs_ );
        VectorD x = pardisoLDLT.solve( rhs_ );
        rhs_ = x;
    }

#endif

    //--------------------------------------------------------------------------
#ifdef LOAD_UMFPACK
    void umfPackLUSolve()
    {
        Eigen::UmfPackLU<Eigen::SparseMatrix<number> > umfPackLU( lhs_ );
        VectorD x = umfPackLU.solve( rhs_ );
        rhs_ = x;
    }
#endif
    

    int cgSolve()
    {
        // fill A and b
        Eigen::ConjugateGradient<Eigen::SparseMatrix<number> > cg;
        cg.compute( lhs_ );
        VectorD x = cg.solve( rhs_ );
        // std::cout << "#iterations: " << cg.iterations() << std::endl;
        // std::cout << "estimated error: " << cg.error() << std::endl;
        rhs_ = x;

        return cg.iterations();
    }        


    //--------------------------------------------------------------------------
    number getValue( const std::size_t index ) const
    {
        return rhs_[ index ];
    }

    //--------------------------------------------------------------------------
    void systemInfo( std::ostream& out ) const
    {
        out << lhs_.rows() << " X " << lhs_.cols() << " sparse matrix with "
            << lhs_.nonZeros() << " non-zero entries \n";
    }

    //--------------------------------------------------------------------------
    void debugLHS( std::ostream & out) const
    {
        for ( int k=0; k < lhs_.outerSize(); k++ ) {
            for ( Eigen::SparseMatrix<number>::InnerIterator it(lhs_,k); it; ++it) {

                out << it.row() << " " << it.col() << " " << it.value() << "\n";
            }
        }
    }

    //--------------------------------------------------------------------------
    void debugRHS( std::ostream & out) const
    {
        for ( int i = 0; i < rhs_.size(); i++ )
            out << i << " " << rhs_[i] << "\n";
    }

    //--------------------------------------------------------------------------
    void debugTriplet( std::ostream& out ) const
    {
        for ( std::size_t t = 0; t < triplets_.size(); t++ ) {
            out << triplets_[t].row()   << "  "
                << triplets_[t].col()   << "  "
                << triplets_[t].value() << "\n";
        }
    }
    
private:
    std::vector< Triplet > triplets_;
    VectorD                     rhs_;
    Eigen::SparseMatrix<number> lhs_;
};

#endif
