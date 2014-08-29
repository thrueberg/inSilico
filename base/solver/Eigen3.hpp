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
#include <vector>
#include <set>
#include <boost/unordered_set.hpp>

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
#include <base/io/Format.hpp>
#include <base/solver/TripletContainer.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace solver{
        class Eigen3;
    }
}

//------------------------------------------------------------------------------
/** Solver for linear system based on the Eigen3 library.
 *  This object holds the data and provides solution routines for the
 *  linear system of equations
 *  \f[
 *         A x = b
 *  \f]
 *  where \f$ A \f$ is a \f$ N \times N \f$ matrix, preferably invertible,
 *  \f$ x \f$ the vector of unknowns and \f$ b \f$ the given right-hand side,
 *  both of size \f$ N \f$.
 *  The storage and solution functionality of this object is entirely based on
 *  <a href="http://eigen.tuxfamily.org/">Eigen3</a> which provides the
 *  interface to several solver packages.
 *  
 */
class base::solver::Eigen3
{
public:
    //!
    typedef base::solver::TripletContainer TripletContainer;

    //! Constructor with the size \f$ N \f$ of matrix and vector
    Eigen3( const std::size_t size )
    {
        b_.resize( size );
        b_.fill( 0. );
        A_.resize( static_cast<int>(size),
                   static_cast<int>(size) );
    }

    //--------------------------------------------------------------------------
    //! Insert numbers to matrix storage
    template<typename MATRIX, typename RDOFS, typename CDOFS>
    void insertToLHS( const MATRIX & matrix, 
                      const RDOFS  & rowDoFs,
                      const CDOFS  & colDoFs )
    {
        const std::size_t numRowDoFs = rowDoFs.size();
        const std::size_t numColDoFs = colDoFs.size();
        const std::size_t numTotalDoFs = b_.size();
        
        for ( std::size_t i = 0; i < numRowDoFs; i ++ ) {
            for ( std::size_t j = 0; j < numColDoFs; j ++ ) {
                const std::size_t rowIndex = rowDoFs[i];
                const std::size_t colIndex = colDoFs[j];

                VERIFY_MSG( rowIndex < numTotalDoFs,
                            "Row index out of bound: " + x2s( rowIndex ) );
                VERIFY_MSG( colIndex < numTotalDoFs,
                            "Col index out of bound: " + x2s( colIndex ) );

                tripletContainer_.insert( static_cast<unsigned>( rowIndex ),
                                          static_cast<unsigned>( colIndex ),
                                          matrix( i, j ) );
            }
        }
        return;
    }


    //--------------------------------------------------------------------------
    //! Insert numbers to RHS vector
    template<typename VECTOR, typename DOFS>
    void insertToRHS( const VECTOR & vector,
                      const DOFS   & dofs )
    {
        const std::size_t numNewDoFs = dofs.size();
        const std::size_t numTotalDoFs = b_.size();
        
        for ( std::size_t i = 0; i < numNewDoFs; i ++ ) {
            const std::size_t index = dofs[i];
            VERIFY_MSG( index < numTotalDoFs, x2s( index ) + " out of bound" );
            b_[ dofs[i] ] += vector[i];
        }
        return;
    }

    //--------------------------------------------------------------------------
    //! Give the norm of the rhs/solution vector
    double norm() const
    {
        return ( this -> norm( 0, b_.size() ) );
    }

    double norm( const std::size_t first,
                 const std::size_t last  ) const
    {
        return ( b_.segment( first, last-first).norm() /
                 static_cast<double>( b_.segment( first, last-first).size() ) );
    }

    //--------------------------------------------------------------------------
    //! Convert the triplet to sparse matrix storage
    void finishAssembly( const bool destroyTriplet = true )
    {
        tripletContainer_.prepare();
        

        // fill sparse matrix from tripletContainer
        A_.setFromTriplets( tripletContainer_.begin(), tripletContainer_.end() );

        if ( destroyTriplet ) tripletContainer_.destroy();
        
        return;
    }

    //--------------------------------------------------------------------------
    //! For A s.p.d., solution by a Cholesky method
    void choleskySolve()
    {
        // create a LLT-decomposition from the system matrix
        Eigen::SimplicialLDLT< Eigen::SparseMatrix<number> > chol( A_ );

        /* This is the version, which should always work, but does not give
         * access to the factorisation, as needed in the block solver:
         *
         * const VectorD x = chol.solve( b_ );
         * b_ = x;
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
        
        b_  = chol.permutationP() * b_;                   // 1)
        chol.matrixL().solveInPlace( b_ );                // 2)
        b_ = chol.vectorD().asDiagonal().inverse() * b_;  // 3)
        chol.matrixU().solveInPlace( b_ );                // 4)
        b_ = chol.permutationPinv() * b_;                 // 5)
    }

    //--------------------------------------------------------------------------
    //! @name Solve routines by means of external packages
    //@{
#ifdef LOAD_SUPERLU
    void superLUSolve()
    {
        Eigen::SuperLU< Eigen::SparseMatrix<number> >  superLU( A_ );
        VectorD x = superLU.solve( b_ );
        b_ = x;
    }
#endif

    void luSolve()
    {
#ifdef LOAD_PARDISO
        this -> pardisoLUSolve();
#elif defined(LOAD_UMFPACK)
        this -> umfPackLUSolve();
#elif defined(LOAD_SUPERLU)
        this -> superLUSolve();
#else
        Eigen::SparseLU< Eigen::SparseMatrix<number> > luSolver;
        luSolver.analyzePattern( A_ );
        luSolver.factorize( A_ );
        VectorD x = luSolver.solve( b_ );
        b_ = x;
#endif
        return;
    }
    

#ifdef LOAD_PARDISO
    void pardisoLUSolve( const bool outOfCore = true )
    {
        Eigen::PardisoLU<Eigen::SparseMatrix<number> > pardisoLU;
        //Eigen::PardisoLU<Eigen::SparseMatrix<number> > pardisoLU( A_ );
        if ( not outOfCore ) pardisoLU.pardisoParameterArray()[59] = 0;
        pardisoLU.compute( A_ );
        VectorD x = pardisoLU.solve( b_ );
        b_ = x;
    }

    void pardisoCholeskySolve()
    {
        Eigen::PardisoLDLT<Eigen::SparseMatrix<number> > pardisoLDLT( A_ );
        VectorD x = pardisoLDLT.solve( b_ );
        b_ = x;
    }

#endif

#ifdef LOAD_UMFPACK
    void umfPackLUSolve()
    {
        Eigen::UmfPackLU<Eigen::SparseMatrix<number> > umfPackLU( A_ );
        VectorD x = umfPackLU.solve( b_ );
        b_ = x;
    }
#endif
    //@}
    //--------------------------------------------------------------------------
    
    //! Use a conjugate gradient method for solution (A=A' and A > 0 !!)
    int cgSolve()
    {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<number> > cg;
        cg.compute( A_ );
        VectorD x = cg.solve( b_ );
        // std::cout << "#iterations: " << cg.iterations() << std::endl;
        // std::cout << "estimated error: " << cg.error() << std::endl;
        b_ = x;

        return cg.iterations();
    }

    int biCGStabSolve()
    {
        typedef Eigen::IncompleteLUT<number> PreCond;
        
        Eigen::BiCGSTAB<Eigen::SparseMatrix<number>,PreCond > biCG;
        biCG.compute( A_ );
        VectorD x = biCG.solve( b_ );
        // std::cout << "#iterations: "     << biCG.iterations() << std::endl;
        // std::cout << "estimated error: " << biCG.error() << std::endl;
        b_ = x;

        return biCG.iterations();
    }
    
    //--------------------------------------------------------------------------
    //! Direct access to an entry in the RHS/solution vector
    number getValue( const std::size_t index ) const
    {
        return b_[ index ];
    }

    //--------------------------------------------------------------------------
    //! @name Debug routines for printing
    //@{
    void systemInfo( std::ostream& out ) const
    {
        out << A_.rows() << " X " << A_.cols() << " sparse matrix with "
            << A_.nonZeros() << " non-zero entries \n";
    }

    void debugLHS( std::ostream & out) const
    {
        for ( int k=0; k < A_.outerSize(); k++ ) {
            for ( Eigen::SparseMatrix<number>::InnerIterator it(A_,k); it; ++it) {

                out << it.row() << " " << it.col() << " " << it.value() << "\n";
            }
        }
    }

    void debugRHS( std::ostream & out) const
    {
        for ( int i = 0; i < b_.size(); i++ )
            out << i << " " << b_[i] << "\n";
    }

    std::ostream& debugTriplet( std::ostream& out ) const
    {
        return tripletContainer_.write( out );
    }
    //@}

    
    //--------------------------------------------------------------------------
    //! Delegate registering of test and trial field DoFs to tripletContainer
    template<typename FIELDTUPLEBINDER, typename FIELDBINDER>
    void registerFields( const FIELDBINDER& fieldBinder )
    {
        tripletContainer_.registerFields<FIELDTUPLEBINDER>( fieldBinder );
    }
    
private:
    TripletContainer            tripletContainer_; //!< Temp. storage of triplets
    VectorD                     b_;                //!< Given force vector
    Eigen::SparseMatrix<number> A_;                //!< Sparse matrix
};

#endif
