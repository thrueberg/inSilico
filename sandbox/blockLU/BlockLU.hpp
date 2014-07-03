//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   BlockLU.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef blocklu_hpp
#define blocklu_hpp

//------------------------------------------------------------------------------
// std   includes
#include <ostream>
#include <vector>
#include <set>
#include <boost/unordered_set.hpp>

// Eigen includes
#include <Eigen/Sparse>
#include <Eigen/Core>


#ifdef LOAD_UMFPACK
#include <umfpack.h>
#else
#error Current implementation requires UMFPACK
#endif

// base includes
#include <base/linearAlgebra.hpp>
#include <base/io/Format.hpp>
#include <base/solver/TripletContainer.hpp>

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/** \ingroup thomas
 *  Block-LU solver, does not perform too well.
 */
class BlockLU
{
public:
    //!
    typedef base::solver::TripletContainer TripletContainer;

    //! Constructor with the size \f$ N \f$ of matrix and vector
    BlockLU( const std::size_t size1, const std::size_t size2 )
        : size1_( size1 ), size2_( size2 )
    {
        f_.resize( size1_ ); f_.fill( 0. );
        g_.resize( size2_ ); g_.fill( 0. );
        
        A_.resize( static_cast<int>( size1_ ),
                   static_cast<int>( size1_ ) );
        B_.resize( static_cast<int>( size1_ ),
                   static_cast<int>( size2_ ) );
        Ct_.resize( static_cast<int>( size1_ ),
                    static_cast<int>( size2_ ) );
        D_.resize( static_cast<int>( size2_ ),
                   static_cast<int>( size2_ ) );

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

        const bool rowsFirst = rowDoFs[0] < size1_;
        const bool colsFirst = colDoFs[0] < size1_;

        const std::size_t rowSize = (rowsFirst ? size1_ : size2_ );
        const std::size_t colSize = (colsFirst ? size1_ : size2_ );

        const std::size_t rowShift = (rowsFirst ? 0 : size1_ );
        const std::size_t colShift = (colsFirst ? 0 : size1_ );

        TripletContainer& tmp = (rowsFirst and colsFirst ? aTmp_ :
                                 (rowsFirst and not colsFirst ? bTmp_ :
                                  (colsFirst ? cTmp_ : dTmp_ ) ) );

        const bool trans = (not rowsFirst and colsFirst);
                                  
        
        for ( std::size_t i = 0; i < numRowDoFs; i ++ ) {
            for ( std::size_t j = 0; j < numColDoFs; j ++ ) {
                const std::size_t rowIndex = rowDoFs[i] - rowShift;
                const std::size_t colIndex = colDoFs[j] - colShift;

                VERIFY_MSG( rowIndex < rowSize,
                            "Row index out of bound: " + x2s( rowIndex ) );
                VERIFY_MSG( colIndex < colSize, 
                            "Col index out of bound: " + x2s( colIndex ) );

                const unsigned rowIndex2 = static_cast<unsigned>(
                    trans? colIndex : rowIndex );
                const unsigned colIndex2 = static_cast<unsigned>(
                    trans? rowIndex : colIndex );

                tmp.insert( rowIndex2, colIndex2, matrix( i, j ) );
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
        const bool first = dofs[0] < size1_;
        const std::size_t size = (first ? size1_ : size2_ );
        const std::size_t shift = (first ? 0 : size1_ );

        base::VectorD& tmp = (first ? f_ : g_ );
                                  
        for ( std::size_t i = 0; i < numNewDoFs; i ++ ) {
            const std::size_t index = dofs[i] - shift;
            VERIFY_MSG( index < size, x2s( index ) + " out of bound" );
            tmp[ index ] += vector[i];
        }
        return;
    }

    //--------------------------------------------------------------------------
    //! Give the norm of the rhs/solution vector
    double norm() const
    {
        const double normF = (size1_ > 0 ? f_.norm() / size1_ : 0.);
        const double normG = (size2_ > 0 ? g_.norm() / size2_ : 0.);
        return
            std::sqrt( normF*normF + normG*normG );
    }


    //--------------------------------------------------------------------------
    //! Convert the triplet to sparse matrix storage
    void finishAssembly( const bool destroyTriplet = true )
    {
        aTmp_.prepare();
        bTmp_.prepare();
        cTmp_.prepare();
        dTmp_.prepare();
        

        // fill sparse matrix from tripletContainer
        A_.setFromTriplets( aTmp_.begin(), aTmp_.end() );
        B_.setFromTriplets( bTmp_.begin(), bTmp_.end() );
        Ct_.setFromTriplets( cTmp_.begin(), cTmp_.end() );
        D_.setFromTriplets( dTmp_.begin(), dTmp_.end() );

        if ( destroyTriplet ) {
            aTmp_.destroy();
            bTmp_.destroy();
            cTmp_.destroy();
            dTmp_.destroy();
        }
        
        return;
    }


    //--------------------------------------------------------------------------
    void umfPackLUSolve()
    {
        void* numericA;
        double control [UMFPACK_CONTROL];
        umfpack_di_defaults(control);
        control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
        
        // Factorise A
        {
            void* symbolicA;
            const int errorCodeA1 =
                umfpack_di_symbolic( size1_, size1_, 
                                     A_.outerIndexPtr(),
                                     A_.innerIndexPtr(),
                                     A_.valuePtr(),
                                     &symbolicA, control, 0 );
            VERIFY_MSG( errorCodeA1==0,
                        x2s( "Symbolic factorisation of A failed" ) +
                        x2s( errorCodeA1 ) );

            const int errorCodeA2 =
                umfpack_di_numeric( A_.outerIndexPtr(),
                                    A_.innerIndexPtr(),
                                    A_.valuePtr(),
                                    symbolicA,
                                    &numericA, control, 0 );
            VERIFY_MSG( errorCodeA2==0,
                        x2s( "Numeric factorisation of A failed: " ) +
                        x2s( errorCodeA2 ) );
            
            umfpack_di_free_symbolic( &symbolicA );
        }

        // Solve for UB and LCt
        Eigen::SparseMatrix<number> UB, LCt;
        UB.resize(  size1_, size2_ );
        LCt.resize( size1_, size2_ );
        
        for ( std::size_t c = 0; c < size2_; c++ ) {

            base::VectorD colB = B_.col(c);
            base::VectorD colUB; colUB.resize( size1_ );

            const int errorCodeB1 =
                umfpack_di_scale( &colB[0],
                                  &colB[0],
                                  numericA );
            
            const int errorCodeB2 =
                umfpack_di_solve( UMFPACK_Pt_L,
                                  A_.outerIndexPtr(),
                                  A_.innerIndexPtr(),
                                  A_.valuePtr(),
                                  &colUB[0],
                                  &colB[0],
                                  numericA,
                                  control, 0 );

            // insert to LB
            for ( std::size_t r = 0; r < size1_; r++ ) {
                if ( std::abs( colUB[r] > 1.e-10 ) ) UB.insert( r, c ) = colUB[r];
            }

            base::VectorD colCt = Ct_.col(c);
            base::VectorD colLCt; colLCt.resize( size1_ );
            
            const int errorCodeC2 =
                umfpack_di_solve( UMFPACK_Q_Ut,
                                  A_.outerIndexPtr(),
                                  A_.innerIndexPtr(),
                                  A_.valuePtr(),
                                  &colLCt[0],
                                  &colCt[0],
                                  numericA,
                                  control, 0 );

            // insert to LCt;
            for ( std::size_t r = 0; r < size1_; r++ ) {
                if ( std::abs( colLCt[r] > 1.e-10 ) ) LCt.insert( r, c ) = colLCt[r];
            }

        }

        // update Schur complement
        D_ -= LCt.transpose() * UB;

        std::cout << "dim(S) = " << size2_ << " nnz(S) = " << D_.nonZeros()
                  << std::endl;
        
        // Factorise S = D
        void* numericS;
        {
            void* symbolicS;
            const int errorCodeS1 =
                umfpack_di_symbolic( size2_, size2_, 
                                     D_.outerIndexPtr(),
                                     D_.innerIndexPtr(),
                                     D_.valuePtr(),
                                     &symbolicS, 0, 0 );
            VERIFY_MSG( errorCodeS1==0,
                        x2s( "Symbolic factorisation of S failed " ) +
                        x2s( errorCodeS1 ) );

            const int errorCodeS2 =
                umfpack_di_numeric( D_.outerIndexPtr(),
                                    D_.innerIndexPtr(),
                                    D_.valuePtr(),
                                    symbolicS,
                                    &numericS, 0, 0 );
            VERIFY_MSG( errorCodeS2==0,
                        x2s( "Numeric factorisation of S failed: " ) +
                        x2s( errorCodeS2 ) );
            
            umfpack_di_free_symbolic( &symbolicS );
        }

        // update rhs terms
        {
            base::VectorD f2( size1_ );

            const int errorCodeF1 =
                umfpack_di_scale( &f_[0],
                                  &f_[0],
                                  numericA );

            const int errorCodeF2 = 
                umfpack_di_solve( UMFPACK_Pt_L,
                                  A_.outerIndexPtr(),
                                  A_.innerIndexPtr(),
                                  A_.valuePtr(),
                                  &f2[0],
                                  &f_[0],
                                  numericA,
                                  control, 0 );
            f_ = f2;
                                  
        }
        g_ -= LCt.transpose() * f_;

        // Solve Schur complement system
        {
            base::VectorD y; y.resize( size2_ );
            const int errorCodeS =
                umfpack_di_solve( UMFPACK_A,
                                  D_.outerIndexPtr(),
                                  D_.innerIndexPtr(),
                                  D_.valuePtr(),
                                  &y[0],
                                  &g_[0],
                                  numericS, 0, 0 );
            VERIFY_MSG( errorCodeS==0,
                        x2s("Solve Sy=g failed: ") + x2s(errorCodeS) );
            g_ = y;

            umfpack_di_free_numeric( &numericS );
        }

        // update RHS
        f_ -= UB * g_;

        // Solve rest
        base::VectorD x; x.resize( size1_ );
        const int errorCodeA3 =
            umfpack_di_solve( UMFPACK_U_Qt,
                              A_.outerIndexPtr(),
                              A_.innerIndexPtr(),
                              A_.valuePtr(),
                              &x[0],
                              &f_[0],
                              numericA, control, 0 );
        VERIFY_MSG( errorCodeA3==0,
                    x2s("Solve Ax=b failed: ") + x2s(errorCodeA3) );
        f_ = x;

        umfpack_di_free_numeric(  &numericA );
    }



    //@}
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    //! Direct access to an entry in the RHS/solution vector
    number getValue( const std::size_t index ) const
    {
        if ( index < size1_ ) return f_[ index ];
        return g_[ index-size1_ ];
    }

    //--------------------------------------------------------------------------
    //! @name Debug routines for printing
    //@{

    //void debugLHS( std::ostream & out) const
    //{
    //    for ( int k=0; k < A_.outerSize(); k++ ) {
    //        for ( Eigen::SparseMatrix<number>::InnerIterator it(A_,k); it; ++it) {
    //            out << it.row() << " " << it.col() << " " << it.value() << "\n";
    //        }
    //    }
    //}

    void debugRHS( std::ostream & out) const
    {
        for ( int i = 0; i < f_.size(); i++ ) out << i << " " << f_[i] << "\n";
        for ( int i = 0; i < g_.size(); i++ ) out << i << " " << g_[i] << "\n";
    }

    //@}

    
    //--------------------------------------------------------------------------
    //! Delegate registering of test and trial field DoFs to tripletContainer
    //template<typename FIELDTUPLEBINDER, typename FIELDBINDER>
    //void registerFields( const FIELDBINDER& fieldBinder )
    //{
    //    tripletContainer_.registerFields<FIELDTUPLEBINDER>( fieldBinder );
    //}
    
private:
    const std::size_t size1_;
    const std::size_t size2_;
    
    TripletContainer       aTmp_, bTmp_, cTmp_, dTmp_;
    base::VectorD               f_, g_;            //!< Given force vector
    Eigen::SparseMatrix<number> A_, B_, Ct_, D_;    //!< Matrix blocks
};

#endif
