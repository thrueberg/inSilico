//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Eigen.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef mat_eigen_hpp
#define mat_eigen_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/verify.hpp>
#include <base/linearAlgebra.hpp>
// boost includes
#include <boost/array.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>

//------------------------------------------------------------------------------
namespace mat{

    //--------------------------------------------------------------------------
    /** Compute the three eigenvalues of a symmetric tensor.
     *  This method is taken from
     *  http://en.wikipedia.org/wiki/Eigenvalue_algorithm under the subsection
     *  of 3x3 matrices. The basis of the method is that if \f$ A = p B + q I\f$,
     *  the eigenvalues \f$ \lambda \f$ of \f$ A \f$ relate to the eigenvalues
     *  \f$ \mu \f$ by \f$ \lambda = p \mu + q \f$. It is know that
     *  the eigenvalues of a 3x3 matrix are the roots of a cubic equation.
     *  Now, the values of \f$ p \f$ and \f$ q \f$ are chosen such that this
     *  cubic equation for the matrix \f$ B \f$ is easier to solve.
     *  The proposed choice is \f$ p = \sqrt{ tr( (A - qI )^2 ) } / 6 \f$ and
     *  \f$ q = tr(A) / 3 \f$. Then we get
     *  \f[
     *        \mu = 2 \cos \left( \frac{1}{3} acos( (\det B)/2 ) +
     *                            \frac{2 k \pi}{3} \right), k = 0,1,2
     *  \f]
     *  and from these the values of \f$ \lambda \f$. Note that for efficieny
     *  reasons the sum of the off-diagonal squares is compared to the square
     *  of the trace in the beginning. If they are significant smaller then
     *  the matrix is assumed to be diagonal and the computation is trivial.
     */
    inline boost::array<double,3> eigenValues( const Tensor& A )
    {
        boost::array<double,3> eigenVals;
        const double q = trace( A ) / 3.;
        
        const double offDiagSquares =
            (A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2));

        if ( ( offDiagSquares/ (q*q) ) <  1.e-10 ) {
            // Matrix is practically diagonal
            eigenVals[0] = A(0,0);
            eigenVals[1] = A(1,1);
            eigenVals[2] = A(2,2);
        }
        else {
            // Matrix is not diagonal

            // computation of value of p
            const double aux =
                (A(0,0)-q) * (A(0,0)-q) +
                (A(1,1)-q) * (A(1,1)-q) +
                (A(2,2)-q) * (A(2,2)-q) +
                2. * offDiagSquares;
            const double p = std::sqrt( aux / 6. );

            // shifted matrix with easier eigenvalue equation
            const Tensor B = (1./p) * (A - q * Tensor::Identity() );

            // aux variable
            const double r = determinant( B ) / 2.;

            // compute cos^{-1} of r (make sure -1<=r<=1)
            const double phi =
                ( r <= -1.0 ? M_PI/3. :
                  ( r >= 1.0 ? 0. :
                    std::acos( r ) /3. ) );

            // eigenvalues of A
            eigenVals[0] = q + 2. * p * std::cos( phi );
            eigenVals[2] = q + 2. * p * std::cos( phi + 2.*M_PI/3. );
            eigenVals[1] = 3. * q - eigenVals[0] - eigenVals[2];
        }
        
        return eigenVals;
        
    }

    //--------------------------------------------------------------------------
    namespace detail_{

        /** Construct an eigenvector corresponding to a given eigenvalue.
         *  The requirements are (i) the matrix A is symmetric and (ii) the
         *  given eigenvalue has multiplicity one.
         *  Then, the matrix
         *  \f[
         *       B = A - \lambda I
         *  \f]
         *  has rank two. That means that its image space is two-dimensional
         *  and the desired eigenvector is orthogonal to that plane. Here,
         *  all Cartesian basis vectors are tested as pre-images and the vector
         *  \f[
         *        v = (B e_{i}) x (B e_{j})
         *  \f]
         *  with maximal length is chosen as output.
         */
        Vector constructEigenVector( const Tensor& A,
                                     const double eVal )
        {
            // identity matrix
            const Tensor I = Tensor::Identity();
            
            // rank-deficient matrix
            const Tensor B = A - eVal * I;

            // array of eigenvector candidates
            boost::array<Vector,3> candidates;
            double maxNorm = 0.;
            unsigned index = base::invalidInt;
            
            for ( unsigned i = 0; i < 3; i ++ ) {

                // span of the image space
                Vector b1, b2;
                b1.noalias() = B * I.col( i );
                b2.noalias() = B * I.col( (i+1)%3 );

                // orthogonal vector
                const Vector v  = b1.cross( b2 );
                candidates[i] = v;
                const double vNorm = base::norm( v );

                // check the norm of the orthongal vector
                if ( vNorm > maxNorm ) {
                    maxNorm = vNorm;
                    index   = i;
                }
            }

            // eigenvector
            Vector v = candidates[ index ];

            // normalise
            return v / base::norm( v );
        }
        
    }

    //--------------------------------------------------------------------------
    /** Given a symmetric tensor A, compute the eigen values and vectors.
     *  There are three essentially different cases for the set of eigen-values
     *  and the corresponding vectors. Let
     *  \f[
     *      B^i = A - \lambda^i I
     *  \f]
     *  \f$ v^i \f$ the corresponding eigen-vector. Then we have the cases:
     *    1.   All three eigenvalues are identical, \f$ B^i \f$ has rank 0,
     *         every vector is an eigenvector
     *    2.   All three eigenvalues are different. \f$ B^i \f$ has rank 2,
     *         an eigenvector can be safely constructed for all eigen values.
     *    3.   Two eigenvalues are equal, one is different. \f$ B^i \f$ has
     *         rank 1. This case is the most complicated and handled as follows:
     *         -  find different eigenvalue and compute its vector (standard)
     *         -  generate the rank-one matrix \f$ B^i \f$ for the equal values,
     *            generate the image space by taking the image \f$ B^i e_j \f$
     *            with the largest length (\f$ e_j \f$ are the Cartesian base
     *            vectors).
     *         -  construct a vector orthogonal to that image and enother one
     *            orthogonal to the image and the first one.
     *  
     */
    inline void eigenPairs( const Tensor& A,
                            boost::array<double,3>& eVal,
                            Tensor& eVec )
    {
        // compute eigen values
        eVal = eigenValues( A );

        // tolerance necessary
        const double tol = 1.e-8;

        // maximal eigenvalue
        double maxEval = 0.;
        for ( unsigned i = 0; i < 3; i++ ) {
            // sanity check
            VERIFY_MSG( (eVal[i] >= 0.),
                        "Implementation requires positive values" );
            maxEval = std::max( maxEval, eVal[i] );
        }
        
        // if all the eigenvalues are the same, return identity matrix
        if ( ( std::abs( eVal[0] - eVal[1] ) / maxEval < tol ) and
             ( std::abs( eVal[0] - eVal[2] ) / maxEval < tol ) ) {
            // v_i = e_i
            eVec = Tensor::Identity();
            return;
        }

        // check if two eigenvalues are equal and remember indices
        unsigned diffEvalIndex = base::invalidInt;
        unsigned equalIndex    = base::invalidInt;
        for ( unsigned i = 0; i < 3; i ++ ) {
            if ( std::abs(eVal[i] - eVal[(i+1)%3])/maxEval < tol ) {
                diffEvalIndex = (i+2)%3;
                equalIndex    = i;
            }
        }

        if ( diffEvalIndex < 3 ) { // 2 equal eigenvalues

            // 1) eigenvector corresponding to different eigen value
            eVec.col(0) =
                detail_::constructEigenVector( A, eVal[ diffEvalIndex ] );

            // construct rank-1 matrix
            const Tensor I = Tensor::Identity();
            const Tensor B = A - eVal[ equalIndex ] * I;

            // 2) construct vector generating the image space
            unsigned index = base::invalidInt;
            Vector x;
            {
                double maxNorm = 0.;
                boost::array<Vector,3> candidates;
                for ( unsigned i = 0; i < 3; i ++ ) {
                    // multiply rank-deficient matrix with Cartesian basis
                    Vector v;
                    v.noalias()= B * I.col( i );
                    candidates[i] = v;
                    const double vNorm = base::norm( v );

                    // memorise the maximal norm
                    if ( vNorm > maxNorm ) {
                        maxNorm = vNorm;
                        index   = i;
                    }
                }
                // image space is spanned by this vector
                x = candidates[index];
                // normalise
                x /= base::norm( x );
            }
                
            // 3) construct second vector orthogonal to x
            index = base::invalidInt;
            double minEntry = base::invalidReal();
            // find entry with smallest absolute value
            for ( unsigned i; i < 3; i++ ) {
                if ( std::abs( x[i] ) < minEntry ) {
                    index = i;
                    minEntry = std::abs( x[i] );
                }
            }

            // construct vector with the larger entries swapped
            Vector v1 = x;
            v1[ (index+1)%3 ] = - x[ (index+2)%3 ];
            v1[ (index+2)%3 ] =   x[ (index+1)%3 ];

            // 4) construct third vector with cross-product
            Vector v2 = x.cross( v1 );

            // normalise output
            eVec.col(1) = v1 / base::norm( v1 );
            eVec.col(2) = v2 / base::norm( v2 );

            return;
        }

        // use auxiliary method to construct the eigen-space
        for ( unsigned i = 0; i < 3; i ++ ) {
            eVec.col(i) = detail_::constructEigenVector( A, eVal[i] );
        }

        return;
    }
}

#endif
