//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   linearAlgebra.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_linearalgebra_hpp
#define base_linearalgebra_hpp

//------------------------------------------------------------------------------
// eigen3 includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
// base   includes
#include <base/numbers.hpp>
#include <base/verify.hpp>

//------------------------------------------------------------------------------
namespace base{

    //--------------------------------------------------------------------------
    //! Definition of a fixed-size vector
    template<unsigned SIZE, typename NUM=double>
    struct Vector
    {
        typedef Eigen::Matrix<NUM,SIZE,1> Type;
    };

    //! Generate a vector with constant entries
    template<unsigned SIZE>
    static typename Vector<SIZE>::Type constantVector( const double value )
    {
        return Vector<SIZE,double>::Type::Constant( value  );
    }

    //! Function which generates an invalidated vector
    template<unsigned SIZE>
    static typename Vector<SIZE>::Type invalidVector()
    {
        return Vector<SIZE,double>::Type::Constant( base::invalidReal() );
    }

    //--------------------------------------------------------------------------
    //! Definition of a fixed-size matrix
    template<unsigned SIZE1, unsigned SIZE2, typename NUM=double>
    struct Matrix
    {
        typedef Eigen::Matrix<NUM,SIZE1,SIZE2> Type;
    };

    //! Generate matrix with constant value in all entries
    template<unsigned SIZE1, unsigned SIZE2>
    static typename Matrix<SIZE1,SIZE2>::Type constantMatrix( const double value )
    {
        return Matrix<SIZE1,SIZE2,double>::Type::Constant( value );
    }

    //--------------------------------------------------------------------------
    //! Definition of a dynamic vectors and matrices
#ifdef INSILICO_COMPLEX
    typedef Eigen::VectorXcd  VectorD;
    typedef Eigen::MatrixXcd  MatrixD;
#else
    typedef Eigen::VectorXd   VectorD;
    typedef Eigen::MatrixXd   MatrixD;
#endif

    //--------------------------------------------------------------------------
    //! Cross product function for the columns of 3 x 2 matrix
    inline base::Vector<3>::Type
    crossProduct( const base::Matrix<3,2>::Type & surfJac )
    {
        const base::Vector<3>::Type result =
            (surfJac.col(0)).cross( surfJac.col(1) );

        return result;
    }

    //! Cross product function for two vectors
    inline base::Vector<3>::Type
    crossProduct( const base::Vector<3>::Type & a,
                  const base::Vector<3>::Type & b )
    {
        return a.cross( b );
    }

    //! 2D cross product: gives the normal vector to the column of a 2x1 matrix
    inline base::Vector<2>::Type
    crossProduct( const base::Matrix<2,1>::Type & surfJac )
    {
        base::Vector<2>::Type result;
        result[0] =  surfJac( 1, 0 );
        result[1] = -surfJac( 0, 0 );
        return result;
    }

    //--------------------------------------------------------------------------
    //! @name  Block operations
    //@{
    template<unsigned SIZE,typename VEC>
    typename base::Vector<SIZE,typename VEC::Scalar>::Type tail( const VEC& vec )
    {
        return vec.template tail<SIZE>();
    }

    template<unsigned SIZE,typename VEC>
    typename base::Vector<SIZE,typename VEC::Scalar>::Type head( const VEC& vec )
    {
        return vec.template head<SIZE>();
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Traits
    //@{
    //! Compile-time size of matrix
    template<typename MAT>
    struct MatSize { static const int value = MAT::SizeAtCompileTime; };
    //! Compile-time number of rows of matrix
    template<typename MAT>
    struct MatRows { static const int value = MAT::RowsAtCompileTime; };
    //! Compile-time number of columns of matrix
    template<typename MAT>
    struct MatCols { static const int value = MAT::ColsAtCompileTime; };
    //@}

    //--------------------------------------------------------------------------
    //! Convert matrix into one row; useful for output
    template<typename MAT>
    typename base::Matrix<1,MAT::SizeAtCompileTime, typename MAT::Scalar>::Type
    makeRow( const MAT& mat )
    {
        typename base::Vector<MAT::SizeAtCompileTime, typename MAT::Scalar>::Type
            result;
        for ( unsigned r = 0; r < MAT::RowsAtCompileTime; r++ ) {
            for ( unsigned c = 0; c < MAT::ColsAtCompileTime; c++ ) {
                result[ r * MAT::ColsAtCompileTime + c ] = mat( r, c );
            }
        }

        return result;
    }

    //--------------------------------------------------------------------------
    //! Compute the inner product of two vectors / matrices
    template<typename VEC1, typename VEC2>
    double dotProduct( const VEC1& a, const VEC2& b )
    {
        double result = 0.;
        for ( int i = 0; i < base::MatRows<VEC1>::value; i++ )
            for ( int j = 0; j < base::MatCols<VEC1>::value; j++ )
                result += a(i,j) * b(i,j);
        
        return result;
    }

    //--------------------------------------------------------------------------
    /** Compute the Euclidean norm of a vector / Frobenius norm of a matrix
     *  Definition
     *  \f[
     *      \| a \|_2 = ( \sum_{i=1}^d a_i^2 )^{1/2}
     *  \f]
     *  \tparam VEC  Type of vector
     *  \param[in] a Input vector
     *  \return      Euclidean norm of a
     */
    template<typename VEC>
    double norm( const VEC& a )
    {
        return std::sqrt( dotProduct( a, a ) );
    }


    //--------------------------------------------------------------------------
    /** Compute all eigenvalues of a given symmetric matrix.
     *  \tparam MAT Type of matrix to operate on (size <= 3)
     */
    template<typename MAT>
    typename base::Vector<MatRows<MAT>::value>::Type eigenValues( const MAT& A )
    {
        static const unsigned numRows = MatRows<MAT>::value;
        
        STATIC_ASSERT_MSG( numRows == MatCols<MAT>::value, "Matrix must be square" );
        STATIC_ASSERT_MSG( numRows <= 3, "Direct computation only for N<=3" );
        
        Eigen::SelfAdjointEigenSolver<MAT> es( numRows );
        es.computeDirect( A, Eigen::EigenvaluesOnly );
        const typename base::Vector<numRows>::Type result = es.eigenvalues();
        return result;
    }

    //--------------------------------------------------------------------------
    /** Compute all eigenvalues and eigenvectors of a given symmetric matrix.
     *  \tparam MAT Type of matrix to operate on (size <= 3)
     */
    template<typename MAT>
    typename base::Vector<MatRows<MAT>::value>::Type eigenPairs( const MAT& A,
                                                                 MAT& X )
    {
        static const unsigned numRows = MatRows<MAT>::value;
        
        STATIC_ASSERT_MSG( numRows == MatCols<MAT>::value, "Matrix must be square" );
        STATIC_ASSERT_MSG( numRows <= 3, "Direct computation only for N<=3" );
        
        Eigen::SelfAdjointEigenSolver<MAT> es( numRows );
        es.computeDirect( A, Eigen::ComputeEigenvectors );
        const typename base::Vector<numRows>::Type result = es.eigenvalues();
        X = es.eigenvectors();
        return result;
    }

    //--------------------------------------------------------------------------
    namespace detail_{

        template<unsigned SIZE> struct ComputeAngleAxis;
        
        //! Helper to compute the axis-angle representation of 2x2 rotation matrix
        template<>
        struct ComputeAngleAxis<2>
        {
            template<typename MAT>
            static std::pair<double,base::Vector<3>::Type> apply( const MAT& R )
            {
                // trivial axis: in 2D rotate always around z-axis
                base::Vector<3>::Type axis;
                axis[0] = axis[1] = 0.; axis[2] = 1.;
                // determine angle from first column
                double angle = std::atan2( R(1,0), R(0,0) );
                // ensure positive angle
                if ( angle < 0. ) angle += M_PI;
                return std::make_pair( angle, axis );
            }
        };

        //! Helper to compute the axis-angle representation of 3x3 rotation matrix
        template<>
        struct ComputeAngleAxis<3>
        {
            template<typename MAT>
            static std::pair<double,base::Vector<3>::Type> apply( const MAT& R )
            {
                // http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
                base::Vector<3>::Type axis;
                axis[0] = R( 2, 1 ) - R( 1, 2 );
                axis[1] = R( 0, 2 ) - R( 2, 0 );
                axis[2] = R( 1, 0 ) - R( 0, 1 );
                const double r = axis.norm();

                double angle = 0.;
                if ( r >= 1.e-10 ) {
                    // in case of a zero axis (diagonal matrix, hence no rotation)
                    axis /=  r;
                    angle = std::atan2( r, (R(0,0) + R(1,1) + R(2,2)) - 1. );
                }

                // ensure positive angle
                if ( angle < 0. ) angle += M_PI;
                
                return std::make_pair( angle, axis );
            }
        };

    }

    //--------------------------------------------------------------------------
    /** Given a rotation matrix, determine the angle axis representation.
     *  Let a matrix \f$ R \f$ be given which represents a (2D/3D) rotation.
     *  There is an alternative representation by axis and angle, see
     *  http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
     *  This object generates this representation.
     *  \tparam MAT  Type of the matrix to analyse
     */
    template<typename MAT>
    std::pair<double,base::Vector<3>::Type> angleAxis( const MAT& R )
    {
        static const unsigned numRows = MatRows<MAT>::value;
        STATIC_ASSERT_MSG( (numRows <= 3) and (numRows > 1),
                           "Matrix size does not make sense" );
        return detail_::ComputeAngleAxis<numRows>::apply( R );
    }

}

//------------------------------------------------------------------------------
#endif
