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
// base   includes
#include <base/numbers.hpp>

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
    //! Compute the inner product of two vectors
    template<typename VEC>
    double dotProduct( const VEC& a, const VEC& b )
    {
        return a.dot( b );
    }

    //--------------------------------------------------------------------------
    /** Compute the Euclidean norm of a vector.
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
    //! Cross product function for the columns of 3 x 2 matrix
    base::Vector<3>::Type
    crossProduct( const base::Matrix<3,2>::Type & surfJac )
    {
        const base::Vector<3>::Type result =
            (surfJac.col(0)).cross( surfJac.col(1) );

        return result;
    }

    //! Cross product function for two vectors
    base::Vector<3>::Type
    crossProduct( const base::Vector<3>::Type & a,
                  const base::Vector<3>::Type & b )
    {
        return a.cross( b );
    }

    //! 2D cross product: gives the normal vector to the column of a 2x1 matrix
    base::Vector<2>::Type
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
    template<typename MAT>
    struct MatSize { static const int value = MAT::SizeAtCompileTime; };

    template<typename MAT>
    struct MatRows { static const int value = MAT::RowsAtCompileTime; };

    template<typename MAT>
    struct MatCols { static const int value = MAT::ColsAtCompileTime; };
    //@}
    
}

//------------------------------------------------------------------------------
#endif
