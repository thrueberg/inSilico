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
    struct VectorType
    {
        typedef Eigen::Matrix<NUM,SIZE,1> Type;
    };

    //! Generate a vector with constant entries
    template<unsigned SIZE>//, typename NUM = double>
    static typename VectorType<SIZE>::Type constantVector( const double value )
    {
        typedef double NUM;
        return VectorType<SIZE,NUM>::Type::Constant( value  );
    }

    //! Function which generates an invalidated vector
    template<unsigned SIZE>//, typename NUM = double>
    static typename VectorType<SIZE>::Type invalidVector()
    {
        typedef double NUM;
        return VectorType<SIZE,NUM>::Type::Constant( base::invalidReal() );
    }

    //--------------------------------------------------------------------------
    //! Definition of a fixed-size matrix
    template<unsigned SIZE1, unsigned SIZE2, typename NUM=double>
    struct MatrixType
    {
        typedef Eigen::Matrix<NUM,SIZE1,SIZE2> Type;
    };

    //! Generate matrix with constant value in all entries
    template<unsigned SIZE1, unsigned SIZE2>
    static typename MatrixType<SIZE1,SIZE2>::Type constantMatrix( const double value )
    {
        return MatrixType<SIZE1,SIZE2>::Type::Constant( value );
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
    typename base::VectorType<3>::Type
    crossProduct( const typename base::MatrixType<3,2>::Type & surfJac )
    {
        typename base::VectorType<3>::Type result =
            (surfJac.col(0)).cross( surfJac.col(1) );

        return result;
    }

    //! 2D cross product: gives the normal vector to the column of a 2x1 matrix
    typename base::VectorType<2>::Type
    crossProduct( const typename base::MatrixType<2,1>::Type & surfJac )
    {
        typename base::VectorType<2>::Type result;
        result[0] =  surfJac( 1, 0 );
        result[1] = -surfJac( 0, 0 );
        return result;
    }

    
}

//------------------------------------------------------------------------------
#endif
