//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Lagrange1D.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_sfun_lagrange1d_hpp
#define base_sfun_lagrange1d_hpp

//------------------------------------------------------------------------------
// base/sfun includes
#include <base/sfun/ShapeFunTraits.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace sfun{

        template<unsigned DEGREE>
        class Lagrange1D;

    }
}

//------------------------------------------------------------------------------
//! Implementation of equi-distant Lagrange interpolation functions
template<unsigned DEGREE>
class base::sfun::Lagrange1D
{
public:
    //! Template parameter: polynomial degree
    static const unsigned degree = DEGREE;

    //! For introspection
    static const unsigned dim = 1;

    //! Number of functions
    static const unsigned numFun   = degree + 1;

    //! Traits object for type defintions
    typedef base::sfun::ShapeFunResultArrays<dim,numFun,
                                             base::sfun::ScalarShapeFunResult>  SFRA;
    
    //--------------------------------------------------------------------------
    //! @name Use types from traits object
    //@{
    typedef typename SFRA::VecDim                VecDim;
    typedef typename SFRA::FunArray              FunArray;    
    typedef typename SFRA::GradArray             GradArray;   
    typedef typename SFRA::HessianArray          HessianArray;
    //@}

    //--------------------------------------------------------------------------
    //! Evaluation functions
    //@{
    //! Plain function evaluation
    void fun( const VecDim & xi, FunArray& values ) const;
    //! Evaluation of the functions' gradients
    void gradient( const VecDim & xi, GradArray& values ) const;
    //! Evaluation of the functions' Hessians
    void hessian(  const VecDim & xi, HessianArray& values ) const;
    //@}

    //--------------------------------------------------------------------------
    //! Provide the support points of the functions
    static void supportPoints( boost::array< VecDim, numFun > & supportPoints )
    {
        const double oneOverDegree = 1. / static_cast<double>( DEGREE );
        
        for ( unsigned i = 0; i < numFun; i ++ ) {
            supportPoints[i] =
                base::constantVector<1>( static_cast<double>(i) * oneOverDegree );
        }
    }
    
};

//! Include implementation file
#include <base/sfun/Lagrange1D.ipp>

#endif
