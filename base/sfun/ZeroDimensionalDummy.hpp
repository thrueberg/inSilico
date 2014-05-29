//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ZeroDimensionalDummy.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_sfun_zerodimensionaldummy_hpp
#define base_sfun_zerodimensionaldummy_hpp

//------------------------------------------------------------------------------
// base/sfun includes
#include <base/sfun/ShapeFunTraits.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace sfun{

        template<unsigned DEGREE>
        class ZeroDimensionalDummy;

    }
}

//------------------------------------------------------------------------------
template<unsigned DEGREE>
class base::sfun::ZeroDimensionalDummy
{
public:
    //! Template parameter: polynomial degree
    static const unsigned degree = DEGREE;

    //! For inspection
    static const unsigned dim = 0;

    //! Number of functions
    static const unsigned numFun   = 1;

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
    void fun( const VecDim & xi, FunArray& values ) const
    {
        values[0] = 1.;
        return;
    }
    //! Evaluation of the functions' gradients
    void gradient( const VecDim & xi, GradArray& values ) const
    {
        return;
    }
   
    //! Evaluation of the functions' Hessians
    void hessian(  const VecDim & xi, HessianArray& values ) const
    {
        return;
    }
    //@}

    //--------------------------------------------------------------------------
    //! Provide the support points of the functions
    static void supportPoints( boost::array< VecDim, numFun > & supportPoints )
    {
        return;
    }
    
};

#endif
