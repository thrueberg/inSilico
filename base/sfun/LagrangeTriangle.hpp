//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   LagrangeTriangle.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_sfun_lagrangetriangle_hpp
#define base_sfun_lagrangetriangle_hpp

//------------------------------------------------------------------------------
// base/sfun includes
#include <base/sfun/ShapeFunTraits.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace sfun{

        template<unsigned DEGREE>
        class LagrangeTriangle;

    }
}

//------------------------------------------------------------------------------
//! Implementation of equi-distant Lagrange interpolation functions
//! on reference triangle.
template<unsigned DEGREE>
class base::sfun::LagrangeTriangle
{
public:
    //! Template parameter: polynomial degree
    static const unsigned degree = DEGREE;

    //! For inspection
    static const unsigned dim                  = 2;
    static const base::sfun::Ordering ordering = base::sfun::HIERARCHIC;
    
    //! Number of functions
    static const unsigned numFun   = (degree+1) * (degree+2) / 2;

    //! Traits object for type defintions
    typedef base::sfun::ShapeFunResultArrays<dim,numFun,
                                             base::sfun::ScalarShapeFunResult> SFRA;
    
    //--------------------------------------------------------------------------
    //! @name Use types from traits object
    //!{
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
    static void supportPoints( boost::array< VecDim, numFun > & supportPoints );
    
};

//! Include implementation file
#include <base/sfun/LagrangeTriangle.ipp>

#endif
