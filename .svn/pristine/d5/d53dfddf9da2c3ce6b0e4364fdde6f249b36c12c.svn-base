//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ShapeFunTraits.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_sfun_shapefuntraits_hpp
#define base_sfun_shapefuntraits_hpp

//------------------------------------------------------------------------------
// boost  includes
#include <boost/array.hpp>
// base   includes
#include <base/linearAlgebra.hpp>


namespace base{
    namespace sfun{

        //! Result types for scalar shape functions
        template<unsigned DIM>
        struct ScalarShapeFunResult
        {
            typedef double                                   Fun;
            typedef typename base::Vector<DIM>::Type     Grad;
            typedef typename base::Matrix<DIM,DIM>::Type Hessian;
        };

        //! Result types for vectorial shape functions (e.g. Raviart-Thomas)
        template<unsigned DIM>
        struct VectorShapeFunResult
        {
            typedef typename base::Vector<DIM>::Type     Fun;
            typedef typename base::Matrix<DIM,DIM>::Type Grad;
            typedef typename boost::array<Grad,DIM>          Hessian;
        };

        //----------------------------------------------------------------------
        //! Array types for the shape functions
        template<unsigned DIM, unsigned NFUN,
                 template <unsigned> class SHAPEFUNRESULT>
        struct ShapeFunResultArrays
        {
            //! Evaluation point coordinate
            typedef typename base::Vector<DIM>::Type      VecDim;

            //! Traits object for result type declarations
            typedef  SHAPEFUNRESULT<DIM>                      SFR;

            //! @name Result container types
            //@{
            typedef boost::array<typename SFR::Fun, NFUN>    FunArray;
            typedef boost::array<typename SFR::Grad,NFUN>    GradArray;
            typedef boost::array<typename SFR::Hessian,NFUN> HessianArray;
            //@}
            
        };

    }
}


#endif
