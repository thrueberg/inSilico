//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FunEvaluationPolicy.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_auxi_funevaluationpolicy_hpp
#define base_auxi_funevaluationpolicy_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/geometry.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace auxi{


        //----------------------------------------------------------------------
        /** Evaluate a function which receices a physical coordinate.
         *  Given a function with the signature
         *  \code{.cpp}
         *  result_type fun( const GlobalVecDim& x )
         *  \endcode
         *  this functor wraps this function into the form
         *  \code{.cpp}
         *  result_type fun( const GEOMELEMENT* ep, const LocalVecDim& xi,
         *                   const FUN& fFUn )
         *  \endcode
         *  in which the geometry element and a local coordinate \f$ \xi \f$ are
         *  used to compute the physical coordinate \f$ x = x(\xi) \f$ and
         *  evaluate the given function \a fFun at this coordinate. The types
         *  \a GlobalVecDim and \a LocalVecDim are derived from the type of the
         *  geometry element.
         *  \tparam GEOMELEMENT  Type of geometry element
         *  \tparam FUN          Type of function/functor to evaluate
         */
        template<typename GEOMELEMENT, typename FUN>
        struct EvaluateDirectly
        {
            typedef FUN Fun;
            typedef typename FUN::result_type result_type;

            typedef typename base::GeomTraits<GEOMELEMENT>::LocalVecDim  LocalVecDim;
            
            static result_type apply( const GEOMELEMENT* gep,
                                      const LocalVecDim& xi, 
                                      const Fun& fFun )
            {
                return fFun( base::Geometry<GEOMELEMENT>()( gep, xi ) );
            }
        };

        //----------------------------------------------------------------------
        /** Evaluate a function which receives an element pointer and a local
         *  coordinate. 
         *  In order to have the same signature as EvaluateDirectly, this object
         *  simply performs the evaluation of a given function of type
         *  \code{.cpp}
         *  result_type fun( const GEOMELEMENT* ep, const LocalVecDim& xi )
         *  \endcode
         *  \tparam GEOMELEMENT  Type of geometry element
         *  \tparam FUN          Type of function/functor to evaluate
         */
        template<typename GEOMELEMENT, typename FUN>
        struct EvaluateViaElement
        {
            typedef FUN Fun;
            typedef typename FUN::result_type result_type;

            typedef typename base::GeomTraits<GEOMELEMENT>::LocalVecDim LocalVecDim;

            static result_type apply( const GEOMELEMENT* gep,
                                      const LocalVecDim& xi, 
                                      const Fun& fFun )
            {
                return fFun( gep, xi );
            }
        };

    }
}
#endif
