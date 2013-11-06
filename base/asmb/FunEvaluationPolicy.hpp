//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FunEvaluationPolicy.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_funevaluationpolicy_hpp
#define base_asmb_funevaluationpolicy_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace asmb{


        //----------------------------------------------------------------------
        template<typename GEOMELEMENT, typename FUN>
        struct EvaluateDirectly
        {
            typedef FUN Fun;
            typedef typename FUN::result_type result_type;

            typedef typename base::GeomTraits<GEOMELEMENT>::LocalVecDim  LocalVecDim;
            typedef typename base::GeomTraits<GEOMELEMENT>::GlobalVecDim GlobalVecDim;
            
            static result_type apply( const GEOMELEMENT* gep,
                                      const LocalVecDim& xi, 
                                      const Fun& fFun )
            {
                return fFun( base::Geometry<GEOMELEMENT>()( gep, xi ) );
            }
        };

        //----------------------------------------------------------------------
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
