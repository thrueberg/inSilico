//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   kernel/KernelFun.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_kernel_kernelfun_hpp
#define base_kernel_kernelfun_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace kernel{

        //----------------------------------------------------------------------
        /** Define the type of function a kernel comprises.
         *  This function type is useful for cases of kernels with operator()
         *  overloading. 
         */
        template<typename FIELDTUPLE, typename RESULT>
        struct KernelFun
        {
            typedef typename FIELDTUPLE::GeomElement GeomElement;
            typedef typename base::GeomTraits<GeomElement>::LocalVecDim LocalVecDim;

            typedef boost::function<void( const FIELDTUPLE&,const LocalVecDim&,
                                          const double,RESULT& )>          Type;
        };
        
    }
}

#endif
