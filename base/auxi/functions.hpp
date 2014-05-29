//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   functions.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_auxi_functions_hpp
#define base_auxi_functions_hpp

//------------------------------------------------------------------------------
namespace base{
    namespace auxi{

        template<typename RETURNTYPE>
        class ConstantFun;
    }
}

//------------------------------------------------------------------------------
//! Function which always returns a constant value
template<typename RETURNTYPE>
class base::auxi::ConstantFun
{
public:
    ConstantFun( const RETURNTYPE & constant )
    : constant_( constant ) { }

    template<typename ARG>
    RETURNTYPE operator()( ARG arg ) const
    {
        return constant_;
    }

    template<typename ARG1, typename ARG2>
    RETURNTYPE operator()( ARG1 arg1, ARG2 arg2 ) const
    {
        return constant_;
    }

    template<typename ARG1, typename ARG2, typename ARG3>
    RETURNTYPE operator()( ARG1 arg1, ARG2 arg2, ARG3 arg3 ) const
    {
        return constant_;
    }

private:
    const RETURNTYPE constant_;
};



#endif
