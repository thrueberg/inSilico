#include <iostream>
#include <boost/lexical_cast.hpp>
#include <base/sfun/Lagrange1D.hpp>

#include <docs/plotFun/FunEvalPolicy.hpp>

//------------------------------------------------------------------------------
template<template<unsigned> class FUN, int CTR>
struct EvaluateDegrees
{
    static void apply( const unsigned degree,     const unsigned derivative,
                       const unsigned resolution, std::ostream& out )
    {
        // Go through available degrees
        if ( degree == CTR )
            docs::plotFun::EvaluateFunDerivative<FUN<CTR> >::apply(derivative,
                                                                   resolution, out );
        else
            EvaluateDegrees<FUN,CTR-1>::apply( degree, derivative,
                                               resolution, out );
    }
};

//------------------------------------------------------------------------------
template<template<unsigned> class FUN>
struct EvaluateDegrees<FUN,-1>
{
    static void apply( const unsigned degree,     const unsigned derivative,
                       const unsigned resolution, std::ostream& out )
    {
        std::cerr << "Degree = " << degree << std::endl;
        VERIFY_MSG( false, "You tried to call with a wrong degree" );
    }
};

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    if ( argc != 4 ) {
        std::cout << "Usage:  " << argv[0] << " degree  derivative resolution \n";
        return 0;
    }
    
    const unsigned degree     = boost::lexical_cast<unsigned>( argv[1] );
    const unsigned derivative = boost::lexical_cast<unsigned>( argv[2] );
    const unsigned resolution = boost::lexical_cast<unsigned>( argv[3] );

    const int maxDegree = 3;

    EvaluateDegrees< base::sfun::Lagrange1D, maxDegree >::apply( degree, derivative,
                                                                 resolution, std::cout );

    return 0;
}
