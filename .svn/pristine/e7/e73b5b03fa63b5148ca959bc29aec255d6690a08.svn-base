#include <iostream>
#include <boost/lexical_cast.hpp>
#include <base/sfun/BSpline.hpp>

#include <docs/plotFun/FunEvalPolicy.hpp>

//------------------------------------------------------------------------------
// Evaluate the desired spline regularity
template<template<unsigned,unsigned> class FUN, unsigned DEGREE, int CTR>
struct EvaluateRegularity
{
    static void apply( const unsigned derivative,
                       const unsigned regularity, const unsigned resolution,
                       std::ostream& out )
    {
        // Go through available degrees
        if ( regularity == CTR )
            docs::plotFun::EvaluateFunDerivative<FUN<DEGREE,CTR> >::apply(derivative,
                                                                          resolution, out );
        else
            EvaluateRegularity<FUN,DEGREE,CTR-1>::apply( derivative, regularity, resolution, out );
    }
};

//------------------------------------------------------------------------------
// Error dummy not to be called
template<template<unsigned,unsigned> class FUN, unsigned DEGREE>
struct EvaluateRegularity<FUN,DEGREE,-1>
{
    static void apply( const unsigned derivative,
                       const unsigned regularity, const unsigned resolution,
                       std::ostream& out )
    {
        std::cerr << "Regularity = " << regularity << ", degree = " << DEGREE << std::endl;
        VERIFY_MSG( false, "You tried to call with a wrong regularity" );
    }
};

//------------------------------------------------------------------------------
// Evaluate the desired spline degree
template<template<unsigned,unsigned> class FUN, int CTR>
struct EvaluateDegrees
{
    static void apply( const unsigned degree,     const unsigned derivative,
                       const unsigned regularity, const unsigned resolution,
                       std::ostream& out )
    {
        // Go through available degrees
        if ( degree == CTR )
            EvaluateRegularity<FUN,CTR,CTR>::apply( derivative, regularity,
                                                    resolution, out );
        else
            EvaluateDegrees<FUN,CTR-1>::apply( degree, derivative,
                                               regularity, resolution, out );
    }
};

//------------------------------------------------------------------------------
// Error dummy not to be called
template<template<unsigned,unsigned> class FUN>
struct EvaluateDegrees<FUN,-1>
{
    static void apply( const unsigned degree,     const unsigned derivative,
                       const unsigned regularity, const unsigned resolution,
                       std::ostream& out )
    {
        std::cerr << "Degree = " << degree << std::endl;
        VERIFY_MSG( false, "You tried to call with a wrong degree" );
    }
};

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    if ( argc != 5 ) {
        std::cout << "Usage:  " << argv[0]
                  << " degree  regularity derivative resolution \n";
        return 0;
    }

    // read from user input
    const unsigned degree     = boost::lexical_cast<unsigned>( argv[1] );
    const unsigned regularity = boost::lexical_cast<unsigned>( argv[2] );
    const unsigned derivative = boost::lexical_cast<unsigned>( argv[3] );
    const unsigned resolution = boost::lexical_cast<unsigned>( argv[4] );

    // maximally available degree pre-defined
    const int maxDegree = 4;

    // act
    EvaluateDegrees< base::sfun::BSpline, maxDegree >::apply( degree, derivative,
                                                              regularity, resolution,
                                                              std::cout );

    return 0;
}
