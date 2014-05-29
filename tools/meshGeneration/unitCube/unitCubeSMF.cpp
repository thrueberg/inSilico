//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   unitCubeSMF.cpp
//! @author Thomas Rueberg
//! @date   2013

// std includes
#include <iostream>
#include <vector>
// local include
#include <tools/meshGeneration/unitCube/unitCube.hpp>

//------------------------------------------------------------------------------
/** Generate an equi-distant mesh on a unit hypercube.
 *  The user provides the number of elements per direction (the number of these
 *  numbers implies the dimension) and this tool writes the generated mesh in
 *  SMF format to stdout.
 */
int main( int argc, char* argv[] )
{
    namespace unitCube = tools::meshGeneration::unitCube;

#ifdef SIMPLEX
    const bool makeSimplices = true;
#else
    const bool makeSimplices = false;
#endif

#ifdef QUADRATIC
    const unsigned degree = 2;
#else
    const unsigned degree = 1;
#endif

    // dimension and number of elements per direction
    unsigned dim, e1, e2, e3;
    const bool input = 
        unitCube::userInput( argc, argv, dim, e1, e2, e3 );
    if ( not input ) return 0;
    
    unitCube::unitCubeSMF<makeSimplices,degree>( dim, e1, e2, e3, std::cout );
    return 0;
}
