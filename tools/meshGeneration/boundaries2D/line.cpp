//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   line.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <vector>
// boost includes
#include <boost/lexical_cast.hpp>
// tools include
#include <tools/meshGeneration/meshGeneration.hpp>
#include <tools/meshGeneration/boundaries2D/boundaries2D.hpp>

//------------------------------------------------------------------------------
namespace line{
    bool userInput( const int argc, char* argv[],
                    tools::meshGeneration::Point& p1,
                    tools::meshGeneration::Point& p2,
                    std::size_t& numElements )
    {
        if ( argc != 6 ) {
            std::cerr << "Usage: " << argv[0]
                      << " x1 y1  x2 y2  nInt \n\n"
                      << "To span nInt line elements between the points "
                      << "(x1,y1) and (x2,y2)\n\n";
            return false;
        }
    
        // starting point
        p1[0] = boost::lexical_cast<double>( argv[1] );
        p1[1] = boost::lexical_cast<double>( argv[2] );
        p1[2] = 0.;
        // end point
        p2[0] = boost::lexical_cast<double>( argv[3] );
        p2[1] = boost::lexical_cast<double>( argv[4] );
        p2[2] = 0.;
        // number of elements
        numElements = boost::lexical_cast<std::size_t>( argv[5] );

        return true;
    }
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // container for points and connectivities
    typedef tools::meshGeneration::Point    Point;
    typedef tools::meshGeneration::Element  Element;
    std::vector<Point>   points;
    std::vector<Element> elements;

    Point X1, X2;
    std::size_t nElem;
    const bool input = line::userInput( argc, argv, X1, X2, nElem );

    if ( not input ) return 0;

    // compute lengths
    const double lengthX = X2[0] - X1[0];
    const double lengthY = X2[1] - X1[1];
    const double hX = lengthX / static_cast<double>( nElem );
    const double hY = lengthY / static_cast<double>( nElem );

    // generate points
    for ( unsigned i = 0; i < nElem + 1; i ++ ) {
        const tools::meshGeneration::Point p = {{ X1[0] + i * hX, X1[1] + i * hY, 0. }};
        points.push_back( p );
    }

    // generate elements
    const bool closed = false;
    tools::meshGeneration::boundaries2D::generateElements( points.size(),
                                                           elements,
                                                           closed );

    // write output to stream
    tools::meshGeneration::writeSMFHeader( "line", 2,
                                           points.size(), elements.size(), std::cout );
    tools::meshGeneration::writePoints(   points,   std::cout );
    tools::meshGeneration::writeElements( elements, std::cout );

    return 0;
}
