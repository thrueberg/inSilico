//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   torus.cpp
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


//------------------------------------------------------------------------------
namespace torus{
    
    bool userInput( const int argc, char* argv[],
                    double& radius1, double& radius2, 
                    tools::meshGeneration::Point& centre,
                    std::size_t& numElements1,
                    std::size_t& numElements2,
                    bool& makeTriangles,
                    std::string& message )
    {
        if ( ( argc != 8 ) and ( argc != 9 ) ){
            std::cerr << "Usage: " << argv[0]
                      << " r1 r2 cx cy cz n1 n2 [mt=0]\n\n"
                      << "To create a torus with radii r1 and r2 (r1 > r2) around "
                      << "centre (cx,cy,cz) with a grid of n1 x n2 elements.\n"
                      << "The optional flag mt (default=false=0) forces the "
                      << "generation of triangles instead of quadrilaterals\n\n";
            return false;
        }

        // radii
        radius1 = boost::lexical_cast<double>( argv[1] );
        radius2 = boost::lexical_cast<double>( argv[2] );

        // a word to the user
        if ( not( radius1 > radius2) ) {
            std::cerr << "(WW) The major radius is NOT greater than the minor radius.\n"
                      << "(WW) The torus will be degenerate. \n\n";
        }
    
        // centre point
        centre[0] = boost::lexical_cast<double>( argv[3] );
        centre[1] = boost::lexical_cast<double>( argv[4] );
        centre[2] = boost::lexical_cast<double>( argv[5] );
    
        // number of elements
        numElements1 = boost::lexical_cast<std::size_t>( argv[6] );
        numElements2 = boost::lexical_cast<std::size_t>( argv[7] );

        // make triangles
        makeTriangles = ( argc == 8 ? false :
                          boost::lexical_cast<bool>( argv[8] ) );

        // concatenate the input arguments
        for ( int c = 0; c < argc; c++ ) {
            message += argv[c];
            message += " ";
        }

        return true;
    }
}


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // container for points and connectivities
    typedef tools::meshGeneration::Point    Point;
    typedef tools::meshGeneration::Element  Element;

    double radius1, radius2;    // radii of the torus
    Point centre;               // centre of the torus
    std::size_t nElem1, nElem2; // grid of nElem1 x nElem2
    bool makeTriangles;         // flag to decide element shape
    std::string message;        // message for trace back of output
    const bool input = torus::userInput( argc, argv, radius1, radius2, centre,
                                         nElem1, nElem2, makeTriangles, message );
    if ( not input ) return 0;  // unsuccessful input

    // angle increments
    const double dPhi   = 2. * M_PI / static_cast<double>( nElem1 );
    const double dTheta = 2. * M_PI / static_cast<double>( nElem2 );

    // generate points
    std::vector<Point>   points;
    for ( std::size_t j = 0; j < nElem2; j++ ) {
        for ( std::size_t i = 0; i < nElem1; i++ ) {
            Point p = centre;
            const double iD = static_cast<double>( i );
            const double jD = static_cast<double>( j );
            p[0] += ( radius1 + radius2 * std::cos(jD * dTheta) )* cos(iD * dPhi );
            p[1] += ( radius1 + radius2 * std::cos(jD * dTheta) )* sin(iD * dPhi );
            p[2] +=             radius2 * std::sin(jD * dTheta);
            points.push_back( p );
        }
    }

    // generate elements
    std::vector<Element> elements;

    for ( std::size_t j = 0; j < nElem2; j++ ) {
        for ( std::size_t i = 0; i < nElem1; i++ ) {

            // following indices with boundary periodicity
            const std::size_t iNext = (i < nElem1 - 1? i+1 : 0);
            const std::size_t jNext = (j < nElem2 - 1? j+1 : 0);
                
            // the four vertices of a quadrilateral
            const std::size_t i1 = j     * nElem1 + i;
            const std::size_t i2 = j     * nElem1 + iNext;
            const std::size_t i4 = jNext * nElem1 + i;
            const std::size_t i3 = jNext * nElem1 + iNext;

            if ( makeTriangles ) {
                // generate two triangles
                tools::meshGeneration::Element tri1, tri2;
                tri1.push_back( i1 );
                tri1.push_back( i2 );
                tri1.push_back( i3 );
                tri2.push_back( i3 );
                tri2.push_back( i4 );
                tri2.push_back( i1 );
                elements.push_back( tri1 );
                elements.push_back( tri2 );
                
            }
            else {
                // generate a quadrilateral
                tools::meshGeneration::Element quad;
                quad.push_back( i1 );
                quad.push_back( i2 );
                quad.push_back( i3 );
                quad.push_back( i4 );
                elements.push_back( quad );
            }
        }
    }

    // write output to stream
    const std::string elementShape = ( makeTriangles ? "triangle" : "quadrilateral" );
    const unsigned    elementNumPoints = ( makeTriangles ? 3 : 4 );

    tools::meshGeneration::writeSMFComment( message, std::cout );
    tools::meshGeneration::writeSMFHeader( elementShape, elementNumPoints, 
                                           points.size(), elements.size(), std::cout );
    tools::meshGeneration::writePoints(   points,   std::cout );
    tools::meshGeneration::writeElements( elements, std::cout );

    return 0;
}
