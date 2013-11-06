//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ellipse.cpp
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
namespace circle{
    bool userInput( const int argc, char* argv[],
                    double& radius1, double& radius2,
                    tools::meshGeneration::Point& centre,
                    std::size_t& numElements,
                    std::string& message )
    {
        if ( argc != 6 ) {
            std::cerr << "Usage: " << argv[0]
                      << " r1 r2 cx cy nInt \n\n"
                      << "To create an ellipse radii r1 & r2 with nInt line "
                      << "elements around the centre (cx,cy)\n\n";
            return false;
        }

        // radius
        radius1 = boost::lexical_cast<double>(   argv[1] );
        radius2 = boost::lexical_cast<double>(   argv[2] ); 
    
        // centre point
        centre[0] = boost::lexical_cast<double>( argv[3] );
        centre[1] = boost::lexical_cast<double>( argv[4] );
        centre[2] = 0.;
        // number of elements
        numElements = boost::lexical_cast<std::size_t>( argv[5] );

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

    double radius1, radius2;
    Point centre;
    std::size_t nElem;
    std::string message;
    const bool input =
        circle::userInput( argc, argv, radius1, radius2, centre, nElem, message);
    if ( not input ) return 0;

    // angle increment
    const double dPhi = 2. * M_PI / static_cast<double>(nElem);

    // generate points
    std::vector<Point>   points;
    const std::size_t nNodes    = nElem;
    for ( std::size_t i = 0; i < nNodes; i++ ) {
        Point p;
        p[0] = centre[0] + radius1 * cos( static_cast<double>(i) * dPhi );
        p[1] = centre[1] + radius2 * sin( static_cast<double>(i) * dPhi );
        p[2] =  0.0;
        points.push_back( p );
    }

    // write calling string as comment to output
    tools::meshGeneration::writeSMFComment( message, std::cout );

#if 0
    //--------------------------------------------------------------------------
    // let user know geometric items
    const double area = radius * radius * M_PI;
    const double polar = M_PI/2. * radius * radius * radius * radius;

    std::cout << "# --------------------------------------------------" << std::endl
              << "#  Area:         " << area << std::endl
              << "#  Centroid:     (" << centre[0] << ", " << centre[1] << ")" << std::endl
              << "#  Polar moment: " << polar << std::endl
              << "# --------------------------------------------------" << std::endl;
#endif

    // generate elements
    const bool closed = true;
    std::vector<Element> elements;
    tools::meshGeneration::boundaries2D::generateElements( points.size(), elements, closed );

    // write output to stream
    tools::meshGeneration::writeSMFHeader( "line", 2,
                                           points.size(), elements.size(), std::cout );
    tools::meshGeneration::writePoints(   points,   std::cout );
    tools::meshGeneration::writeElements( elements, std::cout );

    return 0;
}
