//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   circle/circle.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <vector>
#include <string>
// boost includes
#include <boost/lexical_cast.hpp>
// base includes
#include <base/verify.hpp>
// tools include
#include <tools/meshGeneration/meshGeneration.hpp>

//------------------------------------------------------------------------------
namespace circle{

    // get data from the command line
    bool userInput( const int argc, char* argv[],
                    double& radius,
                    double& length,
                    std::size_t& N,
                    std::size_t& M,
                    std::string& message )
    {
        if ( argc != 5 ) {
            std::cerr << "Usage: "
                      << argv[0] << " R L N M \n\n"
                      << "R - Radius of circle around (0,0) \n"
                      << "L - Section length for transition square-circle\n"
                      << "N - Number of elements along the length L\n"
                      << "M - Number of elements along a quarter circle\n"
                      << std::endl;
            return false;
        }

        // read from input arguments
        radius = boost::lexical_cast<double     >( argv[1] );
        length = boost::lexical_cast<double     >( argv[2] );
        N      = boost::lexical_cast<std::size_t>( argv[3] );
        M      = boost::lexical_cast<std::size_t>( argv[4] );

        // assert positivity of radius and length
        VERIFY_MSG( radius > 0., "Enter positive radius" );
        VERIFY_MSG( length > 0., "Enter positive section length" );

        // assert positivity of radius and length
        VERIFY_MSG( (M > 0) and (N > 0),
                    "Enter positive number of elements" );

        // section length must be less than R / sqrt(2)
        if ( length >= radius / std::sqrt(2.) ) {
            std::cerr << ""
                      << std::endl;
            return false;
        }

        // concatenate the input arguments
        for ( int c = 0; c < argc; c++ ) {
            message += argv[c];
            message += " ";
        }


        return true;
    }

}

//------------------------------------------------------------------------------
//  Generate a quadrilateral mesh in a circle of radius R around (0,0).
//  The mesh is structured into 5 regions:
//  - a square of size 2L x 2L centered around (0,0)
//  - four segments between each of the square's sides and a quarter circle
//  The square is triangulated by 2N x 2N squares and the segments are
//  triangulated by M elements in radial and 2N elements in tangential
//  direction.
//  In polar coordinates the outer boundary of the square is something like
//  L / cos(phi) and this radius is linear interpolated in radial direction
//  to the true radius of the outer circle boundary.
//  By varying the parameter L (0 < L < R/sqrt(2) ) and the numbers N and M,
//  the user can try to adjust the element distortion to his taste.
//
int main( int argc, char* argv[] )
{
    // get all data from caller
    double radius, length;
    std::size_t N, M;
    std::string message;
    if ( not circle::userInput( argc, argv, radius, length, N, M, message ) )
        return 1;

    // point and element containers
    typedef tools::meshGeneration::Point   Point;
    typedef tools::meshGeneration::Element Element;
    std::vector<Point>   points;
    std::vector<Element> elements;

    //--------------------------------------------------------------------------
    // central square

    // coordinates
    const double h = length / static_cast<double>( N );

    for ( std::size_t j = 0; j <= 2*N; j++ ) {
        for ( std::size_t i = 0; i <= 2*N; i++ ) {
            Point p;
            p[0] = -length + h * static_cast<double>( i );
            p[1] = -length + h * static_cast<double>( j );
            p[2] = 0.0;
            points.push_back( p );
        }
    }

    // connectivity
    for ( std::size_t j = 0; j < 2*N; j++ ) {
        for ( std::size_t i = 0; i < 2*N; i++ ) {

            Element element;
            const std::size_t v = i + j * (2 * N + 1);

            element.push_back( v     );
            element.push_back( v + 1 );
            element.push_back( v + 1 + (2*N+1) );
            element.push_back( v +     (2*N+1) );

            elements.push_back( element );
        }
    }

    //--------------------------------------------------------------------------
    // circle segments

    // coordinates
    const double dPhi = M_PI / 4. / static_cast<double>( N );

    for ( unsigned s = 0; s < 4; s ++ ) {

        // starting angle of the segment
        const double phi0 = M_PI / 4. * static_cast<double>( 2.*s -1. );
        for ( std::size_t j = 0; j < 2*N; j++ ) {

            // current angle
            const double psi = static_cast<double>( j ) * dPhi;
            const double phi = phi0 + psi;
            
            // radius at the square boundary
            const double rho = length / std::cos( psi - M_PI/4. );
            for ( std::size_t i = 0; i < M; i++ ) {
                
                // interpolation coordinate
                const double xi = static_cast<double>(i+1) / static_cast<double>(M);

                // radius linearly interpolated
                const double r  = (1.-xi) * rho + xi * radius;

                // angle at the square boundary
                const double phi2 =
                    atan2( -length + h * static_cast<double>(j),
                           length ) + s * M_PI/2.;

                // interpolate angles
                const double alpha =
                    (1.+xi*xi-2.*xi) * phi2 + (2.*xi-xi*xi) * phi;

                // point via polar coordinates
                Point p;
                p[0] = r * std::cos( alpha );
                p[1] = r * std::sin( alpha );
                p[2] = 0.0;
                
                points.push_back( p );
            }
        }
    }

    // connectivity

    // total number of points in the square
    const std::size_t last = (2*N+1) * (2*N+1);

    // corner indices and increments for the square boundary points
    const boost::array<std::size_t,4> corners =
        {{ 2*N, last-1, last-1-2*N, 0 }};
    const boost::array<std::size_t,4> shifts = {{ 2*N+1, 1, 2*N+1, 1 }};
    const boost::array<int,4> signs = {{ 1, -1, -1, 1 }};

    // go through segments
    for ( unsigned s = 0; s < 4; s++ ) {

        // first index (-1) of new segment node
        const std::size_t first = last + s * (2 * M * N) - 1;

        for ( std::size_t j = 0; j < 2*N; j++ ) {
            for ( std::size_t i = 0; i < M; i++ ) {

                Element element;

                // element vertices
                const std::size_t v1 =
                    ( i == 0 ?
                      corners[s] + j * signs[s] * shifts[s] :
                      first + i + j * M );

                const std::size_t v2 = first + i + j * M + 1;

                // trial indices (changed if end of circle is reached)
                const std::size_t v3t = v2 + M;

                const std::size_t v4t =
                    ( i == 0 ?
                      corners[s] + (j+1) * signs[s] * shifts[s] :
                      v3t - 1);

                // special case: close the ring
                const std::size_t v4 =
                    ( s==3 and j==2*N-1 and i>0 ? last - 1 + i : v4t );
                
                const std::size_t v3 =
                    ( s==3 and j==2*N-1 ?         last     + i : v3t );

                // store vertex indices
                element.push_back( v1 );
                element.push_back( v2 );
                element.push_back( v3 );
                element.push_back( v4 );

                elements.push_back( element );
                
                
            }
        }
    }
    
    // write output to stream
    tools::meshGeneration::writeSMFComment( message, std::cout );
    tools::meshGeneration::writeSMFHeader( "quadrilateral", 4,
                                           points.size(), elements.size(), std::cout );
    tools::meshGeneration::writePoints(    points,   std::cout );
    tools::meshGeneration::writeElements(  elements, std::cout );

    return 0;
}
