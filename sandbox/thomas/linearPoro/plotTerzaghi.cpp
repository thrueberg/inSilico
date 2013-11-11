#include <iostream>
#include <fstream>
#include <string>

#include <base/verify.hpp>
#include <base/io/PropertiesParser.hpp>

#include "Terzaghi.hpp"

int main( int argc, char* argv[] )
{
    const std::string inputFile = "terz.inp";

    
    double E, nu, alpha, c0, k, H, F, tMax;
    unsigned numSteps, numIntervals;
    {
        base::io::PropertiesParser pp;
        pp.registerPropertiesVar( "E",     E     );
        pp.registerPropertiesVar( "nu",    nu    );
        pp.registerPropertiesVar( "alpha", alpha );
        pp.registerPropertiesVar( "c0",    c0    );
        pp.registerPropertiesVar( "k",     k     );
        pp.registerPropertiesVar( "H",     H     );
        pp.registerPropertiesVar( "F",     F     );

        pp.registerPropertiesVar( "tMax",         tMax  );
        pp.registerPropertiesVar( "numSteps",     numSteps );
        pp.registerPropertiesVar( "numIntervals", numIntervals );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str() );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        pp.readValues( inp );
        inp.close( );

        // Make sure all variables have been found
        if ( not pp.isEverythingRead() ) {
            pp.writeUnread( std::cerr );
            VERIFY_MSG( false, "Could not find above variables" );
        }

    }

    Terzaghi terzaghi( E, nu, alpha, c0, k, H, F );

    const double dx =    H / static_cast<double>( numIntervals );
    const double dt = tMax / static_cast<double>( numSteps );

    std::ofstream uOut( "u.dat" );
    std::ofstream pOut( "p.dat" );

    uOut << "# x ";
    pOut << "# x ";
    for ( unsigned it = 0; it <= numSteps; it ++ ) {
            const double t = dt * it;
            uOut << t << "  ";
            pOut << t << "  ";
    }
    uOut << '\n';
    pOut << '\n';
    

    for ( unsigned ix = 0; ix <= numIntervals; ix ++ ) {

        const double x = dx * ix;

        uOut << x << " ";
        pOut << x << " ";

        for ( unsigned it = 0; it <= numSteps; it ++ ) {

            const double t = dt * it;

            const double p = ( it == 0 ? terzaghi.initialPressure() :
                               terzaghi.pressure( x, t ) );
            const double u = ( it == 0 ? terzaghi.initialDisplacement( x ) :
                               terzaghi.displacement( x, t ) );

            uOut << u << "  ";
            pOut << p << "  ";
        }
        
        uOut << '\n';
        pOut << '\n';
    }
    
    uOut.close();
    pOut.close();
    
    return 0;
}








