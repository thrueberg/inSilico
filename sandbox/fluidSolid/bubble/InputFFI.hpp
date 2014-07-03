#ifndef inputffi_h
#define inputffi_h

//------------------------------------------------------------------------------
#include <string>
#include <fstream>

#include <base/verify.hpp>
#include <base/linearAlgebra.hpp>
#include <base/io/PropertiesParser.hpp>

//------------------------------------------------------------------------------
/** Input parameters for FFI
 */
template<unsigned DIM>
struct InputFFI
{
    static const unsigned dim = DIM;
    typedef typename base::Vector<dim,double>::Type   VecDim;
    typedef typename base::Vector<dim,unsigned>::Type VecIntDim;

    InputFFI( const std::string& inputFile )
    {
        double bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ;
        double cX, cY, cZ;
        unsigned NX, NY, NZ;
        
         //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "bbminX",           bbminX );
        prop.registerPropertiesVar( "bbminY",           bbminY );
        prop.registerPropertiesVar( "bbminZ",           bbminZ );
        prop.registerPropertiesVar( "bbmaxX",           bbmaxX );
        prop.registerPropertiesVar( "bbmaxY",           bbmaxY );
        prop.registerPropertiesVar( "bbmaxZ",           bbmaxZ );
        prop.registerPropertiesVar( "NX",               NX );
        prop.registerPropertiesVar( "NY",               NY );
        prop.registerPropertiesVar( "NZ",               NZ );
        
        prop.registerPropertiesVar( "R",                R );
        prop.registerPropertiesVar( "cX",               cX );
        prop.registerPropertiesVar( "cY",               cY );
        prop.registerPropertiesVar( "cZ",               cZ );
        
        prop.registerPropertiesVar( "viscosity1",       viscosity1 );
        prop.registerPropertiesVar( "rho1",             rho1 );
        prop.registerPropertiesVar( "viscosity2",       viscosity2 );
        prop.registerPropertiesVar( "rho2",             rho2 );
        prop.registerPropertiesVar( "sigma",            sigma );


        prop.registerPropertiesVar( "numLoadSteps",     numLoadSteps );
        prop.registerPropertiesVar( "alpha",            alpha );
        prop.registerPropertiesVar( "penaltyFac",       penaltyFac );
        prop.registerPropertiesVar( "dt",               dt );
        prop.registerPropertiesVar( "readMeshFromFile", readMeshFromFile );
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "fluidIter",        fluidIter );

        prop.registerPropertiesVar( "bc",               bc );
        prop.registerPropertiesVar( "ubar",             ubar );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(),                  "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );

        N[0] = NX; bbmin[0] = bbminX; bbmax[0] = bbmaxX; center[0] = cX;
        if ( dim > 1 ) { N[1] = NY; bbmin[1] = bbminY; bbmax[1] = bbmaxY; center[1] = cY;}
        if ( dim > 2 ) { N[2] = NZ; bbmin[2] = bbminZ; bbmax[2] = bbmaxZ; center[2] = cZ;}
    }

public: // on purpose

    // read from input file
    double R, rho1, rho2, ubar, viscosity1, viscosity2, penaltyFac, dt, sigma;
    double alpha;
    VecDim center;
    VecDim bbmin, bbmax;
    unsigned numLoadSteps, fluidIter;
    VecIntDim N;
    bool readMeshFromFile;
    BC   bc;
    std::string meshFile;

};

#endif
