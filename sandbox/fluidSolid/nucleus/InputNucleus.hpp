#ifndef inputfsi_h
#define inputfsi_h

//------------------------------------------------------------------------------
#include <string>
#include <fstream>

#include <base/verify.hpp>
#include <base/linearAlgebra.hpp>
#include <base/io/PropertiesParser.hpp>

//------------------------------------------------------------------------------
/** Input parameters for FSI
 */
template<unsigned DIM>
struct InputNucleus
{
    static const unsigned dim = DIM;
    typedef typename base::Vector<dim,double>::Type   VecDim;
    typedef typename base::Vector<dim,unsigned>::Type VecIntDim;

    InputNucleus( const std::string& inputFile )
    {
        double bbminX, bbminY, bbminZ, bbmaxX, bbmaxY, bbmaxZ;
        double cMX, cMY, cMZ, cNX, cNY, cNZ;
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
        
        prop.registerPropertiesVar( "RN",               RN );
        prop.registerPropertiesVar( "cNX",              cNX );
        prop.registerPropertiesVar( "cNY",              cNY );
        prop.registerPropertiesVar( "cNZ",              cNZ );

        prop.registerPropertiesVar( "RM",               RM );
        prop.registerPropertiesVar( "cMX",              cMX );
        prop.registerPropertiesVar( "cMY",              cMY );
        prop.registerPropertiesVar( "cMZ",              cMZ );

        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "viscosityGel",     viscosityGel );
        prop.registerPropertiesVar( "viscosityCS",      viscosityCS );
        prop.registerPropertiesVar( "sigma",            sigma );

        prop.registerPropertiesVar( "numLoadSteps",     numLoadSteps );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "tolerance",        tolerance );
        prop.registerPropertiesVar( "alpha",            alpha );
        prop.registerPropertiesVar( "penaltyFac",       penaltyFac );
        prop.registerPropertiesVar( "dt",               dt );
        prop.registerPropertiesVar( "findTolerance",    findTolerance );
        prop.registerPropertiesVar( "readMeshFromFile", readMeshFromFile );
        prop.registerPropertiesVar( "meshFile",         meshFile );

        prop.registerPropertiesVar( "bc",               bc );
        prop.registerPropertiesVar( "ubar",             ubar );
        prop.registerPropertiesVar( "forceValue",       forceValue );


        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(),                  "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );

        N[0] = NX; bbmin[0] = bbminX; bbmax[0] = bbmaxX; centerM[0] = cMX; centerN[0] = cNX;
        if ( dim > 1 ) {
            N[1] = NY; bbmin[1] = bbminY; bbmax[1] = bbmaxY; centerM[1] = cMY; centerN[1] = cNY;
        }
        if ( dim > 2 ) {
            N[2] = NZ; bbmin[2] = bbminZ; bbmax[2] = bbmaxZ; centerM[2] = cMZ; centerN[2] = cNZ;
        }
    }

public: // on purpose

    // geomtry
    double RM, RN;
    VecDim centerM, centerN;
    VecDim bbmin, bbmax;

    // solid
    double E, nu;

    // fluids
    double viscosityGel, viscosityCS;
    double sigma;

    // stabilisation
    double alpha;

    // 
    double ubar, tolerance, penaltyFac, dt, findTolerance;
    unsigned numLoadSteps, maxIter;
    VecIntDim N;
    bool readMeshFromFile;
    BC   bc;
    double forceValue;
    std::string meshFile;
};

#endif
