#ifndef helper_hpp
#define helper_hpp

#include <string>
#include <fstream>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <base/linearAlgebra.hpp>
#include <base/io/PropertiesParser.hpp>


//------------------------------------------------------------------------------
//  Bock of material, occupying (0,1)^DIM, fix x_1 = 0 and pull at x_1 = 1.
//  Optionally, at x_1 = 1, a surface traction is applied or a normal
//  displacement.
template<unsigned DIM>
class PulledSheetProblem
{
public:
    typedef typename base::Vector<DIM>::Type VecDim;

    // Fix x_0=0 and optionally pull at x_1=1
    template<typename DOF>
    static void dirichletBC( const VecDim& x, DOF* doFPtr,
                             const double value ) 
    {
        // tolerance for coordinate identification
        const double tol = 1.e-5;

        // location at x_1 = 0 or x_1 = 1
        const bool onLeftBdr = ( std::abs( x[0] -  0. ) < tol );
        const bool onRightBdr = ( std::abs( x[0] -  1. ) < tol );

        // Fix left boundary at x_0 = 0
        if ( onLeftBdr ) {
            for ( unsigned d = 0; d < DOF::size; d++ ) {
                if ( doFPtr -> isActive(d) )
                    doFPtr -> constrainValue( d, 0.0 );
            }
        }

        // If assked for, apply normal displacement at x_1=1
        if (  onRightBdr ) { 
            if ( doFPtr -> isActive(0) )
                doFPtr -> constrainValue( 0, value );
            
        }
        
        return;
    }

};


//------------------------------------------------------------------------------
bool userInput( const int argc, char * argv[],
                std::string& meshFile, std::string& baseName, 
                double& E, double& nu, double& pull, double& tolerance,
                unsigned& maxIter, unsigned& loadSteps,
                bool& updated )
{
     // usage message
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0] << "mesh.smf  input.dat \n";
        return false;
    }

    // read name of input file
    meshFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string inputFile = boost::lexical_cast<std::string>( argv[2] );

    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "pull",             pull );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "loadSteps",        loadSteps );
        prop.registerPropertiesVar( "tolerance",        tolerance );
        prop.registerPropertiesVar( "updated",          updated );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        prop.readValues( inp );
        inp.close( );

        // Make sure all variables have been found
        if ( not prop.isEverythingRead() ) {
            prop.writeUnread( std::cerr );
            VERIFY_MSG( false, "Could not find above variables" );
        }
    }

    // find base name from mesh file
    baseName = base::io::baseName( meshFile, ".smf" );

    return true;;
}

#endif
