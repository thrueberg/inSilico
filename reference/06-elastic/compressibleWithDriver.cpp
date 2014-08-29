#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/smf/Reader.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include <base/BoundaryValueProblem.hpp>
#include <solid/CompressibleDriver.hpp>

const double coordTol = 1.e-6;
#include "PulledSheet.hpp"

namespace ref06{
    
    int compressibleWithDriver( int argc, char * argv[] );
}


//------------------------------------------------------------------------------
int ref06::compressibleWithDriver( int argc, char * argv[] )
{
    // usage message
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0] << "  mesh.smf input.dat \n";
        return 0;
    }

    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::HyperCubeShape<SPACEDIM>::value;

    // read name of input file
    const std::string meshFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double E, nu, pull, tolerance;
    unsigned maxIter, loadSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "pull",             pull );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "loadSteps",        loadSteps );
        prop.registerPropertiesVar( "tolerance",        tolerance );

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
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
    // define a mesh
    typedef base::Unstructured<shape,geomDeg>    Mesh;
    const unsigned dim = Mesh::Node::dim;

    // create a mesh and read from input
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    //
    typedef base::BoundaryValueProblem<Mesh,dim,fieldDeg> BoundaryValueProblem;
    BoundaryValueProblem bvp( mesh );

    // typedef mat::hypel::StVenant Material;
    typedef mat::hypel::NeoHookeanCompressible Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    solid::CompressibleDriver<BoundaryValueProblem,Material>
        driver( bvp, material );

    // Dirichlet conditions
    const double firstPull = pull / static_cast<double>( loadSteps );

    driver.dirichlet(
        boost::bind(
            &ref06::PulledSheet<dim>::dirichletBC<BoundaryValueProblem::DoF>,
            _1, _2, true, firstPull ) );

    //
    for ( unsigned step = 0; step < loadSteps; step++ ) {
        std::cout << step << "  "
                  << driver.solve( step, maxIter, tolerance )
                  << std::endl;
    }
    
    return 0;
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    return ref06::compressibleWithDriver( argc, argv );
}
