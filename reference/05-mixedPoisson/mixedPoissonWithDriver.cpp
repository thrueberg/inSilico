#include <fstream>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>

#include <base/BoundaryValueProblem.hpp>
#include <heat/PoissonDriver.hpp>

#include <base/io/smf/Reader.hpp>

const double coordTol = 1.e-5;

#include "ReferenceSolution.hpp"
#include "GivenData.hpp"

//------------------------------------------------------------------------------
namespace ref05{
    int mixedPoissonWithDriver( int argc, char * argv[] );
}

//------------------------------------------------------------------------------
/** Solve a boundary value problem with driver.
 *  Same problem is solved as in ref05::mixedPoisson(), but this time the
 *  heat::Driver class is used and the convenience function
 *  heat::solveMixedPoissonProblem().
 *
 *  New features are
 *  - use of the driver class heat::Driver
 *
 */
int ref05::mixedPoissonWithDriver( int argc, char * argv[] )
{
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf \n\n";
        return -1;
    }

    // Input argument: the mesh file in smf format
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    // create the basename 
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    //--------------------------------------------------------------------------
    // Define parameters
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;

    //--------------------------------------------------------------------------
    // Define a mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;

    // creat mesh object, read from file
    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    //--------------------------------------------------------------------------
    // Use driver class
    typedef base::BoundaryValueProblem<Mesh,1,fieldDeg> BoundaryValueProblem;
    BoundaryValueProblem bvp( mesh );

    //--------------------------------------------------------------------------
    // Define reference solution and BVP manager
    typedef ReferenceSolution<BoundaryValueProblem::dim> RefSol;
    typedef GivenData<RefSol>                 GivenData;
    RefSol    refSol( 3., 5., 4. );
    GivenData givenData( refSol );

    heat::PoissonDriver<BoundaryValueProblem> ppd( bvp, 1. );
    ppd.dirichlet( boost::bind( &GivenData::dirichletBC<BoundaryValueProblem::DoF>,
                                &givenData, _1, _2 ) );
    ppd.neumann(   boost::bind( &GivenData::neumannBC, &givenData, _1, _2 ) );
    ppd.force(     boost::bind( &GivenData::forceFun,  &givenData, _1 ) );
    ppd.solve();
    
    //--------------------------------------------------------------------------
    // write solution
    bvp.writeVTKFile( "test" );
    
    return 0;
}

//------------------------------------------------------------------------------
// Delegate
int main( int argc, char * argv[] )
{
    return ref05::mixedPoissonWithDriver( argc, argv );
}
