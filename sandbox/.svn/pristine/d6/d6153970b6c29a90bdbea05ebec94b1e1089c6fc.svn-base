#include <base/io/PropertiesParser.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include "Surface.hpp"
#include "InternalPressure.hpp"

//------------------------------------------------------------------------------
// Statically determined positioning
template<typename FIELD>
void constrainSupports( FIELD& field )
{
    VERIFY_MSG( (FIELD::DegreeOfFreedom::size==2),
                "This method only works for 2D" );
    
    typename FIELD::DoFPtrIter dIter = field.doFsBegin();

    const std::size_t numDoFs = std::distance( dIter, field.doFsEnd() );

    VERIFY_MSG( (numDoFs % 4) == 0, "Trying to fix the quarter points" );

    (*dIter) -> constrainValue( 1, 0. );

    std::advance( dIter, numDoFs/4 );
    (*dIter) -> constrainValue( 0, 0. );
        
    std::advance( dIter, numDoFs/4 );
    (*dIter) -> constrainValue( 1, 0. );

    std::advance( dIter, numDoFs/4 );
    (*dIter) -> constrainValue( 0, 0. );
}


//------------------------------------------------------------------------------
double pressureFun( const double t, const double max )
{
    if ( t < 1.0 )
        return t * max;

    return max;
}
 //------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    const unsigned    dim      = 2;

    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  input.dat \n";
        return 0;
    }

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double A, B, gamma, dt, tolerance, pressure;
    unsigned maxIter, numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "A",                A );
        prop.registerPropertiesVar( "B",                B );
        prop.registerPropertiesVar( "gamma",            gamma );
        prop.registerPropertiesVar( "dt",               dt );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "numSteps",         numSteps );
        prop.registerPropertiesVar( "tolerance",        tolerance );
        prop.registerPropertiesVar( "pressure",         pressure );

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

    // if there is no mass of the surface, have to apply constraints
    const bool constrain = (gamma == 0.);
    
    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    typedef Surface<dim> Surface;
    Surface surface( A, B, gamma );

    surface.readMesh( meshFile );

    if ( constrain ) constrainSupports( surface.accessDisplacements() );

    surface.numberDoFs();

    surface.writeVTKFile( baseName, 0);

    double          volume = surface.enclosedVolume();
    Surface::VecDim centre = surface.moment() / volume;

    std::cout << "# step  time pressure   radius   centre\n"
              << " 0  0.0   0.0  "
              << std::pow( volume, 1./static_cast<double>(dim))
              << "   " << centre.transpose()
              << "\n";
    
    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < numSteps; step++ ) {

        base::dof::clearDoFs( surface.accessForces() );

        // current load
        const double p = pressureFun( (step+1)*dt, pressure );

        // new field for surface forces
        if ( constrain ) constrainSupports( surface.accessForces() );
        computeInternalPressure( surface, p, surface.accessForces() );
        
        // solve structure problem
        const unsigned iter = surface.findEquilibrium( dt, step,
                                                       maxIter, tolerance );
        
        // warning
        if ( iter == maxIter ) {
            std::cout << "# (WW) Step " << step << " has not converged within "
                      << maxIter << " iterations \n";
        }

        surface.computeSurfaceVelocity( dt );

        // write a vtk file
        surface.writeVTKFile( baseName, step+1 );
        
        volume = surface.enclosedVolume();
        centre = surface.moment() / volume;
        std::cout << step+1 << "  " << (step+1)*dt << "  " << p << " "
                  << std::pow( volume, 1./static_cast<double>(dim))
                  << "   " << centre.transpose()
                  << "\n";


        surface.updateHistory();
    }
    // Finished load steps
    //--------------------------------------------------------------------------

    
    return 0;
}
