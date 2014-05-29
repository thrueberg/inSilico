#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/io/PropertiesParser.hpp>
#include <base/linearAlgebra.hpp>
#include <base/verify.hpp>

#include "FluidBox.hpp"
#include "Surface.hpp"
#include "InternalPressure.hpp"

//------------------------------------------------------------------------------
// Driven cavity BC
template<unsigned DIM, typename DOF>
void drivenCavity( const typename base::Vector<DIM>::Type& x,
                   DOF* doFPtr,
                   const std::vector<double>&  supports ) 
{
    const double tol = 1.e-5;

    // if d-th coordinate has the value 1.0
    bool onLid = ( std::abs( x[DIM-1] - 1.0 ) < tol );
    // remove the corner/edge locations
    for ( unsigned d = 0; d < DIM-1; d ++ ) {
        if ( std::abs( x[d] - 0.0 ) < tol ) onLid = false;
        if ( std::abs( x[d] - 1.0 ) < tol ) onLid = false;
    }

    const double supportSize = supports[ doFPtr -> getID() ];

    // boundary condition is either 0 or the e_1 vector 
    for ( unsigned d = 0; d < DIM; d ++ ) {
        const double value = ( ( (d==0) and onLid ) ? 1.0 : 0.0 ) * supportSize;
        
        if ( doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, value );
    }
}

//------------------------------------------------------------------------------
// Box BC
template<unsigned DIM, typename DOF>
void fixedBottom( const typename base::Vector<DIM>::Type& x,
                  DOF* doFPtr,
                  const std::vector<double>&  supports ) 
{
    const double tol = 1.e-5;

    // if d-th coordinate has the value 0.
    const bool onBottom = ( std::abs( x[DIM-1] - 0.0 ) < tol );
    const bool onTop    = ( std::abs( x[DIM-1] - 1.0 ) < tol );

    const double supportSize = supports[ doFPtr -> getID() ];

    // boundary condition is either 0 or the e_1 vector
    if ( onBottom or onTop ) {
        for ( unsigned d = 0; d < DIM; d ++ ) {
            const double value =
                (d == 0 ? 1.0 : 0.0) * (onTop ? 1.0 : -1.0) * supportSize;
        
            if ( doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, value );
        }
    }
}
//[dirichlet]}

template<unsigned DIM>
base::Vector<1>::Type
source( const typename base::Vector<DIM>::Type& x, const double value ) 
{
    base::Vector<1>::Type result;

    result[0] = value;
    return result;
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    const unsigned    dim       = 2;

    //--------------------------------------------------------------------------
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input.dat \n\n";
        return -1;
    }

    const unsigned fsiIter = 1;
    
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    //--------------------------------------------------------------------------
    std::string meshFile, surfFile;
    double viscosity, density, tolerance, penaltyFactor, A, B, dt, robinFac, aux, pressure;
    unsigned maxIter, numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "surfFile",         surfFile );
        prop.registerPropertiesVar( "viscosity",        viscosity );
        prop.registerPropertiesVar( "density",          density );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "tolerance",        tolerance );
        prop.registerPropertiesVar( "penaltyFactor",    penaltyFactor );
        prop.registerPropertiesVar( "robinFac",         robinFac );
        prop.registerPropertiesVar( "aux",              aux );
        prop.registerPropertiesVar( "pressure",         pressure );

        prop.registerPropertiesVar( "A", A );
        prop.registerPropertiesVar( "B", B );
        prop.registerPropertiesVar( "dt",    dt );
        prop.registerPropertiesVar( "numSteps", numSteps );

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

    const std::string baseName = base::io::baseName( meshFile, ".smf" );
    
    //--------------------------------------------------------------------------
    // Construct two fluid solvers
    typedef FluidBox<dim> FluidBox;
    const bool inside = false;
    FluidBox fluidBoxOut( density, viscosity, inside );
    fluidBoxOut.readMesh( meshFile );

    FluidBox fluidBoxIn( density, viscosity, not inside );
    fluidBoxIn.readMesh( meshFile );

    //--------------------------------------------------------------------------
    // Construct a surface solver
    const double h = base::mesh::Size<FluidBox::Mesh::Element>::apply(
        fluidBoxIn.accessMesh().elementPtr(0) );
    const double nitscheWeight = penaltyFactor * viscosity / h;

    const double robin = nitscheWeight * robinFac;
    const double robinForce = robin * aux;

    typedef Surface<dim> Surface;
    Surface surface( A, B, robin );
    surface.readMesh( surfFile );
    surface.numberDoFs();

    // no convection used here
    const bool withConvection = false;

    double volume = surface.enclosedVolume();
    Surface::VecDim centroid = surface.moment() / volume;

    double oldVolume            = volume;
    Surface::VecDim oldCentroid = centroid;

    for ( unsigned step = 0; step < numSteps; step++ ) {

        
        
        std::cout << "Step: "<< step;

        for ( unsigned iter = 0; iter < fsiIter; iter++ ) {
            std::cout << " " << iter; 

            base::dof::clearDoFs( surface.accessForces() );
            
            //------------------------------------------------------------------
            // Get current surface configuration
            Surface::Mesh current;
            surface.currentConfiguration( current );

            // Immerse
            fluidBoxOut.immerseSurface( current );

            // Apply box-boundary conditions to outer domain
            fluidBoxOut.constrainVelocity( boost::bind( &fixedBottom<dim,
                                                        FluidBox::Velocity::DegreeOfFreedom>,
                                                        _1, _2,
                                                        boost::ref( fluidBoxOut.supportSizes() )));

            // Fix a pressure dof
            FluidBox::VecDim fixPPoint = base::constantVector<dim>( 0. );
            //fluidBoxOut.fixPressureDoF( fixPPoint, 1.e-8, 0.0 );

            //------------------------------------------------------------------
            std::pair<std::size_t, std::size_t> numDoFs = fluidBoxOut.numberDoFs();
            //std::cout << "Out: " << numDoFs.first << ", " << numDoFs.second << "\n";


            // solve fluid problem in outer domain
            fluidBoxOut.solveProblem( surface.accessMesh(),
                                      surface.accessVelocity(),
                                      boost::bind( &source<dim>, _1, 0. ),
                                      penaltyFactor,
                                      maxIter, tolerance, withConvection );
            
            fluidBoxOut.computeImmersedSurfaceForce( current,
                                                     surface.accessForces(),
                                                     viscosity, robinForce );
            //------------------------------------------------------------------
#if 0
            // Immerse inner domain
            fluidBoxIn.immerseSurface( current );

            // point where to fix the pressure (better centroid of the surface)
            // Apply a positive value for the surface not to collapse
            fixPPoint = centroid;
            fluidBoxIn.fixPressureDoF( fixPPoint, 2.e-2, 0.0 );

            numDoFs = fluidBoxIn.numberDoFs();
            //std::cout << "In: " << numDoFs.first << ", " << numDoFs.second << "\n";

            // solver inner fluid problem
            fluidBoxIn.solveProblem( surface.accessMesh(),
                                     surface.accessVelocity(),
                                     boost::bind( &source<dim>, _1, 0.00 ),
                                     penaltyFactor,
                                     maxIter, tolerance, withConvection );

            fluidBoxIn.computeImmersedSurfaceForce( current,
                                                    surface.accessForces(),
                                                    viscosity, robinForce );
#endif
            //------------------------------------------------------------------
            // solve the surface

            computeForceField( surface, pressure, surface.accessForces() );

            // solve structure problem
            surface.findEquilibrium( dt, step, maxIter, tolerance );

            surface.computeSurfaceVelocity( dt );


            oldVolume   = volume;
            oldCentroid = centroid;
            
            // write out volume and centroid
            volume = surface.enclosedVolume();
            centroid = (1./volume) * surface.moment();

            const double volConv = std::abs( volume - oldVolume );
            const double cenConv = (centroid - oldCentroid).norm();

            std::cout << " " << volume << " " << centroid.transpose()
                      << " " << volConv << "  " << cenConv
                      << "\n";

            if ( (volConv*volConv) + (cenConv*cenConv) < tolerance*tolerance ) break;

        } // iterations

        surface.updateHistory();
        
        fluidBoxOut.writeVTKFile( baseName + ".out", step );
        fluidBoxIn.writeVTKFile( baseName + ".in", step );
        surface.writeVTKFile( baseName + ".surf", step );
        
    } // time steps
    
    return 0;
}
