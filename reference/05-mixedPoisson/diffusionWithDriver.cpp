#include <fstream>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>

#include <base/BoundaryValueProblem.hpp>
#include <heat/DiffusionDriver.hpp>

#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>

#include <base/time/BDF.hpp>

#include <base/asmb/SimpleIntegrator.hpp>
#include <base/kernel/FieldIntegral.hpp>

#include <base/dof/location.hpp>
#include <base/post/Monitor.hpp>

const double coordTol = 1.e-5;


//------------------------------------------------------------------------------
namespace ref05{

    //--------------------------------------------------------------------------
    /** Prescribe a DIM-dimensional Gaussian with center at 0.5
     */
    template<unsigned DIM, typename DOF>
    void gaussian( const typename base::Vector<DIM,double>::Type& x, DOF* doFPtr )
    {
        // domain centre
        const typename base::Vector<DIM,double>::Type c =
            base::constantVector<DIM>( 0.5 );

        const double sigma = 0.05;

        double value = 1.;
        double norm  = 1.;
        for ( unsigned d = 0; d < DIM; d++ ) {
            const double arg = (x[d] - c[d])*(x[d] - c[d]);
            const double alpha = 0.5 / sigma / sigma;
            value *= std::exp( -arg  / 2. / sigma /sigma );
            norm  *= alpha;
        }

        const double fac =
            std::sqrt( base::Power<DIM>::apply( M_PI ) / norm );

        doFPtr -> setValue( 0, value / fac );
    }

    //--------------------------------------------------------------------------
    //! Integrate over field to get the total mass
    template<typename BVP>
    double totalMass( typename BVP::Mesh& mesh, BVP& bvp )
    {
        typename BVP::VecDoF intU = base::constantVector<BVP::doFSize>( 0. );
        base::Quadrature<3,BVP::shape> quadrature;
        typename BVP::FieldBinder fieldBinder( mesh, bvp.getField() );
        base::asmb::simplyIntegrate<typename BVP::UU>(
            quadrature, intU, fieldBinder,
            base::kernel::FieldIntegral<typename BVP::UU::Tuple>()  );

        return intU[0];
    }

    //--------------------------------------------------------------------------
    int diffusionWithDriver( int argc, char * argv[] );
}

//------------------------------------------------------------------------------
/** Solve a diffusion problem with driver.
 *
 *  \image html concentration.gif  "Diffusion process"
 */
int ref05::diffusionWithDriver( int argc, char * argv[] )
{
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf input.dat \n\n";
        return -1;
    }

    // Input argument: the mesh file in smf format
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    // create the basename 
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    // Input data file
    const std::string inpFile = boost::lexical_cast<std::string>( argv[2] );

    // User input parameter
    unsigned numSteps;  // Number of time steps
    double stepSize;    // Size of a time step
    {
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "numSteps",         numSteps );
        prop.registerPropertiesVar( "stepSize",         stepSize );

        // Read variables from the input file
        std::ifstream inp( inpFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        prop.readValues( inp );
        inp.close( );

        // Make sure all variables have been found
        if ( not prop.isEverythingRead() ) {
            prop.writeUnread( std::cerr );
            VERIFY_MSG( false, "Could not find above variables" );
        }
    }

    //--------------------------------------------------------------------------
    // Define static parameters
    const unsigned    geomDeg  = 1;  // degree of geometry approximation
    const unsigned    fieldDeg = 2;  // degree of field approximation
    const base::Shape shape    = base::QUAD; // element shape
    typedef base::time::BDF<1> MSM;  // type of time integration = Euler BW

    //--------------------------------------------------------------------------
    // Define a mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;

    // create mesh object, read from file
    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    //--------------------------------------------------------------------------
    // Use driver class for the field
    typedef base::BoundaryValueProblem<Mesh,1,fieldDeg,MSM> BVP;
    BVP bvp( mesh );

    //--------------------------------------------------------------------------
    // Diffusion driver
    typedef heat::DiffusionDriver<BVP> Diffusion;
    Diffusion diffusion( bvp );
    diffusion.initialCondition( boost::bind( gaussian<BVP::dim,
                                             BVP::DoF>, _1, _2 ) );

    //--------------------------------------------------------------------------
    // initial values
    bvp.writeVTKFile( baseName + "." + base::io::leadingZeros( 0 ) );

    // find a point in the field
    const std::pair<std::size_t,BVP::VecDim> monitorLocation
        = base::dof::findPointInField( mesh, bvp.getField(),
                                       base::constantVector<BVP::dim>( 0.25 ),
                                       1.e-10 );

    // set up a Monitor 
    base::post::Monitor<Mesh::Element,BVP::Field::Element>
        monitor( mesh.elementPtr(              monitorLocation.first ),
                 bvp.getField().elementPtr( monitorLocation.first ),
                 monitorLocation.second );

    // write initial state
    std::cout << "#time   value(0.25)   total mass \n"
              << "0.    ";
    monitor.solution( std::cout );
    std::cout << " " << ref05::totalMass( mesh, bvp ) << std::endl;

    //--------------------------------------------------------------------------
    // Do a time loop
    for ( unsigned step = 0; step < numSteps; step++ ) {

        // do one time step
        diffusion.advanceInTime( step, stepSize );

        // write a solution file
        bvp.writeVTKFile( baseName + "." +
                                     base::io::leadingZeros( step+1 ) );

        // write some values to stdout
        std::cout << step * stepSize << "  ";
        monitor.solution( std::cout );
        std::cout << " " << ref05::totalMass( mesh, bvp )
                  << std::endl;
    }


    return 0; // bye-bye
}

//------------------------------------------------------------------------------
// Delegate
int main( int argc, char * argv[] )
{
    return ref05::diffusionWithDriver( argc, argv );
}
