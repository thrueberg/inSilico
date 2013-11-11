#include <base/io/PropertiesParser.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include "Surface.hpp"


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
template<typename SURFFIELDTUPLE>
class InternalPressure
{
public:
    
    typedef typename SURFFIELDTUPLE::GeomElement SurfaceElement;
    typedef typename SURFFIELDTUPLE::TrialElement TrialElement;
    
    typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim
    LocalVecDim;

    typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim
    GlobalVecDim;

    static const unsigned numValues = SurfaceElement::numNodes;
    typedef boost::array<GlobalVecDim,numValues> ForceArray;

    // dirty hack
    typedef SURFFIELDTUPLE  arg1_type;
    typedef ForceArray&     arg4_type;

    InternalPressure( const double p,
                      const GlobalVecDim& centre )
        : p_( p ), centre_( centre ) { }
    
    void operator()( const SURFFIELDTUPLE& sft,
                     const LocalVecDim&    xi,
                     const double          weight,
                     ForceArray& forceArray ) const
    {
        // extract surface geometry element of the tuple
        const SurfaceElement* surfEp  = sft.geomElementPtr();
        const TrialElement*   trialEp = sft.trialElementPtr();

        // Get surface metric and normal
        GlobalVecDim normal;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, xi, normal );

        // evaluate shape functions
        typename TrialElement::FEFun::FunArray trialFunValues;
        (trialEp -> fEFun()).evaluate( surfEp, xi, trialFunValues );

        // get coordinate
        const GlobalVecDim x = base::Geometry<SurfaceElement>()( surfEp, xi );

        // angle
        const double phi = std::atan2( (x[1] - centre_[1]), (x[0] - centre_[0]) );

        // force
        const GlobalVecDim force = p_ * normal * (1. + 0.8 * std::cos(2. * phi) );
        
        // assemble to array
        for ( unsigned n = 0; n < numValues; n++ )
            forceArray[n] += force * trialFunValues[n] * detG * weight;
    }
    
private:
    const double p_;
    const GlobalVecDim centre_;
};

//------------------------------------------------------------------------------
template<typename SURFFIELDTUPLE, typename QUADRATURE>
void computeForceField( SURFFIELDTUPLE sft,
                        const QUADRATURE& quadrature,
                        const double p,
                        const typename SURFFIELDTUPLE::GeomElement::Node::VecDim & centre )
{
    typedef InternalPressure<SURFFIELDTUPLE> InternalPressure;
    InternalPressure internalPressure( p, centre );

    static const unsigned dim = SURFFIELDTUPLE::GeomElement::Node::dim;

    // array of nodal forces, initialise to zero
    typename InternalPressure::ForceArray forceArray;
    forceArray.assign( base::constantVector<dim>( 0. ) );

    // compute nodal forces
    quadrature.apply( internalPressure, sft, forceArray );
    
    // add to dofs
    typedef typename SURFFIELDTUPLE::TestElement TestElement;
    TestElement* testEp = sft.testElementPtr();

    // go through degrees of freedom of surface element
    typename TestElement::DoFPtrIter dIter = testEp -> doFsBegin();
    typename TestElement::DoFPtrIter dLast = testEp -> doFsEnd();
    for ( unsigned s = 0; dIter != dLast; ++dIter, s++ ) {

        for ( unsigned d = 0; d < dim; d++ ) {
            const double oldValue = (*dIter) -> getValue( d );
            const double newValue = oldValue + forceArray[s][d];
            (*dIter) -> setValue( d, newValue );
        }
    }

    return;
}

//------------------------------------------------------------------------------
template<typename SURFACE>
void computeForceField( const SURFACE& surface,
                        const bool constrain,
                        const double p,
                        const typename SURFACE::VecDim& centre,
                        typename SURFACE::Field& forces )
{
    // extract surface current configuration
    typename SURFACE::Mesh  current;
    surface.currentConfiguration( current );
    
    // make sure the dof-pattern is the same as for the surface displacements
    if ( constrain ) constrainSupports( forces );
            
    // Number the degrees of freedom
    base::dof::numberDoFsConsecutively( forces.doFsBegin(), forces.doFsEnd() );

    // Bind the surface to the new field
    typedef typename 
        base::asmb::SurfaceFieldBinder<typename SURFACE::Mesh,
                                       typename SURFACE::Field>
        SurfaceFieldBinder;

    // tuple binder: test and trial fields are the force
    typedef typename SurfaceFieldBinder::template TupleBinder<1,1>::Type STBUU;
    SurfaceFieldBinder forceFieldBinder( current, forces );

    // go through all elements and compute the forces
    typename SurfaceFieldBinder::FieldIterator first = forceFieldBinder.elementsBegin();
    typename SurfaceFieldBinder::FieldIterator  last = forceFieldBinder.elementsEnd();
    for ( ; first != last; ++first ) {
        computeForceField( STBUU::makeTuple( *first ),
                           surface.accessQuadrature(), p,
                           centre );
        
    }

    return;
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

    Surface::Field forces, velocity;
    surface.fieldCopy( forces );
    surface.fieldCopy( velocity );

            
    surface.writeVTKFile( baseName, 0, forces, velocity );

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

        base::dof::clearDoFs( forces );

        // current load
        const double p = pressureFun( (step+1)*dt, pressure );

        // new field for surface forces
        computeForceField( surface, constrain, p, centre, forces );
        
        // solve structure problem
        const unsigned iter = surface.findEquilibrium( dt, step,
                                                       maxIter, tolerance, forces,
                                                       1.0 );
        
        // warning
        if ( iter == maxIter ) {
            std::cout << "# (WW) Step " << step << " has not converged within "
                      << maxIter << " iterations \n";
        }

        surface.computeSurfaceVelocity( velocity, dt );

        // write a vtk file
        surface.writeVTKFile( baseName, step+1, forces, velocity );
        
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
