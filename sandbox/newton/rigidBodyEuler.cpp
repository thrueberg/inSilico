#include <string>
#include <boost/lexical_cast.hpp>

#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/Field.hpp>

#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>
#include <base/io/smf/Writer.hpp>

#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/bruteForce.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/ScaledField.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/generateSurfaceMesh.hpp>

#include <base/fe/Basis.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/associateLocation.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/BodyForce.hpp>

#include <base/time/BDF.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>

#include "generateMesh.hpp"
#include "findLocation.hpp"

#include "extrapolateToOutside.hpp"
#include "advectField.hpp"
#include "moveSurface.hpp"

const double coordTol = 1.e-6;

//------------------------------------------------------------------------------
// Level set function for a domain covering the interval a < x_1 < b
template<unsigned DIM>
bool interval( const typename base::Vector<DIM>::Type& x,
                     typename base::Vector<DIM>::Type& xClosest,
               const double a,
               const double b )
{
    if ( x[0] < a ) { xClosest[0] = a; return false; }

    if ( x[0] > b ) { xClosest[0] = b; return false; }

    const double mid = (a + b) / 2.;
    if (x[0] < mid) xClosest[0] = a;
    else            xClosest[0] = b;
    
    return true;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM>::Type forceFun( const typename base::Vector<DIM>::Type& x,
                                           const double value )
{
    typename base::Vector<DIM>::Type f = base::constantVector<DIM>( 0. );        
    f[0] = value;
    return f;
}


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    if ( argc != 3 ) {
        std::cerr << "Usage: " << argv[0] << " N  input.dat\n\n";
        return -1;
    }
    
    // read name of input file
    const unsigned    numElements = boost::lexical_cast<unsigned>(    argv[1] );
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double xmin, xmax, a, b, rho, f, stepSize;
    unsigned numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "xmin",     xmin );
        prop.registerPropertiesVar( "xmax",     xmax );
        prop.registerPropertiesVar( "a",        a );
        prop.registerPropertiesVar( "b",        b );
        prop.registerPropertiesVar( "rho",      rho );
        prop.registerPropertiesVar( "f",        f );
        prop.registerPropertiesVar( "stepSize", stepSize );
        prop.registerPropertiesVar( "numSteps", numSteps );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );
    }

    // spatial dimension
    const unsigned    dim = 1; 
    
    // basic attributes of the computation
    const unsigned             geomDeg  = 1;
    const unsigned             fieldDeg = 1;
    const base::Shape             shape = base::HyperCubeShape<dim>::value;
    const unsigned    kernelDegEstimate = 3;
    const unsigned              doFSize = dim;
    const unsigned              tiOrder = 1;

    // choose a time stepping method
    typedef base::time::BDF<tiOrder> MSM;

    // time stepping method determines the history size
    const unsigned nHist = MSM::numSteps;

    typedef base::Unstructured<shape,geomDeg>  Mesh;
    typedef Mesh::Node::VecDim VecDim;
    Mesh mesh;
    {
        generateMesh( mesh, numElements, xmin, xmax );
    }

    //--------------------------------------------------------------------------
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::analyticLevelSet( mesh,
                                 boost::bind( &interval<dim>, _1, _2, a, b ),
                                 true, levelSet );

    //--------------------------------------------------------------------------
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells;
    base::cut::generateCutCells( mesh, levelSet, cells );

    //--------------------------------------------------------------------------
    // FE
    typedef base::fe::Basis<shape,fieldDeg>               FEBasis;
    typedef base::cut::ScaledField<FEBasis,doFSize,nHist> Field;
    typedef Field::DegreeOfFreedom                        DoF;
    Field displacement, velocity;
    base::dof::generate<FEBasis>( mesh, displacement );
    base::dof::generate<FEBasis>( mesh, velocity );

    // find geometry association for the dofs
    std::vector<std::pair<std::size_t,VecDim> > doFLocation;
    base::dof::associateLocation( displacement, doFLocation );

    // for domain fields
    typedef base::asmb::FieldBinder<Mesh,Field,Field> FieldBinder;
    FieldBinder fieldBinder(  mesh, displacement, velocity );
    typedef FieldBinder::TupleBinder<1,1>::Type UU;
    typedef FieldBinder::TupleBinder<1,2>::Type UV;
    typedef FieldBinder::TupleBinder<2,1>::Type VU;
    typedef FieldBinder::TupleBinder<2,2>::Type VV;

    typedef base::kernel::Mass<UV::Tuple> Mass;
    Mass mass( -rho );

    writeVTKFile( "euler", 0, mesh, displacement, /*velocity,*/ levelSet );

    VecDim left, right;

    //--------------------------------------------------------------------------
    // * Loop over time steps *
    // Solve the system of ODE's
    //
    //     | \dot{u} |   |0 -rho| |u|   |0|
    // rho |         | + |      |*| | = | |
    //     | \dot{v} |   |0    0| |v|   |f|
    //
    //--------------------------------------------------------------------------
    const double supportThreshold = 1.e-10;
    for ( unsigned step = 0; step < numSteps; step++ ) {

        //----------------------------------------------------------------------
        // 0) Solve Lagrangian problem
        displacement.activateAll();
        velocity.activateAll();
    
        typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
        CutQuadrature cutQuadrature( cells, true );

        // compute supports, scale basis
        std::vector<double> supports;
        base::cut::supportComputation( mesh, displacement, cutQuadrature, supports );

        displacement.scaleAndTagBasis( supports,  supportThreshold );
        velocity.scaleAndTagBasis(     supports,  supportThreshold );
    
        // number DoFs (first U, then V)
        const std::size_t activeDoFsU = 
            base::dof::numberDoFsConsecutively( displacement.doFsBegin(),
                                                displacement.doFsEnd() );
        const std::size_t activeDoFsV = 
            base::dof::numberDoFsConsecutively( velocity.doFsBegin(),
                                                velocity.doFsEnd(),
                                                activeDoFsU );
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( activeDoFsU + activeDoFsV );

        // Body force
        const double fFun = f * std::cos( step * stepSize );
        base::asmb::bodyForceComputation<VV>( cutQuadrature, solver, fieldBinder,
                                              boost::bind( &forceFun<dim>, _1, fFun ) );

        // Inertia terms due to time integration
        base::time::computeInertiaTerms<VV,MSM>( cutQuadrature, solver, fieldBinder,
                                                 stepSize, step, rho );
        base::time::computeInertiaTerms<UU,MSM>( cutQuadrature, solver, fieldBinder,
                                                 stepSize, step, rho );

        // Stiffness matrix
        base::asmb::stiffnessMatrixComputation<UV>( cutQuadrature, solver,
                                                    fieldBinder, mass );

        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.superLUSolve();
            
        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, displacement );
        base::dof::setDoFsFromSolver( solver, velocity );

        displacement.unscaleValues( supports, supportThreshold );
        velocity.unscaleValues(     supports, supportThreshold );

        // write a vtk file
        writeVTKFile( "euler", step+1, mesh, displacement, /*velocity,*/ levelSet );


        // 1) Pass Data to complementary domain
        {
            const double h = (xmax - xmin) / numElements; // mesh size
            const double extL = h; // extrapolation length

            // pass extrapolated solution to inactive DoFs
            extrapolateToOutside( mesh, displacement, extL, doFLocation, levelSet );
            extrapolateToOutside( mesh, velocity,     extL, doFLocation, levelSet );
        }
        
        // 2) Geometry update
        {
            
            // Extract Surface mesh
            typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
            SurfaceMesh immersedSurface;
            base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, immersedSurface );

            // Move with displacement solution
            moveSurface( mesh, displacement, immersedSurface );
            
            // get location of domain boundary
            SurfaceMesh::NodePtrIter surfNIter = immersedSurface.nodesBegin();
            (*surfNIter) -> getX( &(left[0] ) );
            surfNIter++;
            (*surfNIter) -> getX( &(right[0] ) );

            // compute new level set
            base::cut::analyticLevelSet( mesh,
                                         boost::bind( &interval<dim>, _1, _2,
                                                      left[0], right[0] ),
                                         true, levelSet );

            // update the cut cell structure
            base::cut::generateCutCells( mesh, levelSet, cells );
        }

        // 3) advect data
        {
            // compute new supports in order to identify the DoFs
            std::vector<double> newSupports;
            base::cut::supportComputation( mesh, displacement,
                                           cutQuadrature, newSupports );

            // Find the location of the DoFs in the previous configuration
            std::vector<std::pair<std::size_t,VecDim> > previousDoFLocation;
            findPreviousDoFLocations( mesh, displacement, newSupports,
                                      doFLocation, previousDoFLocation,
                                      supportThreshold, 1.e-6, 10 );

            // advect displacement field from previous to new location
            advectField( mesh, displacement, previousDoFLocation,
                         newSupports, supportThreshold );
            // advect velocity field from previous to new location
            advectField( mesh, velocity,     previousDoFLocation,
                         newSupports, supportThreshold );

        }

        // push history
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

        // push history
        std::for_each( velocity.doFsBegin(), velocity.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

        std::cout << step * stepSize << "  " << fFun << "  ";
        std::cout << 0.5 * (left[0] + right[0]) - 0.5 * (a+b) << "  "
                  << f * (1. - std::cos( step * stepSize ) )
                  << std::endl;

    }

    return 0;
}
