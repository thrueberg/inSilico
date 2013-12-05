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
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>

//#include <solid/HyperElastic.hpp>
#include "HyperElastic.hpp"
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include "generateMesh.hpp"
#include "findLocation.hpp"

#include "extrapolateToOutside.hpp"
#include "advectField.hpp"
#include "moveSurface.hpp"

const double coordTol = 1.e-6;

//------------------------------------------------------------------------------
// Level set function for a domain covering the interval  0 <= x_1 <= L
template<unsigned DIM>
bool interval( const typename base::Vector<DIM>::Type& x,
                     typename base::Vector<DIM>::Type& xClosest,
               const double L )
{
    xClosest[0] = L;
    if ( x[0] <= L ) return true;

    return false;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM>::Type
tractionFun( const typename base::Vector<DIM>::Type& x,
             const typename base::Vector<DIM>::Type& normal,
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
    double xmax, L, E, nu, f, tolerance;
    unsigned numLoadSteps, maxIter;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "xmax",         xmax );
        prop.registerPropertiesVar( "L",            L );
        prop.registerPropertiesVar( "E",            E );
        prop.registerPropertiesVar( "nu",           nu );
        prop.registerPropertiesVar( "f",            f );
        prop.registerPropertiesVar( "numLoadSteps", numLoadSteps );
        prop.registerPropertiesVar( "maxIter",      maxIter );
        prop.registerPropertiesVar( "tolerance",    tolerance );

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
    const unsigned                nHist = 1;

    typedef base::Unstructured<shape,geomDeg>  Mesh;
    typedef Mesh::Node::VecDim VecDim;
    Mesh mesh;
    {
        generateMesh( mesh, numElements, 0., xmax );
    }

    //--------------------------------------------------------------------------
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::analyticLevelSet( mesh,
                                 boost::bind( &interval<dim>, _1, _2, L ),
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
    Field displacement;
    base::dof::generate<FEBasis>( mesh, displacement );

    // find geometry association for the dofs
    std::vector<std::pair<std::size_t,VecDim> > doFLocation;
    base::dof::associateLocation( displacement, doFLocation );

    // for domain field
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, displacement );
    typedef FieldBinder::TupleBinder<1,1>::Type UU;

    typedef mat::hypel::NeoHookeanCompressible Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );
    typedef HyperElastic<Material,UU::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    writeVTKFile( "euler", 0, mesh, displacement, levelSet );

    VecDim right;

    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    
    //--------------------------------------------------------------------------
    // * Loop over load steps *
    //--------------------------------------------------------------------------
    const double supportThreshold = 1.e-10;
    for ( unsigned step = 0; step < numLoadSteps; step++ ) {

        // for surface field
        typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
        SurfaceMesh immersedSurface;
        {
            base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, immersedSurface );
        }
        
        typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Field> SurfaceFieldBinder;
        SurfaceFieldBinder surfaceFieldBinder( immersedSurface, displacement );
        typedef SurfaceFieldBinder::TupleBinder<1>::Type SU;

        //----------------------------------------------------------------------
        // 0) Solve Lagrangian problem
        displacement.activateAll();

        // fix left end
        Field::DoFPtrIter dIter = displacement.doFsBegin();
        for ( unsigned d = 0; d < dim; d++ ) (*dIter) -> constrainValue( d, 0. );
    
        typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
        CutQuadrature cutQuadrature( cells, true );

        // compute supports, scale basis
        std::vector<double> supports;
        base::cut::supportComputation( mesh, displacement, cutQuadrature, supports );

        displacement.scaleAndTagBasis( supports,  supportThreshold );
    
        // number DoFs 
        const std::size_t activeDoFsU = 
            base::dof::numberDoFsConsecutively( displacement.doFsBegin(),
                                                displacement.doFsEnd() );
        const double fFun = f *
            static_cast<double>(step+1) / static_cast<double>( numLoadSteps );

        //----------------------------------------------------------------------
        // Nonlinear iterations
        for ( unsigned iter = 0; iter < maxIter; iter++ ) {
    
            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( activeDoFsU  );

            base::asmb::computeResidualForces<UU>( cutQuadrature, solver,
                                                   fieldBinder, hyperElastic );
            
            base::asmb::stiffnessMatrixComputation<UU>( cutQuadrature, solver,
                                                        fieldBinder, hyperElastic,
                                                        iter > 0 );

            base::asmb::neumannForceComputation<SU>( surfaceQuadrature,
                                                     solver, surfaceFieldBinder,
                                                     boost::bind( &tractionFun<dim>,
                                                                  _1, _2, fFun ) );


            // Finalise assembly
            solver.finishAssembly();

            // norm of residual
            const double conv1 = solver.norm() / E;

            std::cout << "* " << iter << " " << conv1 << "  ";

            if ( conv1 < tolerance ) {
                std::cout << std::endl;
                break;
            }

            // Solve
            solver.superLUSolve();
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, displacement, iter > 0 );

            const double conv2 = solver.norm();
            std::cout << conv2 << std::endl;

            if ( conv2 < tolerance ) break;
        }


        displacement.unscaleValues( supports, supportThreshold );

        // write a vtk file
        writeVTKFile( "euler", step+1, mesh, displacement, /*velocity,*/ levelSet );


        // 1) Pass Data to complementary domain
        {
            const double h = xmax / static_cast<double>(numElements); // mesh size
            const double extL = h; // extrapolation length

            // pass extrapolated solution to inactive DoFs
            extrapolateToOutside( mesh, displacement, extL, doFLocation, levelSet );
        }
        
        // 2) Geometry update
        {
            // Move with displacement solution
            moveSurface( mesh, displacement, immersedSurface );

            // get location of domain boundary
            SurfaceMesh::NodePtrIter surfNIter = immersedSurface.nodesBegin();
            (*surfNIter) -> getX( &(right[0] ) );

            // compute new level set
            base::cut::analyticLevelSet( mesh,
                                         boost::bind( &interval<dim>, _1, _2,
                                                      right[0] ),
                                         true, levelSet );

            // update the cut cell structure
            base::cut::generateCutCells( mesh, levelSet, cells );
        }

        // write a vtk file
        writeVTKFile( "euler2", step+1, mesh, displacement, /*velocity,*/ levelSet );

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

        }
        
        // write a vtk file
        writeVTKFile( "euler3", step+1, mesh, displacement, /*velocity,*/ levelSet );

        // push history
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

        std::cout << step << "  " << fFun << " ";
        std::cout << right[0] - L << "  " << std::endl;

    }

    return 0;
}
