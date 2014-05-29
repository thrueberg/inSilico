#include <string>
#include <boost/lexical_cast.hpp>

#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/Field.hpp>
#include <base/mesh/Size.hpp>

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

#include <base/cut/stabiliseBasis.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>

#include <base/fe/Basis.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/location.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/SimpleIntegrator.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>

#include "HyperElastic.hpp"
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include "generateMesh.hpp"
#include "findLocation.hpp"

#include "extrapolateToOutside.hpp"
#include "advectField.hpp"
#include "moveSurface.hpp"

//-----------------------------------------------------------------------------
// Applied surface traction
template<unsigned DIM>
typename base::Vector<DIM>::Type
tractionFun( const typename base::Vector<DIM>::Type& x,
             const typename base::Vector<DIM>::Type& normal,
             const double value )
{
    const double phi = std::atan2( x[1], x[0] );
    
    typename base::Vector<DIM>::Type f = value * normal
        * (std::cos(2.*phi));
    return f;
}
//-----------------------------------------------------------------------------
#include "writeVTK.hpp"

const double coordTol = 1.e-4;

//------------------------------------------------------------------------------
// Analytic level set function for a spherical domain
template<unsigned DIM>
bool spherical( const typename base::Vector<DIM>::Type& x,
                typename base::Vector<DIM>::Type& xClosest,
                const double L )
{
    //xClosest[0] = L;
    if ( x.norm() < coordTol ) {
        xClosest = base::constantVector<DIM>(0.);
        xClosest[0] = L;
    }
    else{
        xClosest = (L / x.norm()) * x;
    }

    if ( x.norm() <= L ) return true;
    return false;
}



//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // spatial dimension
    const unsigned    dim = 2; 
    
    if ( argc != 3 ) {
        std::cerr << "Usage: " << argv[0] << " N  input.dat\n"
                  << "(Compiled for dim=" << dim << ")\n\n";
        return -1;
    }
    
    // read name of input file
    const unsigned    numElements = boost::lexical_cast<unsigned>(    argv[1] );
    //const std::string   smfFile   = boost::lexical_cast<std::string>( argv[1] );
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

    // basic attributes of the computation
    const unsigned             geomDeg  = 1;
    const unsigned             fieldDeg = 1;
    const base::Shape             shape = base::HyperCubeShape<dim>::value;
    const unsigned    kernelDegEstimate = 4;
    const unsigned              doFSize = dim;
    const unsigned                nHist = 1;

    typedef base::Unstructured<shape,geomDeg>  Mesh;
    typedef Mesh::Node::VecDim VecDim;
    Mesh mesh;
    {
        const double left = (dim == 1 ? 0. : -xmax );
        generateMesh( mesh, numElements, left, xmax );
        //std::ifstream smf( smfFile.c_str() );
        //base::io::smf::readMesh( smf, mesh );
    }

    const double h = base::mesh::Size<Mesh::Element>::apply( mesh.elementPtr(0) );

    //--------------------------------------------------------------------------
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::analyticLevelSet( mesh,
                                 boost::bind( &spherical<dim>, _1, _2, L ),
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


    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    double initialSurf;
    double surf;
    
    //--------------------------------------------------------------------------
    // * Loop over load steps *
    //--------------------------------------------------------------------------
    
    const double supportThreshold = std::numeric_limits<double>::min();
    for ( unsigned step = 0; step < numLoadSteps; step++ ) {

        // for surface field
        typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
        SurfaceMesh immersedSurface;
        base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, immersedSurface );
        
        typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Field> SurfaceFieldBinder;
        SurfaceFieldBinder surfaceFieldBinder( immersedSurface, displacement );
        typedef SurfaceFieldBinder::TupleBinder<1>::Type SU;

        // compute size of surface
        surf = 0.;
        base::asmb::simplyIntegrate<SU>( surfaceQuadrature, surf, surfaceFieldBinder,
                                         base::kernel::Measure<SU::Tuple>() );
        if ( step == 0 ) initialSurf = surf;
            
        //----------------------------------------------------------------------
        // 0) Solve Lagrangian problem
        displacement.activateAll();

        typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
        CutQuadrature cutQuadrature( cells, true );

        // compute supports, scale basis
        std::vector<double> supports;
        base::cut::supportComputation( mesh, displacement, cutQuadrature, supports );

        //displacement.scaleAndTagBasis( supports,  supportThreshold );
        //displacement.tagBasis( supports,  supportThreshold );
        base::cut::stabiliseBasis( mesh, displacement, supports, doFLocation );

        // fix some points
        {
            VecDim point = base::constantVector<dim>( 0. );
            const std::size_t index = base::dof::findDoFWithLocation( doFLocation, mesh,
                                                                      point, coordTol );
            

            Field::DegreeOfFreedom* dPtr = displacement.doFPtr( index );
            for ( unsigned d = 0; d < dim; d++ ) dPtr -> constrainValue( d, 0. );

            if ( dim > 1 ) {
                // Fix two more points along the x_1 axis
                const unsigned bla = static_cast<unsigned>(L/h) - 2;
                point[0] = bla * h;
                const std::size_t index1 = base::dof::findDoFWithLocation( doFLocation, mesh,
                                                                           point, coordTol );
                
                Field::DegreeOfFreedom* dPtr1 = displacement.doFPtr( index1 );
                for ( unsigned d = 1; d < dim; d++ ) dPtr1 -> constrainValue( d, 0. );

                point[0] = -static_cast<double>(bla) * h;
                const std::size_t index2 = base::dof::findDoFWithLocation( doFLocation, mesh,
                                                                           point, coordTol );
                
                Field::DegreeOfFreedom* dPtr2 = displacement.doFPtr( index2 );
                for ( unsigned d = 1; d < dim; d++ ) dPtr2 -> constrainValue( d, 0. );
            }

        }

        // number DoFs 
        const std::size_t activeDoFsU = 
            base::dof::numberDoFsConsecutively( displacement.doFsBegin(),
                                                displacement.doFsEnd() );

        std::cout << "Num DoFs: " << activeDoFsU << std::endl;

        // account for change in surface area
        const double fFun = f * (initialSurf / surf) *
            static_cast<double>(step+1) / static_cast<double>( numLoadSteps );

        // start from 0 for the increments
        base::dof::clearDoFs( displacement );

        // write a vtk file
        writeVTKFile( "euler", (step+1)*100, mesh, displacement, levelSet, supports );
        writeSurfaceVTKFile( "eulerSurf", (step+1)*100, immersedSurface, fFun );

        //----------------------------------------------------------------------
        // Nonlinear iterations
        unsigned iter = 0;
        for ( ; iter < maxIter; iter++ ) {
    
            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( activeDoFsU  );

            base::asmb::computeResidualForces<UU>( cutQuadrature, solver,
                                                   fieldBinder, hyperElastic );

            // (WW) constraint application not consistent yet
            base::asmb::stiffnessMatrixComputation<UU>( cutQuadrature, solver,
                                                        fieldBinder, hyperElastic,
                                                        iter==0 );

            base::asmb::neumannForceComputation<SU>( surfaceQuadrature,
                                                     solver, surfaceFieldBinder,
                                                     boost::bind( &tractionFun<dim>,
                                                                  _1, _2, fFun ) );

            // Finalise assembly
            solver.finishAssembly();

            // norm of residual
            const double conv1 = solver.norm() / E;

            std::cout << "* " << iter << " " << conv1 << "  " << std::flush;

            if ( isnan( conv1 ) ) return 1;
            
            if ( conv1 < tolerance ) {
                std::cout << std::endl;
                break;
            }

            // Solve
            solver.superLUSolve();
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, displacement );

            writeVTKFile( "euler", (step+1)*100+iter+1, mesh, displacement, levelSet, supports );
            writeSurfaceVTKFile( "eulerSurf", (step+1)*100+iter+1, immersedSurface, fFun );

            const double conv2 = solver.norm();
            std::cout << conv2 << std::endl;

            if ( conv2 < tolerance ) break;
        }

        // add previous solution to current solution: u_{n+1} = \Delta u + u_n
        {
            Field::DoFPtrIter dIter = displacement.doFsBegin();
            Field::DoFPtrIter dEnd  = displacement.doFsEnd();
            for ( ; dIter != dEnd; ++dIter ) {
                for ( unsigned d = 0; d < dim; d++ ) {
                    const double prevU  = (*dIter) -> getHistoryValue<1>( d );
                    const double deltaU = (*dIter) -> getHistoryValue<0>( d );
                    const double currU = prevU + deltaU;
                    (*dIter) -> setValue( d, currU );
                }
            }
        }

        //displacement.unscaleValues( supports, supportThreshold );

        // write a vtk file
        writeVTKFile( "euler", (step+1)*100+iter+2, mesh, displacement, levelSet, supports );
        writeSurfaceVTKFile( "eulerSurf", (step+1)*100+iter+2, immersedSurface, fFun );


        // 1) Pass Data to complementary domain
        {
            std::cout << "* Extrapolate" << std::endl;
            const double extL = h; // extrapolation length

            // pass extrapolated solution to inactive DoFs
            extrapolateToFictitious( mesh, displacement, extL, doFLocation, levelSet );
        }
        
        // write a vtk file
        writeVTKFile( "euler", (step+1)*100+iter+3, mesh, displacement, levelSet, supports );
        writeSurfaceVTKFile( "eulerSurf", (step+1)*100+iter+3, immersedSurface, fFun );

        // 2) Geometry update
        {
            std::cout << "* New geometry" << std::endl;
            
            // Move with displacement solution
            moveSurface( mesh, displacement, immersedSurface );


            base::cut::bruteForce( mesh, immersedSurface, true, levelSet );

            // update the cut cell structure
            base::cut::generateCutCells( mesh, levelSet, cells );

        }


        // 3) advect data
        {
            std::cout << "* Advection" << std::endl;
            
            //std::vector<double> newSupports;
            base::cut::supportComputation( mesh, displacement, cutQuadrature, supports ); 

            // Find the location of the DoFs in the previous configuration
            std::vector<std::pair<std::size_t,VecDim> > previousDoFLocation;
            findPreviousDoFLocations( mesh, displacement, supports, 
                                      doFLocation, previousDoFLocation,
                                      supportThreshold, 1.e-6, 10 );

            // advect displacement field from previous to new location
            advectField( mesh, displacement, previousDoFLocation, supports, 
                         supportThreshold );

        }
        
        // write a vtk file
        writeVTKFile( "euler", (step+1)*100+iter+4, mesh, displacement, levelSet, supports );
        writeSurfaceVTKFile( "eulerSurf", (step+1)*100+iter+4, immersedSurface, fFun );

        // push history
        base::dof::pushHistory( displacement );

        std::cout << step << "  " << fFun << " ";
        std::cout << surf  << "  " << std::endl;

    } // end load-steps

    return 0;
}
