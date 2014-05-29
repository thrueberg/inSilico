//#define WRITEMAT


#include <string>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <base/dof/location.hpp>

#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/Field.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/mesh/Size.hpp>

#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>
#include <base/io/List.hpp>

#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/ScaledField.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>
#include <base/cut/stabiliseBasis.hpp>
#include <base/cut/evaluateOnCutCells.hpp>

#include <base/fe/Basis.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/location.hpp>
#include <base/dof/constrainBoundary.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <heat/Laplace.hpp>

#include <solid/HyperElastic.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>
#include <base/post/ErrorNorm.hpp>

#include "generateMesh.hpp"

const double coordTol = 1.e-6;

//------------------------------------------------------------------------------
//! Helper that converts the status of a DoF into a number
template<typename DOF>
struct DoFStatus
{
    static double apply( const DOF* doFPtr, const unsigned component )
    {
        if ( doFPtr -> isConstrained( component ) ) return 0.;
        if ( doFPtr -> isActive(      component ) ) return 1.;
        return -1.;
    }
};

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM>::Type&x, DOF* doFPtr )
{
    // fix the box face with x_1 = 0
    if ( base::auxi::almostEqualNumbers( x[0], 0. ) ) {

        for ( unsigned d = 0; d < DOF::size; d++ ) {
            if ( doFPtr -> isActive(d) ) {
                doFPtr -> constrainValue( d, 0. );
            }
        }
    }
    return;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
base::Vector<1>::Type
solution( const typename base::Vector<DIM>::Type& x,
          const double kappa, const double L )
{
    const double alpha = 3. * M_PI / 2. / L;
    return base::constantVector<1>( std::sin( alpha * x[0] ) );
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Matrix<DIM,1>::Type
gradient( const typename base::Vector<DIM>::Type& x,
          const double kappa, const double L )
{
    const double alpha = 3. * M_PI / 2. / L;
    typename base::Matrix<DIM,1>::Type result =
        base::constantMatrix<DIM,1>( 0. );

    result(0,0) = alpha * std::cos( alpha * x[0] );

    return result;
}



//------------------------------------------------------------------------------
template<unsigned DIM>
base::Vector<1>::Type
forceFun( const typename base::Vector<DIM>::Type& x,
          const double kappa, const double L )
{
    const double alpha = 3. * M_PI / 2. / L;
    return base::constantVector<1>( kappa * alpha * alpha * std::sin( alpha * x[0] ) );
}

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
void writeVTKFile( const std::string& baseName,
                   const MESH&    mesh,
                   const FIELD&   field,
                   const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet )
{
    const std::string vtkFile = baseName + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );

    std::vector<double> distances;
    std::transform( levelSet.begin(), levelSet.end(),
                    std::back_inserter( distances ),
                    boost::bind( &base::cut::LevelSet<MESH::Node::dim>::getSignedDistance, _1 ) );

    vtkWriter.writeUnstructuredGrid( mesh );
    vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
    base::io::vtk::writePointData( vtkWriter, mesh, field, "u   " );

#if 1
    //--------------------------------------------------------------------------
    // write the status of every DoF
    {
        static const unsigned size = FIELD::DegreeOfFreedom::size;
        std::vector<typename base::Vector<size>::Type> doFStatus;
        
        base::post::evaluateAtNodes(
            mesh, field,
            boost::bind( DoFStatus<typename FIELD::DegreeOfFreedom>::apply, _1, _2 ),
            doFStatus );

        vtkWriter.writePointData( doFStatus.begin(), doFStatus.end(), "status" );
    }
#endif
    vtk.close();
}

//------------------------------------------------------------------------------
// Signed distance to a sphere, outside positive
template<unsigned DIM>
bool plane( const typename base::Vector<DIM>::Type& x,
            typename base::Vector<DIM>::Type& xClosest,
            const double L )
{
    xClosest = x;
    xClosest[0] = L;

    if ( x[0] > L ) return false;
    return true;
}


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // spatial dimension
    const unsigned    dim      = SPACEDIM;
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 3;
    const unsigned    kernelDegEstimate = 5;


    const bool isLaplace = true;
    
    if ( argc != 3 ) {
        std::cerr << "Usage: " << argv[0] << " N  input.dat\n"
                  << "(Compiled for dim=" << dim << ")\n\n";
        return -1;
    }

    // read name of input file
    const unsigned    numElements = boost::lexical_cast<unsigned>(    argv[1] );
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double xmax, E, nu, f, tolerance, cutThreshold, L;
    unsigned numLoadSteps, maxIter;
    std::string meshFile;
    unsigned stabilise; // 0 - not, 1 - Höllig, 2 - scaling
    bool compute;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "xmax",         xmax );
        prop.registerPropertiesVar( "E",            E );
        prop.registerPropertiesVar( "nu",           nu );
        prop.registerPropertiesVar( "f",            f );
        prop.registerPropertiesVar( "numLoadSteps", numLoadSteps );
        prop.registerPropertiesVar( "maxIter",      maxIter );
        prop.registerPropertiesVar( "tolerance",    tolerance );
        prop.registerPropertiesVar( "L",            L );
        
        prop.registerPropertiesVar( "stabilise",    stabilise );
        prop.registerPropertiesVar( "meshFile",     meshFile );
        prop.registerPropertiesVar( "compute",      compute );
        prop.registerPropertiesVar( "cutThreshold", cutThreshold );
        
        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );

    }

    // in case of no computation, do not stabilise
    if ( not compute ) stabilise = 0;
    
    // basic attributes of the computation
    const base::Shape             shape = //base::SimplexShape<dim>::value;
        base::HyperCubeShape<dim>::value;
#if SPACEDIM > 1
    const base::Shape surfShape = base::SimplexShape<dim-1>::value;
#else
    const base::Shape surfShape = base::POINT;
#endif

    const unsigned              doFSize = (isLaplace ? 1 : dim);
    const unsigned              nHist   = 1;

    // Bulk mesh
    typedef base::Unstructured<shape,geomDeg>  Mesh;
    Mesh mesh;
    std::string baseName;
    if ( numElements > 0 ) {
        const double left = 0;
        generateMesh( mesh, numElements, left, xmax );

        baseName = "oneDimensional." + base::io::leadingZeros( numElements );
    }
    else{
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );

        baseName = base::io::baseName( meshFile, ".smf" );
    }

    // Boundary mesh
    typedef base::mesh::BoundaryMeshBinder<Mesh,true>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );
    base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                      meshBoundary.end(),
                                      mesh, boundaryMesh );
        
    typedef Mesh::Node::VecDim VecDim;

    // Cell structures
    typedef base::cut::Cell<shape,geomDeg> Cell;
    std::vector<Cell> cells;

#if SPACEDIM > 1    
    typedef base::cut::Cell<surfShape,geomDeg> SurfCell;
    std::vector<SurfCell> surfCells;
#endif

    // level set
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::analyticLevelSet( mesh,
                                 boost::bind( &plane<dim>, _1, _2, L ), 
                                 true, levelSet );

    //----------------------------------------------------------------------
    base::cut::generateCutCells( mesh, levelSet, cells );
    
    // Make cut cell structure (surface)
#if SPACEDIM > 1
    base::cut::generateCutCells( boundaryMesh, levelSet, surfCells );
#endif


    typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
    SurfaceMesh surfaceMesh;
    base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, surfaceMesh );

    //--------------------------------------------------------------------------
    // FE
    typedef base::fe::Basis<shape,fieldDeg>               FEBasis;
    typedef base::cut::ScaledField<FEBasis,doFSize,nHist> Field;
    Field field;
    base::dof::generate<FEBasis>( mesh, field );

    // for domain field
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type UU;

    // for surface field
    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder boundaryFieldBinder( boundaryMesh, field );
    SurfaceFieldBinder  surfaceFieldBinder(  surfaceMesh, field );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SU;

    //--------------------------------------------------------------------------
    // Quadratures
    typedef base::cut::Quadrature<kernelDegEstimate,shape,Cell> CutQuadrature;
    CutQuadrature cutQuadrature( cells, true );

    typedef base::Quadrature<kernelDegEstimate,surfShape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

#if SPACEDIM > 1
    typedef base::cut::Quadrature<kernelDegEstimate,surfShape,SurfCell> SurfaceCutQuadrature;
    SurfaceCutQuadrature surfaceCutQuadrature( surfCells, true );
#endif

    // compute supports, scale basis
    const std::size_t numDoFs = std::distance( field.doFsBegin(),
                                               field.doFsEnd() );
    std::vector<double> supports;
    supports.resize(  numDoFs );
    
    base::cut::supportComputation( mesh, field, cutQuadrature,  supports );

    if ( stabilise==1 ) {
        std::vector<std::pair<std::size_t,VecDim> > doFLocation;
        base::dof::associateLocation( field, doFLocation );
        base::cut::stabiliseBasis( mesh, field, supports, doFLocation );
    }
    else if ( stabilise==2 ) {
        field.scaleAndTagBasis(  supports,  1.e-10 );
    }
    else field.tagBasis( supports,  1.e-10 );

    typedef base::Vector<1>::Type VecDoF;

    // apply dirichlet constraints
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, field, 
                                           boost::bind( &dirichletBC<dim,
                                                        Field::DegreeOfFreedom>,
                                                        _1, _2 ) );
    // number DoFs
    const std::size_t activeDoFsU = 
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );

    typedef mat::hypel::NeoHookeanCompressible Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );
    typedef solid::HyperElastic<Material,UU::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    typedef heat::Laplace<UU::Tuple> Laplace;
    Laplace laplace( E );

    //unsigned numCGIter;
    if ( compute ) {

        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( activeDoFsU );

        if ( isLaplace ) {
            base::asmb::stiffnessMatrixComputation<UU>( cutQuadrature, solver,
                                                        fieldBinder, laplace, 
                                                        false ); //iter > 0 );
        }
        else {

            base::asmb::computeResidualForces<UU>( cutQuadrature, solver,
                                                   fieldBinder, hyperElastic );
            
            base::asmb::stiffnessMatrixComputation<UU>( cutQuadrature, solver,
                                                        fieldBinder, hyperElastic,
                                                        false ); //iter > 0 );

        }

        
        base::asmb::bodyForceComputation<UU>( cutQuadrature, solver,
                                              fieldBinder,
                                              boost::bind( &forceFun<dim>,
                                                           _1, E, L ) );

        // Finalise assembly
        solver.finishAssembly();

#ifdef WRITEMAT
        {
            std::string matName = baseName + ".mat";
            std::ofstream mat( matName.c_str() );
            solver.debugLHS( mat );

            std::string vecName = baseName + ".vec";
            std::ofstream vec( vecName.c_str() );
            solver.debugRHS( vec );
        }
#endif
        
        // Solve
        solver.superLUSolve();
        //numCGIter = solver.cgSolve();

        
        // distribute results back to dofs
        base::dof::addToDoFsFromSolver( solver, field );

    }

    
    {
        std::vector<double> cutCellDistance;
        base::cut::evaluateCutCellNodeDistances( mesh, levelSet, cells,
                                                 cutCellDistance );

        std::vector<VecDoF> cutCellSolution;
        base::cut::evaluateFieldAtCutCellNodes( mesh, field, cells,
                                                cutCellSolution );

        typedef base::cut::SurfaceSimplexMesh<Mesh,geomDeg> SurfaceSimplexMesh;
        SurfaceSimplexMesh surfaceSimplexMesh;
        base::cut::extractSurfaceMeshFromCutCells( mesh, cells, surfaceSimplexMesh );

        const std::string vtkFileSurf = baseName + ".surf.vtk";
        std::ofstream vtkSurf( vtkFileSurf.c_str() );
        base::io::vtk::LegacyWriter vtkWriterSurf( vtkSurf );

        vtkWriterSurf.writeUnstructuredGrid( surfaceSimplexMesh );
        vtkWriterSurf.writePointData( cutCellDistance.begin(), cutCellDistance.end(), "distances" );
        vtkWriterSurf.writePointData( cutCellSolution.begin(), cutCellSolution.end(), "u" );

    }

    //--------------------------------------------------------------------------
    // write a vtk file
    writeVTKFile( baseName, mesh, field, levelSet );

    double hmin, hmax;
    {
        std::vector<double> h;
        Mesh::ElementPtrConstIter eIter = mesh.elementsBegin();
        Mesh::ElementPtrConstIter eEnd  = mesh.elementsEnd();
        for ( ; eIter != eEnd; ++eIter )  {
            h.push_back( base::mesh::elementSize( *eIter ) );
        }
     
        hmin = *std::min_element( h.begin(), h.end() );
        hmax = *std::max_element( h.begin(), h.end() );
    }

    std::cout << hmin << "  " << hmax << "  ";
    //std::cout << numElements << "  ";
    
    //--------------------------------------------------------------------------
    // compute L2-error
    std::cout << base::post::errorComputation<0>(
        cutQuadrature, mesh, field,
        boost::bind( &solution<dim>, _1, E, L  ) )
        << "  ";

    //--------------------------------------------------------------------------
    // compute H1-error
    std::cout << base::post::errorComputation<1>(
        cutQuadrature, mesh, field,
        boost::bind( &gradient<dim>, _1, E, L ) )
              << "  ";
    
#if 1
    //--------------------------------------------------------------------------
    // Compute mesh volume
    double volume = 0.;
    base::asmb::simplyIntegrate<UU>( cutQuadrature, volume, fieldBinder,
                                     base::kernel::Measure<UU::Tuple>() );
    std::cout << volume << "  ";

    //--------------------------------------------------------------------------
    // Compute surface area
    double area = 0.;
    base::asmb::simplyIntegrate<SU>( surfaceQuadrature, area, surfaceFieldBinder,
                                     base::kernel::Measure<SU::Tuple>() );

#if SPACEDIM > 1    
    base::asmb::simplyIntegrate<SU>( surfaceCutQuadrature, area, boundaryFieldBinder,
                                     base::kernel::Measure<SU::Tuple>() );
#endif


    std::cout << area << "  ";
    
#endif

    //std::cout << numCGIter;
    std::cout << std::endl;

    return 0;
}
