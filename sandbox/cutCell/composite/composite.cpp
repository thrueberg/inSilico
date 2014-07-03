//#define LINEAR

#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/Unstructured.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/io/Format.hpp>
#include <base/io/PropertiesParser.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>

#include <base/Quadrature.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/solver/Eigen3.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>
#include <solid/HyperElastic.hpp>
#include <solid/Stress.hpp>

#include <base/cut/LevelSet.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/stabiliseBasis.hpp>

#include <base/post/ErrorNorm.hpp>

#include <base/asmb/SimpleIntegrator.hpp>

#include <base/auxi/BoundingBox.hpp>
#include <base/cut/ImplicitFunctions.hpp>

#include "../../generateMesh.hpp"
#include "../HyperElastic.hpp"
#include "../InterfaceField.hpp"
#include "../SurfaceField.hpp"
#include "../ImplicitGeometry.hpp"

//------------------------------------------------------------------------------
const double coordTol = 1.e-6;


//------------------------------------------------------------------------------
// Dirichlet boundary conditions 
template<unsigned DIM, typename DOF>
void fix( const typename base::Vector<DIM,double>::Type& x, DOF* doFPtr,
          const base::auxi::BoundingBox<DIM>& bbox )
{
    const bool isBot = bbox.isOnLowerBoundary( x, DIM-1, coordTol );
    const bool isTop = bbox.isOnUpperBoundary( x, DIM-1, coordTol );

    if ( (not isBot) and (not isTop) )  return;

    if ( isBot ) {
        for ( unsigned d = 0; d < DIM; d++ ) {
            if ( doFPtr -> isActive(d) ) {
                doFPtr -> constrainValue( d, 0.0 );
            }
        }
    }

    if ( isTop ) {
        if ( doFPtr -> isActive(DIM-1) ) {
            doFPtr -> constrainValue( DIM-1, 0.0 );
        }
    
    }
    
    return;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM,double>::Type
shearTop( const typename base::Vector<DIM,double>::Type& x,
          const typename base::Vector<DIM,double>::Type& normal,
          const base::auxi::BoundingBox<DIM>& bbox,
          const double value )
{
    const bool isTop = bbox.isOnUpperBoundary( x, DIM-1, coordTol );
    
    if ( not isTop)  return base::constantVector<DIM>( 0. );

    typename base::Vector<DIM,double>::Type t = base::constantVector<DIM>( 0. );
    t[0] = value;
    
    
    return t;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM,double>::Type
sideLoad( const typename base::Vector<DIM,double>::Type& x,
          const typename base::Vector<DIM,double>::Type& normal,
          const base::auxi::BoundingBox<DIM>& bbox,
          const double value )
{
    const bool isSide =
        ( bbox.isOnLowerBoundary( x, 0, coordTol ) or
          bbox.isOnUpperBoundary( x, 0, coordTol ) );
    
    if ( not isSide)  return base::constantVector<DIM>( 0. );

    typename base::Vector<DIM,double>::Type t = base::constantVector<DIM>( 0. );
    t[0] = value;
    
    return t;
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    Eigen::initParallel();
    
    //--------------------------------------------------------------------------
    // User input
    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0]
                  << " file.smf N input.dat \n\n";
        return -1;
    }

    const unsigned dim = 3;
    typedef base::Vector<dim>::Type VecDim;


        // read name of input file
    const unsigned    numElements = boost::lexical_cast<unsigned>(    argv[1] );
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double E1, nu1, E2, nu2, thickness;
    double forceVal;
    unsigned numFibres;
    double   alpha, gamma;
    double   tolerance;
    unsigned maxIter, numLoadSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "E1",           E1 );
        prop.registerPropertiesVar( "nu1",          nu1 );
        prop.registerPropertiesVar( "E2",           E2 );
        prop.registerPropertiesVar( "nu2",          nu2 );
        prop.registerPropertiesVar( "gamma",        gamma );
        
        prop.registerPropertiesVar( "thickness",    thickness );
        
        prop.registerPropertiesVar( "forceVal",     forceVal );

        prop.registerPropertiesVar( "numFibres",    numFibres );
        prop.registerPropertiesVar( "alpha",        alpha );

        prop.registerPropertiesVar( "tolerance",    tolerance );
        prop.registerPropertiesVar( "maxIter",      maxIter );
        prop.registerPropertiesVar( "numLoadSteps", numLoadSteps );
        
        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );

#ifdef LINEAR
        numLoadSteps = 1;
        maxIter      = 1;
#endif
        
    }

    // basic attributes of the computation
    const unsigned             geomDeg  = 1;
    const unsigned             fieldDeg = 1;
    const base::Shape             shape = base::HyperCubeShape<dim>::value;
    const base::Shape         surfShape = base::SimplexShape<dim-1>::value;
    const unsigned    kernelDegEstimate = 5;

    // Bulk mesh
    typedef base::Unstructured<shape,geomDeg>  Mesh;
    typedef Mesh::Node::VecDim VecDim;
    Mesh mesh;
    const std::string baseName =
        std::string("composite.") + base::io::leadingZeros( numElements );

    VecDim a, b;
    a[0] = a[2] = -.5;
    a[1] = -thickness/2.; b[1] = thickness/2.;
    b[0] = b[2] =  .5;

    base::auxi::BoundingBox<dim> bbox( a, b );
    
    base::Vector<dim,unsigned>::Type N;
    N[0] = N[2] = numElements;
    N[1] = static_cast<unsigned>( numElements * thickness );
    generateMesh<dim>( mesh, N, a, b );

    // Boundary mesh
    typedef base::mesh::BoundaryMeshBinder<Mesh,true>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );
    base::mesh::generateBoundaryMesh( meshBoundary.begin(), meshBoundary.end(),
                                      mesh, boundaryMesh );

    
    //--------------------------------------------------------------------------
    // go through immersed surfaces
    ImplicitGeometry<Mesh> geometry( mesh, boundaryMesh );

    // cylinder thickness
    const double rc = thickness / 6.;
    // cylinder separation
    const double sc = 1. / static_cast<double>( numFibres );
    // cylinder axis
    VecDim ca = base::constantVector<dim>( 0. );
    ca[0] = std::sin( alpha / 180. * M_PI );
    ca[2] = std::cos( alpha / 180. * M_PI );

    for ( unsigned c = 0; c < 2 * numFibres; c ++ ) {
    
        // cylinder axis point
        VecDim cp = base::constantVector<dim>( 0. );
        cp[0] = -1. + c * sc;
    
        base::cut::Cylinder<dim> cyl( rc, cp, ca, false );
    
        geometry.intersectAnalytical( cyl, 0. );
    }

    //VecDim cp = base::constantVector<dim>( 0. );
    //ca[0] = 1.0; ca[1] = 0.; ca[2] = 0.;
    //base::cut::Cylinder<dim> cyl( 2.*rc, cp, ca, false );
    //geometry.intersectAnalytical( cyl, 1.e-6 );
    // Plane
    //base::cut::Plane<dim> pl( ca, 0.033333333, false );
    //geometry.intersectAnalytical( pl, 0. );

    // Cell structures
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells = geometry.getCells();

    typedef base::cut::Cell<surfShape> SurfCell;
    std::vector<SurfCell> surfCells = geometry.getSurfCells();

    // intersection of all level sets
    typedef base::cut::LevelSet<dim> LevelSet;
    const std::vector<LevelSet> levelSet = geometry.getLevelSet();

    // Generate a mesh from the immersed surface
    typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
    SurfaceMesh surfaceMesh;
    base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, surfaceMesh );

    //--------------------------------------------------------------------------
    // Quadratures
    typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
    CutQuadrature cutQuadrature1( cells, true  );
    CutQuadrature cutQuadrature2( cells, false );

    typedef base::Quadrature<kernelDegEstimate,surfShape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    typedef base::cut::Quadrature<kernelDegEstimate,surfShape> SurfaceCutQuadrature;
    SurfaceCutQuadrature surfaceCutQuadrature1( surfCells, true );
    SurfaceCutQuadrature surfaceCutQuadrature2( surfCells, false );

    // FE
#ifdef LINEAR    
    typedef mat::hypel::StVenant Material;
#else
    typedef mat::hypel::NeoHookeanCompressible Material;
#endif
    typedef HyperElastic<Mesh,Material,fieldDeg> HyperElastic;
    HyperElastic hyperElastic1( mesh, E1, nu1 );
    HyperElastic hyperElastic2( mesh, E2, nu2 );

    typedef SurfaceField<BoundaryMesh,HyperElastic::Field> BoundaryField;
    BoundaryField boundaryField1( boundaryMesh, hyperElastic1.getField() );
    BoundaryField boundaryField2( boundaryMesh, hyperElastic2.getField() );

    typedef InterfaceField<SurfaceMesh,HyperElastic::Field> InterfaceField;
    InterfaceField interfaceField( surfaceMesh, hyperElastic1.getField(),
                                   hyperElastic2.getField() );

    // support sizes
    std::vector<double> supports1, supports2;
    base::cut::supportComputation( mesh, hyperElastic1.getField(), cutQuadrature1, supports1 );
    base::cut::supportComputation( mesh, hyperElastic2.getField(), cutQuadrature2, supports2 );

    std::vector<std::pair<std::size_t,VecDim> > doFLocation;
    base::dof::associateLocation( hyperElastic1.getField(), doFLocation );

    // Essential BCs
    base::dof::constrainBoundary<HyperElastic::FEBasis>(
        meshBoundary.begin(), meshBoundary.end(), mesh, hyperElastic1.getField(),
        boost::bind( &fix<dim,HyperElastic::DoF>, _1, _2, bbox ) );

    base::dof::constrainBoundary<HyperElastic::FEBasis>(
        meshBoundary.begin(), meshBoundary.end(), mesh, hyperElastic2.getField(),
        boost::bind( &fix<dim,HyperElastic::DoF>, _1, _2, bbox ) );

    // stabilisation
    base::cut::stabiliseBasis( mesh, hyperElastic1.getField(), supports1, doFLocation );
    base::cut::stabiliseBasis( mesh, hyperElastic2.getField(), supports2, doFLocation );


    // DoF numbering
    const std::size_t activeDoFs1 = 
        base::dof::numberDoFsConsecutively( hyperElastic1.getField().doFsBegin(),
                                            hyperElastic1.getField().doFsEnd() );
    const std::size_t activeDoFs2 = 
    base::dof::numberDoFsConsecutively( hyperElastic2.getField().doFsBegin(),
                                        hyperElastic2.getField().doFsEnd(), activeDoFs1 );


#ifndef LINEAR
    hyperElastic1.writeVTKFile(    baseName + "-1", 0, levelSet, cells, true );
    hyperElastic1.writeVTKFileCut( baseName + "-1", 0, levelSet, cells, true );
    hyperElastic2.writeVTKFile(    baseName + "-2", 0, levelSet, cells, false );
    hyperElastic2.writeVTKFileCut( baseName + "-2", 0, levelSet, cells, false );
#endif
    
    //--------------------------------------------------------------------------
    // Computation
    for ( unsigned step = 0; step < numLoadSteps; step++ ) {

        const double loadFactor =
            static_cast<double>( step+1 ) / static_cast<double>( numLoadSteps );

        unsigned iter = 0;
        for ( ; iter < maxIter; iter++ ) {
        
            // Solver and non-zero pattern
            base::solver::Eigen3 solver( activeDoFs1 + activeDoFs2 );
            hyperElastic1.registerInSolver(  solver );
            hyperElastic2.registerInSolver(  solver );
            interfaceField.registerInSolver( solver );

            // a(u_i,v_i)
            hyperElastic1.assembleBulk( cutQuadrature1, solver, iter );
            hyperElastic2.assembleBulk( cutQuadrature2, solver, iter );

            // traction BC
            boundaryField1.applyNeumannBoundaryConditions(
                surfaceCutQuadrature1, solver,
                boost::bind( &sideLoad<dim>, _1, _2, bbox, forceVal*loadFactor ) );
            boundaryField2.applyNeumannBoundaryConditions(
                surfaceCutQuadrature2, solver,
                boost::bind( &sideLoad<dim>, _1, _2, bbox, forceVal*loadFactor ) );

            // interface conditions
            interfaceField.applyInterfaceConditionsWeakly( surfaceQuadrature, solver,
                                                           hyperElastic1.getKernel(),
                                                           hyperElastic2.getKernel(),
                                                           E1, E2, gamma );

            solver.finishAssembly();

#ifndef LINEAR
            const double conv11 = solver.norm( 0, activeDoFs1 );
            const double conv12 = solver.norm( activeDoFs1+1, activeDoFs1+activeDoFs2 );
            std::cout << step << "  " << iter << "  "
                      << conv11 << "  " << conv12 << " ";
            if ( conv11+conv12 < tolerance * std::min( E1, E2 ) ) {
                std::cout << std::endl;
                break;
            }
#endif
            
            solver.luSolve();

            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, hyperElastic1.getField() );
            base::dof::addToDoFsFromSolver( solver, hyperElastic2.getField() );

#ifndef LINEAR
            const double conv21 = solver.norm( 0, activeDoFs1 );
            const double conv22 = solver.norm( activeDoFs1+1, activeDoFs1+activeDoFs2 );
            std::cout << conv21 << "  " << conv22  << std::endl;
            if ( conv21+conv22 < tolerance ) break;
#endif

        }

        hyperElastic1.writeVTKFile(    baseName + "-1", step+1, levelSet, cells, true );
        hyperElastic1.writeVTKFileCut( baseName + "-1", step+1, levelSet, cells, true );
        hyperElastic2.writeVTKFile(    baseName + "-2", step+1, levelSet, cells, false );
        hyperElastic2.writeVTKFileCut( baseName + "-2", step+1, levelSet, cells, false );

        // point at upper right corner
        VecDim point = base::constantVector<dim>(0.);
        point[0] = point[2] = 0.5;
        const std::size_t measureDoFID =
            base::dof::findDoFWithLocation( doFLocation, mesh, point, coordTol );

        // write out value
        std::cout << "Measured displacement: ";
        for ( unsigned d = 0; d < 1; d++  ) { // dim
            std::cout << ( hyperElastic1.getField().doFPtr( measureDoFID ) ) -> getValue(d)
                      << "  "
                      << ( hyperElastic2.getField().doFPtr( measureDoFID ) ) -> getValue(d)
                      << "  ";
        }
        std::cout << std::endl;
    }


    return 0;
}
