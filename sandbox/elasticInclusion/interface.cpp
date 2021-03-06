#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/Unstructured.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/io/Format.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/cut/ScaledField.hpp>

#include <base/Quadrature.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/solver/Eigen3.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/Lame.hpp>
#include <solid/HyperElastic.hpp>
#include <solid/Stress.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>
#include <base/nitsche/Parameters.hpp>

#include <base/cut/LevelSet.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/extractMeshFromCutCells.hpp>

#include <base/post/ErrorNorm.hpp>

#include <base/asmb/SimpleIntegrator.hpp>

#include <base/auxi/Timer.hpp>




template<unsigned DIM>
bool interface( const typename base::Vector<DIM>::Type& x,
                const double radius,
                typename base::Vector<DIM>::Type& xClosest )
{
    const typename base::Vector<DIM>::Type centre = base::constantVector<DIM>( 0.5 );
    const typename base::Vector<DIM>::Type diff   = x - centre;
    const double r      = diff.norm();

    // catch the centre of the sphere
    if ( r    > 1.e-10 )
        xClosest    = centre + (radius / r) * diff;
    else 
        xClosest    = centre + (radius /centre.norm()) * centre;
    
    return ( r < radius );
}

// description of the boundary (left or right in x_1-direction)
template<unsigned DIM>
bool dirichletBoundary( const typename base::Vector<DIM>::Type& x )
{
    const double eps = 1.e-8;
    
    if ( ( std::abs( x[0] - 0.0 ) < eps ) or
         ( std::abs( x[0] - 1.0 ) < eps ) ) return true;
    
    return false;
}

// Dirichlet function on boundary (zero anyway)
template<unsigned DIM>
typename base::Vector<DIM>::Type dirichlet( const typename base::Vector<DIM>::Type& x )
{
    typename base::Vector<DIM>::Type result =
        base::constantVector<DIM>( 0. );
    
    const double eps = 1.e-8;
    
    if ( std::abs( x[0] - 1.0 ) < eps ) result[0] = 1.0;

    return result;
}

//#define VERBOSE

void message( const std::string& msg )
{
#ifdef VERBOSE
    std::cout << msg << "\n";
#endif
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    base::auxi::Timer total, detailed;
    
    //--------------------------------------------------------------------------
    // User input
    if ( argc != 5 ) {
        std::cout << "Usage:  " << argv[0]
                  << " file.smf radius E1 E2 \n\n";
        return -1;
    }

    // input mesh file
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    // problem parameter
    const double radius = boost::lexical_cast<double>( argv[2] );
    const double E1     = boost::lexical_cast<double>( argv[3] );
    const double E2     = boost::lexical_cast<double>( argv[4] );

    const double penaltyFactor = 20.0;

    // material stuff
    const double nu = 0.3;
    const double lambda1 = mat::Lame::lambda( E1, nu );
    const double lambda2 = mat::Lame::lambda( E2, nu );
    const double mu1     = mat::Lame::mu( E1, nu );
    const double mu2     = mat::Lame::mu( E2, nu );

    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
    const unsigned    dim       = SPACEDIM;
    const base::Shape shape     = base::SimplexShape<dim>::value;
    const bool        isSigned  = true;
    
    const unsigned fieldDeg = 1;
    const unsigned doFSize  = dim;

    // use Nitsche
    const bool useNitscheTerms = true;

    //--------------------------------------------------------------------------
    // Domain mesh
    typedef base::Unstructured<shape,geomDeg> Mesh;

    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
    }

    //--------------------------------------------------------------------------
    // Compute the level set data
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::AnalyticSurface<dim>::Type as =
        boost::bind( &interface<dim>, _1, radius, _2 );
    base::cut::analyticLevelSet( mesh, as, isSigned, levelSet );


    //--------------------------------------------------------------------------
    // FE
    typedef base::fe::Basis<shape,fieldDeg>         FEBasis;
    typedef base::cut::ScaledField<FEBasis,doFSize> Field;
    Field fieldIn, fieldOut;

    base::dof::generate<FEBasis>( mesh, fieldIn  );
    base::dof::generate<FEBasis>( mesh, fieldOut );
    

    //--------------------------------------------------------------------------
    // Cut
    typedef base::cut::Cell<shape> Cell;
    std::vector<Cell> cells;
    base::cut::generateCutCells( mesh, levelSet, cells );

    //--------------------------------------------------------------------------
    //  surface meshes
    typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
    SurfaceMesh boundaryMesh, interfaceMesh;

    // from boundary
    {
        // identify list of element boundary faces
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

        // generate a mesh from that list (with a filter)
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh,
                                          boost::bind( &dirichletBoundary<dim>, _1 ) );
    }

    // from interface
    base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, interfaceMesh );

    //--------------------------------------------------------------------------
    // Quadratures
    const unsigned kernelDegEstimate = 5;
    // for domain
    typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
    CutQuadrature cutQuadratureIn(  cells, true  );
    CutQuadrature cutQuadratureOut( cells, false );
    // for surface
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    //--------------------------------------------------------------------------
    // Bind fields

    // for domain fields
    typedef base::asmb::FieldBinder<Mesh,Field,Field> FieldBinder;
    typedef FieldBinder::TupleBinder<1,1>::Type FTB1;
    typedef FieldBinder::TupleBinder<2,2>::Type FTB2;
    FieldBinder fieldBinder(  mesh, fieldIn, fieldOut );

    // for surface fields
    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Field,Field> SurfaceFieldBinder;
    typedef SurfaceFieldBinder::TupleBinder<1,1>::Type STB11;
    typedef SurfaceFieldBinder::TupleBinder<2,2>::Type STB22;
    typedef SurfaceFieldBinder::TupleBinder<1,2>::Type STB12;
    typedef SurfaceFieldBinder::TupleBinder<2,1>::Type STB21;
    SurfaceFieldBinder   boundaryFieldBinder(  boundaryMesh, fieldIn, fieldOut );
    SurfaceFieldBinder interfaceFieldBinder(  interfaceMesh, fieldIn, fieldOut );

    // compute supports, scale basis
    const std::size_t numDoFs = std::distance( fieldIn.doFsBegin(), fieldIn.doFsEnd() );
    std::vector<double> supportsIn, supportsOut;
    supportsIn.resize(  numDoFs );
    supportsOut.resize( numDoFs );
    
    base::cut::supportComputation( mesh, fieldIn,  cutQuadratureIn,  supportsIn );
    base::cut::supportComputation( mesh, fieldOut, cutQuadratureOut, supportsOut );

    fieldIn.scaleAndTagBasis(  supportsIn,  1.e-10 );
    fieldOut.scaleAndTagBasis( supportsOut, 1.e-10 );
    //fieldIn.tagBasis( supportsIn, 1.e-10 );
    //fieldOut.tagBasis( supportsOut, 1.e-10 );
    
    // number DoFs, create solver
    const std::size_t activeDoFsIn = 
        base::dof::numberDoFsConsecutively( fieldIn.doFsBegin(), fieldIn.doFsEnd() );
    const std::size_t activeDoFsOut = 
        base::dof::numberDoFsConsecutively( fieldOut.doFsBegin(),
                                            fieldOut.doFsEnd(), activeDoFsIn );

    typedef base::solver::Eigen3 Solver;
    Solver solver( activeDoFsIn + activeDoFsOut );

    message( "Preprocessing time = " + detailed.print() );
    detailed.reset();
    

    //--------------------------------------------------------------------------
    // Stiffness matrix
    typedef mat::hypel::StVenant Material;
    Material material1( lambda1, mu1 );
    Material material2( lambda2, mu2 );

    // matrix kernel
    typedef solid::HyperElastic<Material,FTB1::Tuple> HyperElastic1;
    HyperElastic1 hyperElastic1( material1 );
    typedef solid::HyperElastic<Material,FTB2::Tuple> HyperElastic2;
    HyperElastic2 hyperElastic2( material2 );

    message( "Stiffness" );
    base::asmb::stiffnessMatrixComputation<FTB1>( cutQuadratureIn, solver,
                                                  fieldBinder, hyperElastic1 );
    base::asmb::stiffnessMatrixComputation<FTB2>( cutQuadratureOut, solver,
                                                  fieldBinder, hyperElastic2 );

    // boundary conditions
#if 1
    {
        message("Boundary Penalty");
        base::nitsche::OuterBoundary ob1( E1);
        base::nitsche::OuterBoundary ob2( E2);
        
        base::nitsche::penaltyLHS<STB11>( surfaceQuadrature, solver,
                                          boundaryFieldBinder, ob1, penaltyFactor );
        base::nitsche::penaltyLHS<STB22>( surfaceQuadrature, solver,
                                          boundaryFieldBinder, ob2, penaltyFactor );
        base::nitsche::penaltyRHS<STB11>( surfaceQuadrature, solver, boundaryFieldBinder, 
                                          boost::bind( &dirichlet<dim>, _1), ob1,
                                          penaltyFactor );
        base::nitsche::penaltyRHS<STB22>( surfaceQuadrature, solver, boundaryFieldBinder,
                                          boost::bind( &dirichlet<dim>, _1), ob2,
                                          penaltyFactor );

        if ( useNitscheTerms ) {
        
            message("Boundary Nitsche");
            base::nitsche::energyLHS<STB11>( hyperElastic1, surfaceQuadrature, solver,
                                             boundaryFieldBinder, ob1 );
            base::nitsche::energyLHS<STB22>( hyperElastic2, surfaceQuadrature, solver,
                                             boundaryFieldBinder, ob2 );

            base::nitsche::energyRHS<STB11>( hyperElastic1, surfaceQuadrature, solver,
                                             boundaryFieldBinder,
                                             boost::bind( &dirichlet<dim>, _1), ob1 );
            base::nitsche::energyRHS<STB22>( hyperElastic2, surfaceQuadrature, solver,
                                             boundaryFieldBinder,
                                             boost::bind( &dirichlet<dim>, _1), ob2 );
        }
        
    }

    // interface conditions
    {
        message("Interface Penalty");
        base::nitsche::ImmersedInterface<Cell> ip(  E1,  E2, cells );

        base::nitsche::penaltyLHS<STB11>( surfaceQuadrature, solver,
                                          interfaceFieldBinder, ip, penaltyFactor );
        base::nitsche::penaltyLHS<STB22>( surfaceQuadrature, solver,
                                          interfaceFieldBinder, ip, penaltyFactor );
        
        base::nitsche::penaltyLHS<STB12>( surfaceQuadrature, solver, interfaceFieldBinder,
                                          ip, -penaltyFactor );
        base::nitsche::penaltyLHS<STB21>( surfaceQuadrature, solver, interfaceFieldBinder,
                                          ip, -penaltyFactor );

        if ( useNitscheTerms ) {
            message("Interface Nitsche");
            
            base::nitsche::energyLHS<STB11>(
                hyperElastic1, surfaceQuadrature, solver,
                interfaceFieldBinder, ip, true, true  ); //  kappa1
        
            base::nitsche::energyLHS<STB21>(
                hyperElastic1, surfaceQuadrature, solver,
                interfaceFieldBinder, ip, true, false ); //  -kappa1
            
            base::nitsche::energyLHS<STB12>(
                hyperElastic2, surfaceQuadrature, solver,
                interfaceFieldBinder, ip, false, true );  //  kappa2
        
            base::nitsche::energyLHS<STB22>(
                hyperElastic2, surfaceQuadrature, solver,
                interfaceFieldBinder, ip, false, false ); // -kappa2
            
        }
    }
#endif
    
    //--------------------------------------------------------------------------
    // Solve and distribute
    solver.finishAssembly();

    message( "Assembly time = " + detailed.print() );
    detailed.reset();

#ifdef VERBOSE
    {
        solver.systemInfo( std::cout );

        std::ofstream mat( "matrix" );
        solver.debugLHS( mat );

        std::ofstream vec( "vector" );
        solver.debugRHS( vec );
    }
#endif

    
#if 1 // decide to solve and post-process

    
    //solver.choleskySolve();
    const unsigned numCGIter = solver.cgSolve();

    //message( "Solve time = " + detailed.print() );
    const double solveTime = detailed.seconds(); 
    detailed.reset();

    base::dof::setDoFsFromSolver( solver, fieldIn  );
    base::dof::setDoFsFromSolver( solver, fieldOut );


    //--------------------------------------------------------------------------
    // Extract distances, closestPoints and location flags from level set data
    {
    }
    
    //--------------------------------------------------------------------------
    {
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        vtkWriter.writeUnstructuredGrid( mesh );
        {
            std::vector<double> distances;
            std::transform( levelSet.begin(), levelSet.end(),
                            std::back_inserter( distances ),
                            boost::bind( &LevelSet::getSignedDistance, _1 ) );
            vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
        }
        
        {
            std::vector<bool>   location;
            std::transform( levelSet.begin(), levelSet.end(),
                            std::back_inserter( location ),
                            boost::bind( &LevelSet::isInterior, _1 ) );
            vtkWriter.writePointData( location.begin(), location.end(), "location" );
        }

        vtkWriter.writePointData( supportsIn.begin(),  supportsIn.end(),  "suppIn" );
        vtkWriter.writePointData( supportsOut.begin(), supportsOut.end(), "suppOut" );

        base::io::vtk::writePointData( vtkWriter, mesh, fieldIn,  "fieldIn"  );
        base::io::vtk::writePointData( vtkWriter, mesh, fieldOut, "fieldOut" );
        
        vtk.close();
    }

#ifdef VERBOSE   // decide error verbosity
    
    // integrate
    double volumeIn = 0.;
    base::asmb::simplyIntegrate<FTB1>( cutQuadratureIn, volumeIn, fieldBinder,
                                       base::kernel::Measure<FTB1::Tuple>() );


    double volumeOut = 0.;
    base::asmb::simplyIntegrate<FTB2>( cutQuadratureOut, volumeOut, fieldBinder,
                                      base::kernel::Measure<FTB2::Tuple>() );
    
    std::cout << "Volume of mesh: " << volumeIn << " + " << volumeOut
              << " = " << volumeIn + volumeOut
              << '\n';

    

#else
    // for convergence analysis
    std::cout << solveTime << "  " << numCGIter << "\n";
#endif


#endif // decide to solve and write

    message( "Post-process time = " + detailed.print() );
    message( "---------------------------------------" );
    
    message( "Total time = " + total.print() );
    
    return 0;
}
