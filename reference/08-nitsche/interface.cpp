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
#include <base/cut/tagBasis.hpp>

#include <base/Quadrature.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/kernel/Laplace.hpp>

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

#include "OneDimensionalProblem.hpp"
#include "SphericalProblem.hpp"
#include "Helper.hpp"

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
                  << " file.smf delta alpha1 alpha2 \n\n";
        return -1;
    }

    // input mesh file
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    // problem parameter
    const double delta  = boost::lexical_cast<double>( argv[2] );
    const double alpha1 = boost::lexical_cast<double>( argv[3] );
    const double alpha2 = boost::lexical_cast<double>( argv[4] );

    const double penaltyFactor = 20.0;

    //--------------------------------------------------------------------------
    const unsigned    geomDeg   = 1;
    const unsigned    dim       = SPACEDIM;
    const bool        isSigned  = true;
    
    const unsigned fieldDeg = 1;
    const unsigned doFSize  = 1;

    // use Nitsche
    const bool useNitscheTerms = true;

    //
    //typedef apps::nitsche::OneDimensionalProblem<dim> Reference;
    typedef   apps::nitsche::SphericalProblem<dim> Reference;

    //--------------------------------------------------------------------------
    // Domain mesh
    typedef apps::nitsche::UnstructuredHelper<dim,geomDeg,fieldDeg> Helper;
    typedef Helper::Mesh Mesh;

    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        Helper::readFromFile( smf, mesh );
    }

    //--------------------------------------------------------------------------
    // Compute the level set data
    typedef base::cut::LevelSet<dim> LevelSet;
    std::vector<LevelSet> levelSet;
    base::cut::AnalyticSurface<dim>::Type as =
        boost::bind( &Reference::interface, _1, delta, _2 );
    base::cut::analyticLevelSet( mesh, as, isSigned, levelSet );


    //--------------------------------------------------------------------------
    // FE
    typedef Helper::FEBasis FEBasis;
    typedef base::cut::ScaledField<FEBasis,doFSize> Field;
    Field fieldIn, fieldOut;

    base::dof::generate<FEBasis>( mesh, fieldIn  );
    base::dof::generate<FEBasis>( mesh, fieldOut );
    

    //--------------------------------------------------------------------------
    // Cut
    typedef base::cut::Cell<Helper::shape> Cell;
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
        Helper::meshBoundary( meshBoundary, mesh );

        // generate a mesh from that list (with a filter)
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh,
                                          boost::bind( &Reference::dirichletBoundary, _1 ) );
    }

    // from interface
    base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, interfaceMesh );

    //--------------------------------------------------------------------------
    // Quadratures
    const unsigned kernelDegEstimate = 5;
    // for domain
    typedef base::cut::Quadrature<kernelDegEstimate,Helper::shape> CutQuadrature;
    CutQuadrature cutQuadratureIn(  cells, true  );
    CutQuadrature cutQuadratureOut( cells, false );
    // for surface
    typedef base::SurfaceQuadrature<kernelDegEstimate,Helper::shape> SurfaceQuadrature;
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
    base::cut::supportComputation( mesh, fieldOut, cutQuadratureOut, supportsOut);

    fieldIn.scaleAndTagBasis(  supportsIn,  1.e-10 );
    fieldOut.scaleAndTagBasis( supportsOut, 1.e-10 );
    // base::cut::tagBasis( fieldIn,  supportsIn, 1.e-10 );
    // base::cut::tagBasis( fieldOut, supportsOut, 1.e-10 );
    
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
    // Body force
    base::asmb::bodyForceComputation<FTB1>( cutQuadratureIn,  solver, fieldBinder,
                                            boost::bind( &Reference::forceFun, _1 ) );
    base::asmb::bodyForceComputation<FTB2>( cutQuadratureOut, solver, fieldBinder,
                                            boost::bind( &Reference::forceFun, _1 ) );


    //--------------------------------------------------------------------------
    // Stiffness matrix
    typedef  base::kernel::Laplace<FTB1::Tuple> LaplaceIn;
    LaplaceIn  laplaceIn(  alpha1  );
    typedef  base::kernel::Laplace<FTB2::Tuple> LaplaceOut;
    LaplaceOut laplaceOut( alpha2 );

    message( "Stiffness" );
    base::asmb::stiffnessMatrixComputation<FTB1>( cutQuadratureIn, solver,
                                                  fieldBinder, laplaceIn );
    base::asmb::stiffnessMatrixComputation<FTB2>( cutQuadratureOut, solver,
                                                  fieldBinder, laplaceOut );

    // boundary conditions
#if 1
    {
        message("Boundary Penalty");
        base::nitsche::OuterBoundary ob1( alpha1);
        base::nitsche::OuterBoundary ob2( alpha2);
        
        base::nitsche::penaltyLHS<STB11>( surfaceQuadrature, solver,
                                          boundaryFieldBinder, ob1, penaltyFactor );
        base::nitsche::penaltyLHS<STB22>( surfaceQuadrature, solver,
                                          boundaryFieldBinder, ob2, penaltyFactor );
        base::nitsche::penaltyRHS<STB11>( surfaceQuadrature, solver, boundaryFieldBinder, 
                                          boost::bind( &Reference::dirichlet, _1,
                                                       delta, alpha1, alpha2 ), ob1,
                                          penaltyFactor );
        base::nitsche::penaltyRHS<STB22>( surfaceQuadrature, solver, boundaryFieldBinder,
                                          boost::bind( &Reference::dirichlet, _1,
                                                       delta, alpha1, alpha2 ), ob2,
                                          penaltyFactor );

        if ( useNitscheTerms ) {

            message("Boundary Nitsche");
            base::nitsche::energyLHS<STB11>( laplaceIn, surfaceQuadrature, solver,
                                             boundaryFieldBinder, ob1 );
            base::nitsche::energyLHS<STB22>( laplaceOut, surfaceQuadrature, solver,
                                             boundaryFieldBinder, ob2 );

            base::nitsche::energyRHS<STB11>( laplaceIn, surfaceQuadrature, solver,
                                             boundaryFieldBinder,
                                             boost::bind( &Reference::dirichlet, _1,
                                                          delta, alpha1, alpha2 ), ob1 );
            
            base::nitsche::energyRHS<STB22>( laplaceOut, surfaceQuadrature, solver,
                                             boundaryFieldBinder,
                                             boost::bind( &Reference::dirichlet, _1,
                                                          delta, alpha1, alpha2 ), ob2 );
        }
        
    }

    // interface conditions
    {
        message("Interface Penalty");
        base::nitsche::ImmersedInterface<Cell> ip(  alpha1,  alpha2, cells );

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
            
            //   kappa1
            base::nitsche::energyLHS<STB11>( laplaceIn, surfaceQuadrature, solver,
                                             interfaceFieldBinder, ip, true, true  );
            //  -kappa1
            base::nitsche::energyLHS<STB21>( laplaceIn, surfaceQuadrature, solver,
                                             interfaceFieldBinder, ip, true, false ); 

            //   kappa2
            base::nitsche::energyLHS<STB12>( laplaceOut, surfaceQuadrature, solver,
                                             interfaceFieldBinder, ip, false, true );
            //  -kappa2
            base::nitsche::energyLHS<STB22>( laplaceOut, surfaceQuadrature, solver,
                                             interfaceFieldBinder, ip, false, false ); 

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

    //--------------------------------------------------------------------------
    // Compute the L2-error
    const double leftErrorL2 = 
        base::post::errorComputation<0>( cutQuadratureIn, mesh, fieldIn,
                                         boost::bind( &Reference::sol, _1,
                                                      delta, alpha1, alpha2 ) );

    const double rightErrorL2 =
        base::post::errorComputation<0>( cutQuadratureOut, mesh, fieldOut,
                                         boost::bind( &Reference::sol, _1,
                                                      delta, alpha1, alpha2 ) );

    const double errorL2 = std::sqrt( leftErrorL2  * leftErrorL2 +
                                      rightErrorL2 * rightErrorL2);

    //--------------------------------------------------------------------------
    // Compute the H1-error
    const double leftErrorH1 = 
        base::post::errorComputation<1>( cutQuadratureIn, mesh, fieldIn,
                                         boost::bind( &Reference::solDeriv, _1,
                                                      delta, alpha1, alpha2 ) );

    const double rightErrorH1 = 
        base::post::errorComputation<1>( cutQuadratureOut, mesh, fieldOut,
                                         boost::bind( &Reference::solDeriv, _1,
                                                      delta, alpha1, alpha2 ) );

    const double errorH1 = std::sqrt( leftErrorH1  * leftErrorH1 +
                                      rightErrorH1 * rightErrorH1 );


    //--------------------------------------------------------------------------
    // Surface errors
    const double boundError =
        base::post::errorComputation<0>( surfaceQuadrature, boundaryMesh, fieldOut,
                                         boost::bind( &Reference::sol, _1,
                                                      delta, alpha1, alpha2 ) );

    const double interfaceError =
        base::post::errorComputation<0>( surfaceQuadrature, interfaceMesh, fieldOut,
                                         boost::bind( &Reference::sol, _1,
                                                      delta, alpha1, alpha2 ) );



#ifdef VERBOSE   // decide error verbosity
    
    // explicit output
    std::cout << "L2-Error = " << errorL2
              << " = sqrt("
              << leftErrorL2 << "^2 + " << rightErrorL2 << "^2 )\n"
              << "H1-Error = " << errorH1
              << " = sqrt("
              << leftErrorH1 << "^2 + " << rightErrorH1 << "^2 )\n";


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
    std::cout << "  " << errorL2 << "  " << errorH1 << "  "
              << boundError << "  " << interfaceError << "  "
              << solveTime << "  " << numCGIter << "\n";
#endif


#endif // decide to solve and write

    message( "Post-process time = " + detailed.print() );
    message( "---------------------------------------" );
    
    message( "Total time = " + total.print() );
    
    return 0;
}
