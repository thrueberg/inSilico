#include <string>
#include <boost/lexical_cast.hpp>

#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/Field.hpp>

#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/cut/Cell.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/cut/analyticLevelSet.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/ScaledField.hpp>
#include <base/cut/Quadrature.hpp>

#include <base/fe/Basis.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <solid/HyperElastic.hpp>
#include <mat/hypel/CompNeoHookean.hpp>
#include <mat/Lame.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>

#include "generateMesh.hpp"

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
    const unsigned              nHist   = 1;

    typedef base::Unstructured<shape,geomDeg>  Mesh;
    Mesh mesh;
    {
        generateMesh( mesh, numElements, 0, xmax );
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

    typedef base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
    SurfaceMesh boundaryMesh;
    {
        base::cut::generateSurfaceMesh<Mesh,Cell>( mesh, cells, boundaryMesh );
    }

    //--------------------------------------------------------------------------
    // FE
    typedef base::fe::Basis<shape,fieldDeg>               FEBasis;
    typedef base::cut::ScaledField<FEBasis,doFSize,nHist> Field;
    typedef Field::DegreeOfFreedom                        DoF;
    Field displacement;
    base::dof::generate<FEBasis>( mesh, displacement );

    // for domain field
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, displacement );
    typedef FieldBinder::TupleBinder<1,1>::Type UU;

    // for surface field
    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, displacement );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SU;

    //--------------------------------------------------------------------------
    // Quadratures
    typedef base::cut::Quadrature<kernelDegEstimate,shape> CutQuadrature;
    CutQuadrature cutQuadrature( cells, true );

    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // compute supports, scale basis
    const std::size_t numDoFs = std::distance( displacement.doFsBegin(),
                                               displacement.doFsEnd() );
    std::vector<double> supports;
    supports.resize(  numDoFs );
    
    base::cut::supportComputation( mesh, displacement, cutQuadrature,  supports );
    displacement.scaleAndTagBasis( supports,  1.e-10 );

    // fix left end
    Field::DoFPtrIter dIter = displacement.doFsBegin();
    for ( unsigned d = 0; d < dim; d++ ) (*dIter) -> constrainValue( d, 0. );
    
    // number DoFs
    const std::size_t activeDoFsU = 
        base::dof::numberDoFsConsecutively( displacement.doFsBegin(), displacement.doFsEnd() );

    typedef mat::hypel::CompNeoHookean Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );
    typedef solid::HyperElastic<Material,UU::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    writeVTKFile( "lagrange", 0, mesh, displacement, levelSet );

    // point to check
    typedef Mesh::Node::VecDim VecDim;
    VecDim x = base::constantVector<dim>( 0. );
    x[0] = L;

    // find point in mesh
    std::pair<std::size_t,VecDim> probe;
    const bool found = base::post::findLocationInMesh( mesh, x, coordTol, 10, probe );
    VERIFY_MSG( found, "Could not find point" );

    // prepare a monitor
    base::post::Monitor<Mesh::Element,Field::Element>
        monitorD( mesh.elementPtr(         probe.first ),
                  displacement.elementPtr( probe.first ),
                  probe.second );

    //--------------------------------------------------------------------------
    // Loop over time steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < numLoadSteps; step++ ) {

        const double fFun = f *
            static_cast<double>(step+1) / static_cast<double>( numLoadSteps );

        for ( unsigned iter = 0; iter < maxIter; iter++ ) {

            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( activeDoFsU );

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
            std::cout << "* " << iter << " " << conv1 << " ";

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

        // write a vtk file
        writeVTKFile( "lagrange", step+1, mesh, displacement, levelSet );


        std::cout << step << "  " << fFun << "  ";
        monitorD.solution( std::cout );
        std::cout << std::endl;

    }

    return 0;
}
