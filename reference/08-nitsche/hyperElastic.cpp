#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
#include <boost/lexical_cast.hpp> // lexical cast between objects

//------------------------------------------------------------------------------
#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/Quadrature.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/post/ErrorNorm.hpp>

#include <base/nitsche/Parameters.hpp>
#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>

#include <base/io/Format.hpp>

#include <mat/Lame.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/hypel/StVenant.hpp>
#include <solid/HyperElastic.hpp>

#include "Helper.hpp"

STATIC_ASSERT_MSG( (SPACEDIM>1) and (SPACEDIM<4), "Inapt choice of dimension" );

const double coordTol = 1.e-6;

//----------------------------------------------------------------------
/** Convenience function for the LHS term of the penalty method */

//--------------------------------------------------------------------------
template<unsigned DIM>
bool boundaryFilter( const typename base::Vector<DIM>::Type& x )
{
    // extract top and bottom boundaries
    if ( std::abs( x[DIM-1] - 0. ) < coordTol ) return true;
    if ( std::abs( x[DIM-1] - 1. ) < coordTol ) return true;
    return false;
}

//--------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM>::Type 
dirichlet( const typename base::Vector<DIM>::Type& x,
           const double value )
{
    typename base::Vector<DIM>::Type u = base::constantVector<DIM>( 0. );
    if ( x[DIM-1] > 0.5 ) u[DIM-1] = value;
    
    return u;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    const unsigned dim      = SPACEDIM;
    const unsigned geomDeg  = 1;
    const unsigned fieldDeg = 1;

    const unsigned maxIter  = 100;
    const double  tolerance = 1.e-8;

    const double   E  = 1000.;
    const double   nu = 0.3;

    typedef mat::hypel::NeoHookeanCompressible Material;
    //typedef mat::hypel::StVenant Material;

#ifdef STRUCTURED
    typedef apps::nitsche::StructuredHelper<dim,geomDeg,fieldDeg>   Helper;
#else
    typedef apps::nitsche::UnstructuredHelper<dim,geomDeg,fieldDeg> Helper;
#endif
    
    // Check the number of input arguments
    if ( argc != 3 ) { 
        std::cout << "Usage:  " << argv[0] << " file" << Helper::suffix()
                  << "   appplied-disp \n";
        return 0;
    }
    
    // input arguments
    const std::string fileName = boost::lexical_cast<std::string>( argv[1] );
    const double      disp     = boost::lexical_cast<double>(      argv[2] );
        
    // extract the basename of the file
    const std::string baseName =
        fileName.substr( 0, fileName.find( Helper::suffix() ) );

    //--------------------------------------------------------------------------
    const double penaltyFactor = 50.;
    const unsigned doFSize     = dim;

    //--------------------------------------------------------------------------
    typedef Helper::Mesh Mesh;
    Mesh mesh;
    {
        std::ifstream smgf( fileName.c_str() );
        Helper::readFromFile( smgf, mesh );
    }

    //--------------------------------------------------------------------------
    // Quadrature 
    const unsigned kernelDegEstimate = 5;
    typedef base::Quadrature<kernelDegEstimate,Helper::shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,Helper::shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // Finite element basis
    typedef Helper::FEBasis FEBasis;
    
    // DOF handling
    typedef base::Field<FEBasis,doFSize>           Field;
    Field field;
    base::dof::generate<FEBasis>( mesh, field );

    // Number of DoFs 
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    //std::cout << " Number of dofs " << numDofs << std::endl;

    //--------------------------------------------------------------------------
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    Helper::meshBoundary( meshBoundary, mesh );

    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                      meshBoundary.end(),
                                      mesh, boundaryMesh,
                                      boost::bind( &boundaryFilter<dim>, _1 ) );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;


    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, field );
    typedef SurfaceFieldBinder::TupleBinder<1,1>::Type STB;

    // material object
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    // matrix kernel
    typedef solid::HyperElastic<Material,FTB::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );
        
    // create table for writing the convergence behaviour of the nonlinear solves
    base::io::Table<3>::WidthArray widths = {{ 2, 10, 10 }};
    base::io::Table<3> table( widths );
    table % "iteration" % "|F|"  % "|x|";
    std::cout << "#" << table;


    for ( unsigned iter = 0; iter < maxIter; iter++ ) {

        table % iter;

        
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        //----------------------------------------------------------------------
        // Stiffness matrix
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                     fieldBinder, hyperElastic );

        base::asmb::computeResidualForces<FTB>( quadrature, solver,
                                                fieldBinder,
                                                hyperElastic );

        // Dirichlet BCs a la Nitsche
        base::nitsche::OuterBoundary ob( E * penaltyFactor );
        base::nitsche::penaltyLHS<STB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                        ob, penaltyFactor );
        
        base::nitsche::penaltyRHS<STB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                        boost::bind( &dirichlet<dim>, _1, disp ), ob,
                                        penaltyFactor );

        base::nitsche::primalEnergyLHS<STB>( hyperElastic, surfaceQuadrature, solver,
                                             surfaceFieldBinder, ob );
        
        base::nitsche::dualEnergyLHS<STB>(   hyperElastic, surfaceQuadrature, solver,
                                             surfaceFieldBinder, ob );
        
        base::nitsche::energyRHS<STB>( hyperElastic, surfaceQuadrature, solver,
                                       surfaceFieldBinder,
                                       boost::bind( &dirichlet<dim>, _1, disp ), ob );

        base::nitsche::energyResidual<STB>( hyperElastic, surfaceQuadrature, solver,
                                            surfaceFieldBinder, ob );

        // Finalise assembly
        solver.finishAssembly();

        // norm of residual 
        const double conv1 = solver.norm();
        table % conv1;

        // convergence via residual norm
        if ( conv1 < tolerance * E ) { // note the tolerance multiplier
            std::cout << table;
            break;
        }

        // Solve
        //solver.choleskySolve();
        solver.superLUSolve();
            
        // distribute results back to dofs
        base::dof::addToDoFsFromSolver( solver, field );

        //--------------------------------------------------------------------------
        // VTK file of the input mesh
        {
            const std::string vtkMeshFileName = baseName + "."
                + base::io::leadingZeros( iter ) + ".vtk";
            std::ofstream vtk( vtkMeshFileName.c_str() );
            
            Helper::writeVTK( vtk, mesh, field );
            vtk.close();
        }

        // norm of displacement increment
        const double conv2 = solver.norm();
        table % conv2;
        std::cout << table;
            
        // convergence via increment
        if ( conv2 < tolerance ) break;


    } // end of non-linear iterations
    
    return 0;
}
