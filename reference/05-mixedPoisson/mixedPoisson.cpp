#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>

#include <base/io/Format.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/Field.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>

#include <base/post/ErrorNorm.hpp>
#include <heat/Laplace.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/fe/Basis.hpp>

#include <base/kernel/Measure.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/NeumannForce.hpp>



const double coordTol = 1.e-5;

#include "ReferenceSolution.hpp"
#include "GivenData.hpp"

//------------------------------------------------------------------------------
namespace ref05{
    int mixedPoisson( int argc, char * argv[] );
}


//------------------------------------------------------------------------------
/** Solve a mixed boundary value problem based on Poisson's equation
 *
 *  New features are
 *  - application of Neumann boundary conditions and body forces
 *  - use classes to describe the reference solution and the boundary
 *    value problem
 *  - convergence analysis (if executed for many meshes)
 *
 *  The picture shows the convergence behaviour of the method
 *  for various mesh sizes, a linear finite element basis
 *  \code{.cpp}
 *  const unsigned fieldDeg = 1;
 *  \endcode
 *  and a quadratic finite element basis
 *  \code{.cpp}
 *  const unsigned fieldDeg = 2;
 *  \endcode
 *  in the \f$ L_2 \f$-norm
 *  \f[
 *         \| u - u^h \|_{L_2(\Omega)}^2
 *         = \int_\Omega (u - u^h)^2 d x
 *  \f]
 *  and the \f$ H^1 \f$-semi-norm
 *  \f[
 *         | u - u^h |_{H^1(\Omega)}^2
 *         = \int_\Omega (\nabla u - \nabla u^h)^2 d x
 *  \f]
 *
 *  \image html convergence.png "Convergence diagram"
 *
 */
int ref05::mixedPoisson( int argc, char * argv[] )
{
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf \n\n";
        return -1;
    }

    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    //--------------------------------------------------------------------------
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;
    const unsigned    doFSize  = 1;

    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned dim = Mesh::Node::dim;
    
    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature and surface quadrature
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;
 
    // DOF handling
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize>           Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field field;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, field );

    // Object of reference solution and the Poisson problem
    typedef ReferenceSolution<dim>       RefSol;
    typedef GivenData<RefSol>            GivenData;
    RefSol    refSol(    3., 5., 4. );
    GivenData givenData( refSol );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, field,
                                           boost::bind( &GivenData::dirichletBC<DoF>,
                                                        &givenData, _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << "Number of dofs " << numDofs << std::endl;

    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;


    // Body force
    base::asmb::bodyForceComputation<FTB>( quadrature, solver, fieldBinder,
                                           boost::bind( &GivenData::forceFun,
                                                        &givenData, _1 ) );
#if 0
    base::asmb::bodyForceComputation2<FTB>( quadrature, solver, fieldBinder,
                                            boost::bind( &GivenData::forceFun2<Mesh::Element>,
                                                         &givenData, _1, _2 ) );
#endif 

    // Create a connectivity out of this list
    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
    }


    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, field );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SFTB;

    // Neumann boundary condition
    base::asmb::neumannForceComputation<SFTB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                               boost::bind( &GivenData::neumannBC,
                                                            &givenData,
                                                            _1, _2 ) );

    // compute stiffness matrix
    typedef heat::Laplace<FTB::Tuple> Laplace;
    Laplace laplace( 1. );
    base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                 fieldBinder, laplace );
    
    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.choleskySolve();

    // distribute results back to dofs
    base::dof::setDoFsFromSolver( solver, field );

    //--------------------------------------------------------------------------
    // output to a VTK file
    {
        // VTK Legacy
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        base::io::vtk::writePointData( vtkWriter, mesh, field, "temperature" );
        vtk.close();
    }

    //--------------------------------------------------------------------------
    // compute L2-error
    std::cout << "L2-error = "
              << base::post::errorComputation<0>(
                  quadrature, mesh, field,
                  boost::bind( &ReferenceSolution<dim>::evaluate,
                               &refSol, _1 ) )
              << '\n';

    //--------------------------------------------------------------------------
    // compute H1-error
    std::cout << "H1-error = "
              << base::post::errorComputation<1>(
                  quadrature, mesh, field,
                  boost::bind( &ReferenceSolution<dim>::evaluateGradient,
                               &refSol, _1 ) )
              << '\n';

    //--------------------------------------------------------------------------
    // Compute mesh volume
    double volume = 0.;

    base::asmb::simplyIntegrate<FTB>( quadrature, volume, fieldBinder,
                                      base::kernel::Measure<FTB::Tuple>() );
    std::cout << "Volume of mesh: " << volume << '\n';

    
    return 0;
}

//------------------------------------------------------------------------------
// Delegation
int main( int argc, char * argv[] )
{
    return ref05::mixedPoisson( argc, argv );
}
