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
#include <base/io/PropertiesParser.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>

#include <heat/Laplace.hpp>
#include <base/solver/Eigen3.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/NeumannForce.hpp>

#include "BoundaryValueProblem.hpp"

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    // Static attributes
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::TET;
    const unsigned    doFSize  = 1;
    const unsigned kernelDegEstimate = 3;

    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input \n\n";
        return -1;
    }

    // first command line argument is the input data file
    const std::string inputFile  = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double kappa;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "kappa",            kappa );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck(inp, std::cerr ), "Input error" );
        inp.close( );
    }

    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );


    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned dim = Mesh::Node::dim;
    
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature and surface quadrature
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;
 
    // DOF handling
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize>           Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field temperature;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, temperature );

    
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(), meshBoundary.end(),
                                           mesh, temperature,
                                           boost::bind( &BoundaryValueProblem<dim>::dirichleBC<DoF>,
                                                        _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( temperature.doFsBegin(), temperature.doFsEnd() );
    std::cout << "Number of dofs " << numDofs << std::endl;

    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, temperature );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;


    // Body force
    base::asmb::bodyForceComputation<FTB>( quadrature, solver, fieldBinder,
                                           boost::bind( &BoundaryValueProblem<dim>::forceFun,
                                                        _1 ) );

#if 0
    // Create a connectivity out of this list
    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(), meshBoundary.end(),
                                          mesh, boundaryMesh );

        const std::string name = baseName + ".boundary.smf";
        std::ofstream smf( name.c_str() );
        base::io::smf::writeMesh( boundaryMesh, smf );
        smf.close();
    }


    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, temperature );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SFTB;

    // Neumann boundary condition
    base::asmb::neumannForceComputation<SFTB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                               boost::bind( &BoundaryValueProblem<dim>::neumannBC,
                                                            _1, _2 ) );
#endif 
    
    // compute stiffness matrix
    typedef heat::Laplace<FTB::Tuple> Laplace;
    Laplace laplace( kappa );
    base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                 fieldBinder, laplace );
    
    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.choleskySolve();

    // distribute results back to dofs
    base::dof::setDoFsFromSolver( solver, temperature );

    //--------------------------------------------------------------------------
    // output to a VTK file
    {
        // VTK Legacy
        const std::string vtkFile = baseName + ".vtk";
        writeVTKFile( mesh, temperature, vtkFile );
    }
    
    return 0;
}
