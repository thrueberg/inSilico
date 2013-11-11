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
#include <base/dof/setField.hpp>

#include <heat/Laplace.hpp>
#include <base/solver/Eigen3.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

#include "BoundaryValueProblem.hpp"

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    //--------------------------------------------------------------------------
    // Static attributes
    const unsigned    geomDeg           = 1;
    const unsigned    fieldDeg          = 2;
    const base::Shape shape             = base::TET;
    const unsigned    doFSize           = 1;
    const unsigned    kernelDegEstimate = 3;
    const unsigned    tiOrder           = 1;

    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input \n\n";
        return -1;
    }

    // first command line argument is the input data file
    const std::string inputFile  = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double kappa, rho, numSteps, stepSize;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "kappa",            kappa );
        prop.registerPropertiesVar( "rho",              rho );
        prop.registerPropertiesVar( "numSteps",         numSteps );
        prop.registerPropertiesVar( "stepSize",         stepSize );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck(inp) , "Input error" );
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
 
    // time integration
    typedef base::time::AdamsMoulton<tiOrder> MSM;
    const unsigned nHist = MSM::numSteps;

    // DOF handling
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>     Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field temperature;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, temperature );

    
    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, temperature );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;
    
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );


    // set initial condition
    base::dof::setField( mesh, temperature,
                         boost::bind( &BoundaryValueProblem<dim>::initialState<DoF>,
                                      _1, _2 ) );
    // pass to history 
    std::for_each( temperature.doFsBegin(), temperature.doFsEnd(),
                   boost::bind( &DoF::pushHistory, _1 ) );

    // write initial state
    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( 0 ) + ".vtk";
    writeVTKFile( mesh, temperature, vtkFile );

    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(), meshBoundary.end(),
                                           mesh, temperature,
                                           boost::bind( &BoundaryValueProblem<dim>::dirichleBC<DoF>,
                                                        _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( temperature.doFsBegin(), temperature.doFsEnd() );
    std::cout << "Number of dofs " << numDofs << std::endl;


    //--------------------------------------------------------------------------
    // Time loop
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < numSteps; step++ ) {

        std::cout << "Step " << step << "\n";
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        // compute stiffness matrix
        typedef heat::Laplace<FTB::Tuple> Laplace;
        Laplace laplace( kappa );
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                     fieldBinder, laplace );

        // compute inertia terms, d/dt, due to time integration
        base::time::computeInertiaTerms<FTB,MSM>( quadrature, solver,
                                                  fieldBinder, stepSize, step,
                                                  rho );

        // compute history of residual forces due to time integration
        base::time::computeResidualForceHistory<FTB,MSM>( laplace, 
                                                          quadrature, solver,
                                                          fieldBinder, step );

        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.choleskySolve();

        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, temperature );


        // pass to history 
        std::for_each( temperature.doFsBegin(), temperature.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

        // write output
        const std::string vtkFile =
            baseName + "." + base::io::leadingZeros( step+1 ) + ".vtk";
        writeVTKFile( mesh, temperature, vtkFile );
    }

    return 0;
}
