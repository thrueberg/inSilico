// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

// mesh-related includes
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>

// input/output includes
#include <base/io/Format.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

// quadrature
#include <base/Quadrature.hpp>

// FE Basis
#include <base/fe/Basis.hpp>

// Field and degrees of freedom
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>

// kernel of Laplace equation
#include <heat/Laplace.hpp>

// linear system solver
#include <base/solver/Eigen3.hpp>

// Assembly module
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/NeumannForce.hpp>
// time integration
#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

// local inclusion of the boundary value problem's definitions
#include "BoundaryConditions.hpp"

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // Tell user how to call this app
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input \n\n";
	return -1;
    }
    
    // Static attributes
    const unsigned geomDeg  			= 1;
    const unsigned fieldDeg 			= 2;
    const base::Shape shape    			= base::TET;
    const unsigned doFSize  			= 1;
    const unsigned kernelDegEstimate 	= 3;
    const unsigned tiOrder           	= 2;

    // first command line argument is the input data file
    const std::string inputFile  = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double diff_default;
	double D1, D2, stepSize;
	unsigned numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "Default",     diff_default );
		prop.registerPropertiesVar( "D1",         		D1 );
		prop.registerPropertiesVar( "D2",         		D2 );
        prop.registerPropertiesVar( "numSteps",    numSteps );
        prop.registerPropertiesVar( "stepSize",    stepSize );


       // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck(inp, std::cerr ), "Input error" );
        inp.close( );
    }

    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
    // create a mesh from file
	std::cout << "Create Mesh from file" << std::endl;
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned dim = Mesh::Node::dim;
    
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature and surface quadrature
	std::cout << "Quadrature" << std::endl;

    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;
 
 	// time integration
    typedef base::time::BDF<tiOrder> MSM;
    //typedef base::time::AdamsMoulton<tiOrder> MSM;
    const unsigned nHist = MSM::numSteps;


    // DOF handling
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>     Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field concetration;

    // generate DoFs from mesh
	std::cout << "Generate DoFs from mesh" << std::endl;
    base::dof::generate<FEBasis>( mesh, concetration );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(), meshBoundary.end(),
                                           mesh, concetration,
                                           boost::bind( &BoundaryConditions<dim>::dirichleBC<DoF>,
                                                        _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( concetration.doFsBegin(), concetration.doFsEnd() );
    std::cout << "Number of dofs " << numDofs << std::endl;

    // Bind the fields together
	std::cout << "Bind Fields" << std::endl;
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, concetration );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

// Time Loop

std::cout << "\n\nTime loop:\n" << std::endl;


for ( unsigned step = 0; step < numSteps; step++ ) {
	
    std::cout << "\033[2K\r\t"<< "Time:" << (step+1)*stepSize << "\tStep: " 
				<< step+1 << "/" << numSteps << std::flush;
    
    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );


    // compute stiffness matrix - diffusion

	typedef heat::Laplace<FTB::Tuple> Laplace;
    Laplace laplace( diff_default );
    boost::function< double( const Mesh::Element *,
                             const Mesh::Element::GeomFun::VecDim & ) >
        diffusionFun = boost::bind( &diffusionConstant<Mesh::Element>, _1, _2,
                                    D1, D2);
    laplace.setConductivityFunction( diffusionFun );

    base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                 fieldBinder, laplace );

	// compute inertia terms, d/dt, due to time integration
    base::time::computeInertiaTerms<FTB,MSM>( quadrature, solver,
                                                  fieldBinder, stepSize, step,
                                                  1.0 );


    // compute history of residual forces due to time integration
    base::time::computeResidualForceHistory<FTB,MSM>( laplace, 
                                                      quadrature, solver,
                                                      fieldBinder, step );

    
    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.choleskySolve();

    // distribute results back to dofs
    base::dof::setDoFsFromSolver( solver, concetration );
    
	// pass to history 
    std::for_each( concetration.doFsBegin(), concetration.doFsEnd(),
                    boost::bind( &DoF::pushHistory, _1 ) );


    // output to a VTK file
    const std::string vtkFile = 
		baseName + base::io::leadingZeros( step+1 )+ ".vtk";
    writeVTKFile( mesh, concetration, vtkFile );
};    
    return 0;
}
