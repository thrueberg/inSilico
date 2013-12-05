// system header
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/Quadrature.hpp>

#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/setField.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/FieldBinder.hpp>

#include <base/solver/Eigen3.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <base/post/Monitor.hpp>

#include <mat/thermal/IsotropicConstant.hpp>
#include <heat/Static.hpp>
#include <heat/Convection.hpp>

#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void velocityFun( const typename base::Vector<DIM,double>::Type& x,
                  DOF* doFPtr ) 
{
    typename DOF::ValueArray v;
    for ( unsigned d = 0; d < DOF::size; d++ ) v[d] = 0.;

    v[0] = 2.;
    v[1] = 1.;

    for ( unsigned d = 0; d < DOF::size; d++ ) doFPtr -> setValue( d, v[d] );
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM,double>::Type& x,
                  DOF* doFPtr ) 
{
    const double tol = 1.e-5;
    const bool left   = ( std::abs( x[0] - 0. ) < tol );
    const bool right  = ( std::abs( x[0] - 1. ) < tol );
    const bool bottom = ( std::abs( x[1] - 0. ) < tol );
    const bool top    = ( std::abs( x[1] - 1. ) < tol );

    if ( doFPtr -> isActive(0) ) {

        if ( right or bottom )
            doFPtr -> constrainValue( 0, 0. );
        else if ( left or top )
            doFPtr -> constrainValue( 0, 1. );

    }
}

//------------------------------------------------------------------------------
// output to a VTK file
template<typename MESH, typename TEMP>
void writeVTKFile( const std::string& baseName,
                   const unsigned     step,
                   const MESH&        mesh,
                   const TEMP&        temperature )
{
    // VTK Legacy
    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );
    base::io::vtk::writePointData( vtkWriter, mesh, temperature, "temperature" );
    vtk.close();
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input.dat \n\n";
        return -1;
    }

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double kappa, density, numSteps, stepSize;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "kappa",            kappa );
        prop.registerPropertiesVar( "density",          density );
        prop.registerPropertiesVar( "numSteps",         numSteps );
        prop.registerPropertiesVar( "stepSize",         stepSize );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        prop.readValues( inp );
        inp.close( );

        // Make sure all variables have been found
        if ( not prop.isEverythingRead() ) {
            prop.writeUnread( std::cerr );
            VERIFY_MSG( false, "Could not find above variables" );
        }
    }
    
    const std::string baseName = base::io::baseName( meshFile, ".smf" );
    
    //--------------------------------------------------------------------------
    const unsigned    geomDeg  = 1;   // degree of geometry approximation
    const unsigned    tiOrder  = 2;   // order of time integrator
    const base::Shape shape    = base::QUAD; // shape of the element

    // choose a time stepping method
    //typedef  base::time::BDF<tiOrder> MSM;
    typedef base::time::AdamsMoulton<tiOrder> MSM;

    // time stepping method determines the history size
    const unsigned nHist = MSM::numSteps;
    
    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned dim = Mesh::Node::dim;

    // create a mesh from file
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature 
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // Finite element basis
    typedef base::fe::IsoparametricBasis<Mesh> FEBasis;
    const unsigned    doFSizeT = 1; // temperature
    const unsigned    doFSizeV = dim; // velocity

    // DOF handling temperature
    typedef base::Field<FEBasis,doFSizeT,nHist>        Temperature;
    Temperature temperature;
    base::dof::generate<FEBasis>( mesh, temperature );

    // Field for velocity variable
    typedef base::Field<FEBasis,doFSizeV>              Velocity;
    Velocity velocity;
    base::dof::generate<FEBasis>( mesh, velocity );

    // set field according to a function
    base::dof::setField( mesh, velocity,
                         boost::bind( &velocityFun<dim,
                                                   Velocity::DegreeOfFreedom>, _1, _2 ) );


    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // constrain the boundary degrees of freedom
    base::dof::constrainBoundary<FEBasis>( meshBoundary.boundaryBegin(),
                                           meshBoundary.boundaryEnd(),
                                           mesh, temperature,
                                           boost::bind( &dirichletBC<dim,
                                                                     Temperature::DegreeOfFreedom>,
                                                        _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( temperature.doFsBegin(),
                                            temperature.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs
              << ", number of time steps " << numSteps
              << std::endl;

    typedef base::asmb::FieldBinder<Mesh,Temperature,Velocity> FieldBinder;
    FieldBinder fieldBinder( mesh, temperature, velocity );
    typedef FieldBinder::TupleBinder<1,1,2>::Type FieldTupleBinder;

    // choose a material behaviour
    typedef mat::thermal::IsotropicConstant Material;
    Material material( kappa );

    // Object for static heat transfer
    typedef heat::Static<Material,FieldTupleBinder::Tuple> StaticHeat;
    StaticHeat staticHeat( material );

    // Convection term
    typedef heat::Convection<FieldTupleBinder::Tuple> Convection;
    Convection convection( density );

    // write initial state
    writeVTKFile( baseName, 0, mesh, temperature );

    // Monitor of solution
    const std::size_t numElements = std::distance( mesh.elementsBegin(),
                                                   mesh.elementsEnd() );
    const std::size_t N = static_cast<std::size_t>( std::sqrt( numElements ) );
    const std::size_t midElement = (N * N + N) / 2;
    base::post::Monitor<Mesh::Element,Temperature::Element>
        monitorSol( mesh.elementPtr(        midElement ),
                    temperature.elementPtr( midElement ),
                    base::constantVector<dim>( 0. ) );

    std::cout << "# " << ( monitorSol.location() ).transpose() << std::endl;

    //--------------------------------------------------------------------------
    // Time loop
    for ( unsigned step = 0; step < numSteps; step ++ ) {

        //std::cout << "Time step " << step << std::endl;
        const double time = step * stepSize;
        std::cout << time << " ";
        monitorSol.solution( std::cout );
        monitorSol.gradient( std::cout );
        std::cout << std::endl;
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        // Compute system matrix from Laplacian
        base::asmb::stiffnessMatrixComputation<FieldTupleBinder>( quadrature, solver,
                                                                  fieldBinder,
                                                                  staticHeat );

        // Compute system matrix from Convection
        base::asmb::stiffnessMatrixComputation<FieldTupleBinder>( quadrature, solver, 
                                                                  fieldBinder,
                                                                  convection );

        // compute inertia terms, d/dt, due to time integration
        base::time::computeInertiaTerms<FieldTupleBinder,MSM>( quadrature, solver,
                                                               fieldBinder, stepSize, step,
                                                               density );

        // compute history of residual forces due to time integration
        base::time::computeResidualForceHistory<FieldTupleBinder,MSM>( staticHeat,
                                                                       quadrature, solver,
                                                                       fieldBinder, step );

        // Finalise assembly
        solver.finishAssembly();

        // Solve a possibly non-symmetric system
        solver.superLUSolve();

        // distribute results back to dofs
        {
            base::dof::setDoFsFromSolver( solver, temperature );

            // push history
            std::for_each( temperature.doFsBegin(), temperature.doFsEnd(),
                           boost::bind( &Temperature::DegreeOfFreedom::pushHistory, _1 ) );
        }

        writeVTKFile( baseName, step+1, mesh, temperature );

    } // end time loop

    return 0;
}
