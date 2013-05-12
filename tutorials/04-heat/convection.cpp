// system header
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>

#include <base/fe/Basis.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Element.hpp>
#include <base/dof/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/setField.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

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

    for ( unsigned d = 0; d < DOF::nHist; d++ ) doFPtr -> pushHistory();
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
    const unsigned    fieldDeg = 1;   // degree of field approximation
    const unsigned    tiOrder  = 2;   // order of time integrator
    const base::Shape shape    = base::QUAD; // shape of the element

    // choose a time stepping method
    //typedef  base::time::BDF<tiOrder> MSM;
    typedef base::time::AdamsMoulton<tiOrder> MSM;

    // time stepping method determines the history size
    const unsigned nHist = MSM::numSteps;
    
    //--------------------------------------------------------------------------
    // shape implies spatial dimension (not necessarily)
    const unsigned    dim     = base::ShapeDim<shape>::value;
    typedef base::mesh::Node<dim>                 Node;     // geometry node
    typedef base::LagrangeShapeFun<geomDeg,shape> SFun;     // geometry shape fun
    typedef base::mesh::Element<Node,SFun>        Element;  // geometry element
    typedef base::mesh::Unstructured<Element>     Mesh;     // mesh

    // create a mesh from file
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::Reader<Mesh> smfReader;
        smfReader( mesh, smf ); 
        smf.close();
    }

    // Quadrature 
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // Finite element basis
    const unsigned    doFSizeT = 1; // temperature
    const unsigned    doFSizeV = dim; // velocity
    typedef base::fe::Basis<shape,fieldDeg> FEBasis;

    // DOF handling temperature
    typedef base::LagrangeShapeFun<fieldDeg,shape>     FieldFun;
    typedef base::dof::DegreeOfFreedom<doFSizeT,nHist> DoFT;
    typedef base::dof::Element<DoFT,FEBasis::FEFun>    FieldElementT;
    typedef base::dof::Field<FieldElementT>            Temperature;
    Temperature temperature;
    base::dof::generate<FEBasis>( mesh, temperature );

    // Field for velocity variable
    typedef base::dof::DegreeOfFreedom<doFSizeV,nHist> DoFV;
    typedef base::dof::Element<DoFV,FEBasis::FEFun>    FieldElementV;
    typedef base::dof::Field<FieldElementV>            Velocity;
    Velocity velocity;
    base::dof::generate<FEBasis>( mesh, velocity );

    // set field according to a function
    base::dof::setField( mesh, velocity,
                         boost::bind( &velocityFun<dim,DoFV>, _1, _2 ) );


    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // constrain the boundary degrees of freedom
    base::dof::constrainBoundary<FEBasis>( meshBoundary.boundaryBegin(),
                                           meshBoundary.boundaryEnd(),
                                           mesh, temperature,
                                           boost::bind( &dirichletBC<dim,DoFT>, _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( temperature.doFsBegin(),
                                            temperature.doFsEnd() );
    std::cout << " Number of dofs " << numDofs
              << ", number of time steps " << numSteps
              << std::endl;

    typedef base::asmb::FieldBinder<Mesh,Temperature,Temperature,Velocity> FieldBinder;
    FieldBinder fieldBinder( mesh, temperature, temperature, velocity );
    typedef FieldBinder::ElementPtrTuple FieldTuple;

    // choose a material behaviour
    typedef mat::thermal::IsotropicConstant Material;
    Material material( kappa );

    // Object for static heat transfer
    typedef heat::Static<Material,FieldTuple> StaticHeat;
    StaticHeat staticHeat( material );

    // Convection term
    typedef heat::Convection<FieldTuple> Convection;
    Convection convection( density );

    // write initial state
    writeVTKFile( baseName, 0, mesh, temperature );

    //--------------------------------------------------------------------------
    // Time loop
    for ( unsigned step = 0; step < numSteps; step ++ ) {

        std::cout << "Time step " << step << std::endl;
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );


        // Compute system matrix from Laplacian
        base::asmb::stiffnessMatrixComputation( quadrature, solver,
                                                fieldBinder,
                                                staticHeat );

        // Compute system matrix from Convection
        base::asmb::stiffnessMatrixComputation( quadrature, solver, 
                                                fieldBinder,
                                                convection );

        // Compute system matrix and RHS from time stepping
        typedef base::time::ReactionTerms<Quadrature,Solver,MSM,FieldTuple> RT;
        RT rt( density, quadrature, solver, stepSize, step, 0 );

        // iterate over the fields
        std::for_each( fieldBinder.elementsBegin(), fieldBinder.elementsEnd(), rt );

        // Compute RHS terms from history of forces
        base::time::ResidualForceHistory<StaticHeat,Quadrature,Solver,MSM>
            rfh( staticHeat, quadrature, solver, step );

        // iterate over the fields
        std::for_each( fieldBinder.elementsBegin(), fieldBinder.elementsEnd(),
                       rfh );

        // Finalise assembly
        solver.finishAssembly();

        // Solve a possibly non-symmetric system
        solver.superLUSolve();

        // distribute results back to dofs
        {
            base::dof::Distribute<DoFT,Solver,base::dof::SET> distributeDoF( solver );
            std::for_each( temperature.doFsBegin(),
                           temperature.doFsEnd(), distributeDoF );

            // push history
            std::for_each( temperature.doFsBegin(), temperature.doFsEnd(),
                           boost::bind( &DoFT::pushHistory, _1 ) );
        }

        writeVTKFile( baseName, step+1, mesh, temperature );

    } // end time loop

    return 0;
}
