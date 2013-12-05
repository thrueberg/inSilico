//------------------------------------------------------------------------------
//  Solve the diffusion equation in one dimension
//
//  \dot{u} - D u'' = r
//
//  where D has the distribution
//
//       D1               D2               D1
//   +----------+---------------------+----------+
//  x=0        x=x1                  x=x2       x=x3
//
//  with x1 = L1, x2 = L1 + L2 and x3 = L1 +L2 +L1
//
//  Moreover, the concentration u has the initial distribution
//
//    -----------
//    |         |
//    |  u=u0   |          u=0             u=0
//   +----------+---------------------+----------+
//
//  and boundary conditions u' = 0 at x=0 and x=x3. Between
//  x=x1 and x=x2 a concentration sink r2 is applied such that
//  one has the source terms
//
//       r=0               r=r2            r=0
//   +----------+---------------------+----------+
//


//------------------------------------------------------------------------------
// system includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/bind.hpp>
// mesh related
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
// input/output
#include <base/io/Format.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
// quadrature
#include <base/Quadrature.hpp>
// FE basis
#include <base/fe/Basis.hpp>
// Field and degrees of freedom
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/setField.hpp>
// integral kernel
#include <heat/Laplace.hpp>
// system soler
#include <base/solver/Eigen3.hpp>
// assembly
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
// time integration
#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>
//
#include <base/kernel/FieldIntegral.hpp>

//------------------------------------------------------------------------------
// Generate a mesh with Ni elements per domain with length Li
template<typename MESH>
void generateMesh( MESH& mesh, const double L1, const double L2,
                   const unsigned N1, const unsigned N2 )
{
    static const unsigned degree = MESH::Element::GeomFun::degree;
    
    std::stringstream buffer;
    // write header
    buffer << "! elementShape line \n"
           << "! elementNumPoints  " << degree+1 << "\n";

    // total number of elements
    const unsigned numElements = N1 + N2 + N1;
    // total number of nodes
    const unsigned numNodes = 1 + degree * numElements;

    buffer << numNodes << "  " << numElements << "\n";

    // element lengths
    const double h1 = L1 / static_cast<double>( N1 );
    const double h2 = L2 / static_cast<double>( N2 );

    // generate coordinates
    double x = 0;
    for ( unsigned n = 0; n < N1; n++ ) {
        for ( unsigned d = 0; d < degree; d++ ) {
            buffer << x << "  0.  0. \n";
            x += h1;
        }
    }

    for ( unsigned n = 0; n < N2; n++ ) {
        for ( unsigned d = 0; d < degree; d++ ) {
            buffer << x << "  0.  0. \n";
            x += h2;
        }
    }

    for ( unsigned n = 0; n < N1; n++ ) {
        for ( unsigned d = 0; d < degree; d++ ) {
            buffer << x << "  0.  0. \n";
            x += h1;
        }
    }

    // last node
    buffer << x << "  0.  0. \n";

    // generate connectivity
    unsigned v = 0;
    for ( unsigned n = 0; n < numElements; n++ ) {
        for ( unsigned d = 0; d < degree+1; d++ ) {
            buffer << v+d << " ";
        }
        buffer << "\n";
        v += degree;
    }

    base::io::smf::readMesh( buffer, mesh );

    return;
}

//--------------------------------------------------------------------------
typedef base::Vector<1>::Type VecDim;

template<typename DOF>
void initialState( const VecDim& x, DOF* doFPtr,
                   const double L1, const double L2,
                   const double u0 )
{
    const double value =  (x[0] <= L1 ? u0 : 0.0 );
    //const double value = ( (x[0] >= L1 and x[0] <= L1+L2) ? u0 : 0.0);
    doFPtr -> setValue( 0, value );
    doFPtr -> pushHistory();
}


//--------------------------------------------------------------------------
template<typename ELEMENT, typename FIELD>
base::Vector<1>::Type massProduction( const ELEMENT* geomEp,
                                      const typename ELEMENT::GeomFun::VecDim& xi,
                                      const FIELD& field,
                                      const double rate )
{
    const typename FIELD::Element* fieldEp = field.elementPtr( geomEp -> getID() );
    
    const base::Vector<1>::Type A = base::post::evaluateField( geomEp, fieldEp, xi );

    return -rate * A;
}

//--------------------------------------------------------------------------
// output to a VTK file and add a line of concentrations to the data file
template<typename MESH, typename FIELD>
void writeData( const MESH& mesh, const FIELD& field, 
                const std::string baseName,
                const unsigned step )
{
    typedef base::Vector<1>::Type VecDoF;
    std::vector<VecDoF> nodalValues;
    base::post::evaluateAtNodes( mesh, field, nodalValues );

    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk(  vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );
    vtkWriter.writePointData( nodalValues.begin(), nodalValues.end(), "u" );
    vtk.close();

    const std::string outFile = baseName + ".dat";
    std::ofstream out( outFile.c_str(), std::ofstream::app );
    for ( std::size_t n = 0; n < nodalValues.size(); n++ )
        out << nodalValues[n][0] << "  ";
    out << "\n";
    
    return;
}

//------------------------------------------------------------------------------
template<typename ELEMENT>
double diffusionConstant( const ELEMENT* geomEp,
                          const typename ELEMENT::GeomFun::VecDim& xi,
                          const double L1, const double L2,
                          const double D1, const double D2 )
{
    const typename ELEMENT::Node::VecDim x =
        base::Geometry<ELEMENT>()( geomEp, xi );

    const bool domain2 = (x[0] >= L1 and x[0] <= L1+L2);

    return (domain2 ? D2 : D1);
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    const unsigned    geomDeg           = 1;
    const unsigned    fieldDeg          = 2;
    const base::Shape shape             = base::LINE;
    const unsigned    doFSize           = 1;
    const unsigned    kernelDegEstimate = 3;
    const unsigned    tiOrder           = 2;

    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input \n\n";
        return -1;
    }

    // first command line argument is the input data file
    const std::string inputFile  = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    double L1, L2, stepSize, u0, D1, D2, rate;
    unsigned N1, N2, numSteps;
    std::string outFile;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "L1",          L1 );
        prop.registerPropertiesVar( "L2",          L2 );
        prop.registerPropertiesVar( "N1",          N1 );
        prop.registerPropertiesVar( "N2",          N2 );

        prop.registerPropertiesVar( "numSteps",    numSteps );
        prop.registerPropertiesVar( "stepSize",    stepSize );

        prop.registerPropertiesVar( "u0",          u0 );

        prop.registerPropertiesVar( "D1",          D1 );
        prop.registerPropertiesVar( "D2",          D2 );

        prop.registerPropertiesVar( "rate",        rate );

        prop.registerPropertiesVar( "outFile",     outFile );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck(inp) , "Input error" );
        inp.close( );
    }

    //--------------------------------------------------------------------------
    // create a mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    Mesh mesh;
    generateMesh( mesh, L1, L2, N1, N2 );

    // Quadrature and surface quadrature
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // time integration
    typedef base::time::BDF<tiOrder> MSM;
    //typedef base::time::AdamsMoulton<tiOrder> MSM;
    const unsigned nHist = MSM::numSteps;

    // DOF handling
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>     Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field concentration;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, concentration );

    // set initial condition
    base::dof::setField( mesh, concentration, boost::bind( &initialState<DoF>,
                                                           _1, _2, L1, L2, u0 ) );
    
    // boundary conditions
#if 0 // this piece of code is not active -> do-nothing neumann condition is applied
    Field::DoFPtrIter doFIter = concentration.doFsBegin();
    (*doFIter) -> constrainValue( 0, 0.0 );
    doFIter = concentration.doFsEnd(); doFIter--;
    (*doFIter) -> constrainValue( 0, 0.0 );
#endif
    
    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( concentration.doFsBegin(),
                                            concentration.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs
              << "# Time    total mass \n" << std::endl;


    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, concentration );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // write header to data file
    {
        std::string fileName = outFile + ".dat";
        std::ofstream out( fileName.c_str() );
        out << "# x=";
        Mesh::NodePtrConstIter nIter = mesh.nodesBegin();
        Mesh::NodePtrConstIter nEnd  = mesh.nodesEnd();
        for ( ; nIter != nEnd; ++nIter ) {
            double x[1];
            (*nIter) -> getX( &(x[0]) );
            out << x[0] << " ";
        }
        out << "\n";
    }

    // write initial state
    writeData( mesh, concentration, outFile, 0 );

    {
        std::cout << "0.  ";
        base::Vector<1>::Type totalMass = base::constantVector<1>( 0. );
        base::asmb::simplyIntegrate<FTB>( quadrature, totalMass, fieldBinder,
                                          base::kernel::FieldIntegral<FTB::Tuple>() );
        std::cout << totalMass[0] << std::endl;
    }

    typedef heat::Laplace<FTB::Tuple> Laplace;
    Laplace laplace( 1.0 );
    boost::function< double( const Mesh::Element *,
                             const Mesh::Element::GeomFun::VecDim & ) >
        diffusionFun = boost::bind( &diffusionConstant<Mesh::Element>, _1, _2,
                                    L1, L2, D1, D2);
    laplace.setConductivityFunction( diffusionFun );

    //--------------------------------------------------------------------------
    // Time loop
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < numSteps; step++ ) {

        std::cout << (step+1)*stepSize << "  ";
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        // compute stiffness matrix
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

        //
        base::asmb::bodyForceComputation2<FTB>(
            quadrature, solver, fieldBinder,
            boost::bind( &massProduction<Mesh::Element,Field>,
                         _1, _2, boost::ref( concentration ), rate ) );

        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.choleskySolve();

        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, concentration );


        // pass to history 
        std::for_each( concentration.doFsBegin(), concentration.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

        // write output
        writeData( mesh, concentration, outFile, step+1 );

        base::Vector<1>::Type totalMass = base::constantVector<1>( 0. );
        base::asmb::simplyIntegrate<FTB>( quadrature, totalMass, fieldBinder,
                                          base::kernel::FieldIntegral<FTB::Tuple>() );

        std::cout << totalMass[0] << std::endl;
    }

    return 0;
}
