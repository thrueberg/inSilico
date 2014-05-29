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
#include <base/dof/location.hpp>
// integral kernel
#include <heat/Laplace.hpp>
#include <heat/Convection.hpp>
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

#include "FieldIntegralPart.hpp" 

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
void initialConcentration( const VecDim& x, DOF* doFPtr,
                           const double L1, const double L2,
                           const double u0 )
{
    const double value =  (x[0] <= L1 ? u0 : 0.0 );
    doFPtr -> setValue( 0, value );
    doFPtr -> pushHistory();
}

//------------------------------------------------------------------------------
// Advection velocity according to Darcy's law
template<typename DOF>
void advectionVelocity( const VecDim& x, DOF* doFPtr,
                        const double L1, const double L2,
                        const double p1, const double p2,
                        const double perm )
{
    const bool domain2 = ( x[0] > L1 ) and ( x[0] < L1+L2 );

    typename base::Vector<DOF::size>::Type V = base::constantVector<DOF::size>( 0. );

    if ( domain2 ) V[0] = perm * (p1 - p2);

    for ( unsigned d = 0; d < DOF::size; d++ )
        doFPtr -> setValue( d, V[d] );

    doFPtr -> pushHistory();
}


//--------------------------------------------------------------------------
// Mass sink in domain 2
template<typename ELEMENT, typename FIELD>
base::Vector<1>::Type massProduction( const ELEMENT* geomEp,
                                      const typename ELEMENT::GeomFun::VecDim& xi,
                                      const FIELD& field,
                                      const double rate1, const double rate2, 
                                      const double L1,    const double L2 )
{
    typename base::GeomTraits<ELEMENT>::GlobalVecDim x
        = base::Geometry<ELEMENT>()( geomEp, xi );

    const bool inDomain2 = ( (x[0] > L1) and (x[0] < L1+L2) );
    
    const typename FIELD::Element* fieldEp = field.elementPtr( geomEp -> getID() );
    
    const base::Vector<1>::Type A = base::post::evaluateField( geomEp, fieldEp, xi );

    return (inDomain2 ? -rate2 : -rate1) * A;
}

//--------------------------------------------------------------------------
template<typename FIELD, typename MESH>
void collectProduct( const FIELD& concentration,
                     FIELD& product,
                     const MESH& mesh,
                     const std::vector<std::pair<std::size_t,VecDim> >&
                     doFLocation, 
                     const double rate1, const double rate2, 
                     const double stepSize,
                     const double L1, const double L2 )
{
    typename FIELD::DoFPtrConstIter cIter = concentration.doFsBegin();
    typename FIELD::DoFPtrConstIter cEnd  = concentration.doFsEnd();
    typename FIELD::DoFPtrConstIter outIter = product.doFsBegin();
    for ( ; cIter != cEnd; ++cIter, ++outIter ) {

        const std::pair<std::size_t,VecDim> aux = doFLocation[ (*cIter) -> getID() ];

        typename MESH::Element* geomEp = mesh.elementPtr( aux.first );
        typename MESH::Node::VecDim x  =
            base::Geometry<typename MESH::Element>()( geomEp, aux.second );

        const bool inDomain2 = (x[0] > L1) and (x[0] < L1+L2);

        const double rate = ( inDomain2 ? rate2 : rate1 );

        for ( unsigned d = 0; d < FIELD::DegreeOfFreedom::size; d++ ) {
            const double value    = rate * stepSize * ( (*cIter) -> getValue( d ) );
            const double previous = (*outIter) -> getValue(d);
            (*outIter) -> setValue( d, previous + value );
        }
        
    }
    return;
}


//--------------------------------------------------------------------------
// output to a VTK file and add a line of concentrations to the data file
template<typename MESH, typename FIELD, typename QUADRATURE>
void writeData( const MESH& mesh,
                const FIELD& concentration,
                const FIELD& product,
                const QUADRATURE& quadrature, 
                const std::string baseName,
                const unsigned step,
                const double stepSize,
                const double x1, const double x2)
{
    // evaluate the concentration and the product at the FE nodes
    typedef base::Vector<1>::Type VecDoF;
    std::vector<VecDoF> nodalValuesC, nodalValuesR;
    base::post::evaluateFieldAtNodes( mesh, concentration, nodalValuesC );
    base::post::evaluateFieldAtNodes( mesh, product,       nodalValuesR );

    // write a VTK file
    {
        const std::string vtkFile =
            baseName + "." + base::io::leadingZeros( step ) + ".vtk";
        std::ofstream vtk(  vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        vtkWriter.writePointData( nodalValuesC.begin(), nodalValuesC.end(), "c" );
        vtkWriter.writePointData( nodalValuesR.begin(), nodalValuesR.end(), "r" );
        vtk.close();

    }

    // Write nodal data to ASCII files
    {
        const std::string outFileC = baseName + ".c.dat";
        const std::string outFileR = baseName + ".r.dat";
        std::ofstream outC( outFileC.c_str(), std::ofstream::app );
        std::ofstream outR( outFileR.c_str(), std::ofstream::app );
        for ( std::size_t n = 0; n < nodalValuesC.size(); n++ ) {
            outC << nodalValuesC[n][0] << "  ";
            outR << nodalValuesR[n][0] << "  ";
        }
        outC << "\n"; outC.close();
        outR << "\n"; outR.close();
    }

    // compute the total mass of the concentration and of the product
    {
        typedef base::asmb::FieldBinder<const MESH,const FIELD> FieldBinder;
        typedef typename FieldBinder::template TupleBinder<1,1>::Type FTB;
        FieldBinder fb1( mesh, concentration ); 
        
        base::Vector<1>::Type totalMass = base::constantVector<1>( 0. );
        base::asmb::simplyIntegrate<FTB>( quadrature, totalMass, fb1,
                                          base::kernel::FieldIntegral<typename FTB::Tuple>() );

        FieldBinder fb2( mesh, product );
        base::Vector<1>::Type totalLoss = base::constantVector<1>( 0. );
        base::asmb::simplyIntegrate<FTB>( quadrature, totalLoss, fb2,
                                          base::kernel::FieldIntegral<typename FTB::Tuple>() );


        typedef FieldIntegralPart<typename FTB::Tuple> FIP;
        FIP fip( x1, x2 );
        
        typedef typename FIP::ResultPair ResultPair;
        ResultPair resultPair;
        resultPair.first  = base::constantVector<3>( 0. );
        resultPair.second = base::constantVector<3>( 0. );
        
        base::asmb::simplyIntegrate<FTB>( quadrature, resultPair, fb1, fip );
                                          
        
        std::cout << (step+1)*stepSize << "  "
                  << totalMass[0] << "  "  << totalLoss[0]
                  << "  " << resultPair.second[0] / resultPair.first[0]
                  << "  " << resultPair.second[1] / resultPair.first[1]
                  << "  " << resultPair.second[2] / resultPair.first[2]
                  << std::endl;
    }
    
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
    const unsigned    doFSizeC          = 1;
    const unsigned    kernelDegEstimate = fieldDeg + fieldDeg;
    const unsigned    tiOrder           = 2;

    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " input \n\n";
        return -1;
    }

    // first command line argument is the input data file
    const std::string inputFile  = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    double L1, L2, stepSize, u0, D1, D2, rate1, rate2, p1, p2, perm;
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

        prop.registerPropertiesVar( "rate1",       rate1 );
        prop.registerPropertiesVar( "rate2",       rate2 );

        prop.registerPropertiesVar( "p1",          p1 );
        prop.registerPropertiesVar( "p2",          p2 );
        prop.registerPropertiesVar( "perm",        perm );

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

    // DOF handling (concentration)
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSizeC,nHist>    Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field concentration, product;

    // DOF handling (advection velocity)
    const unsigned doFSizeV = Mesh::Node::dim;
    typedef base::Field<FEBasis,doFSizeV,nHist>    FieldV;
    typedef FieldV::DegreeOfFreedom                DoFV;
    FieldV advection;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, concentration );
    base::dof::generate<FEBasis>( mesh, product   );
    base::dof::generate<FEBasis>( mesh, advection );

    // find geometry association for the dofs
    std::vector<std::pair<std::size_t,VecDim> > doFLocation;
    base::dof::associateLocation( concentration, doFLocation );

    // set initial conditions
    base::dof::setField( mesh, concentration, boost::bind( &initialConcentration<DoF>,
                                                           _1, _2, L1, L2, u0 ) );

    base::dof::setField( mesh, advection,     boost::bind( &advectionVelocity<DoFV>,
                                                           _1, _2, L1, L2, p1, p2, perm ) );
    
    // boundary conditions
#if 0 // this piece of code is not active -> do-nothing neumann condition is applied
    Field::DoFPtrIter doFIter = concentration.doFsBegin();
    (*doFIter) -> constrainValue( 0, 1.0 );
    doFIter = concentration.doFsEnd(); doFIter--;
    (*doFIter) -> constrainValue( 0, 1.0 );
#endif
    
    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( concentration.doFsBegin(),
                                            concentration.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs << "\n "
              << "# Time    total mass    total loss   average1  average2  average3 "
              << std::endl;

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field,FieldV> FieldBinder;
    FieldBinder fieldBinder( mesh, concentration, advection );
    typedef FieldBinder::TupleBinder<1,1,2>::Type FTB;

    // write header to data file
    {
        std::string fileName = outFile + ".c.dat";
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
    writeData( mesh, concentration, product, quadrature, outFile, 0, stepSize, L1, L1+L2 );

    typedef heat::Laplace<FTB::Tuple> Laplace;
    Laplace laplace( 1.0 );
    boost::function< double( const Mesh::Element *,
                             const Mesh::Element::GeomFun::VecDim & ) >
        diffusionFun = boost::bind( &diffusionConstant<Mesh::Element>, _1, _2,
                                    L1, L2, D1, D2);
    laplace.setConductivityFunction( diffusionFun );
    
    // Convection term
    typedef heat::Convection<FTB::Tuple> Convection;
    Convection convection( 1.0 );


    //--------------------------------------------------------------------------
    // Time loop
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < numSteps; step++ ) {

        // collect the concentration of the product in an extra field
        collectProduct( concentration, product, mesh, doFLocation,
                        rate1, rate2, stepSize, L1, L2 );
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        // compute the RHS due to a mass production term
        base::asmb::bodyForceComputation2<FTB>(
            quadrature, solver, fieldBinder,
            boost::bind( &massProduction<Mesh::Element,Field>,
                         _1, _2, boost::ref( concentration ),
                         rate1, rate2, L1, L2 ) );
        
        // compute stiffness matrix
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                     fieldBinder, laplace, false );

        // Compute system matrix from Convection
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver, 
                                                     fieldBinder, convection, false );


        // compute inertia terms, d/dt, due to time integration
        base::time::computeInertiaTerms<FTB,MSM>( quadrature, solver,
                                                  fieldBinder, stepSize, step,
                                                  1.0, false );

        // compute history of residual forces due to time integration
        base::time::computeResidualForceHistory<FTB,MSM>( laplace, 
                                                          quadrature, solver,
                                                          fieldBinder, step );
        
        // compute history of residual forces due to time integration
        base::time::computeResidualForceHistory<FTB,MSM>( convection, 
                                                          quadrature, solver,
                                                          fieldBinder, step );

        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.superLUSolve();

        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, concentration );


        // pass to history 
        std::for_each( concentration.doFsBegin(), concentration.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );

        // write output
        writeData( mesh, concentration, product, quadrature,
                   outFile, step+1, stepSize, L1, L1+L2 );

    }

    return 0;
}
