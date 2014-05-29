//------------------------------------------------------------------------------
//  Solve the diffusion equation in two dimensions
//
//  \dot{u} - D (u_{,xx} + u_{,yy}) = r
//
//  where D has the distribution
//
//       D1         D2          D3         D1
//   +----------+----------+----------+----------+
//  x=x_left   x=x1    (x1+x2)/2     x=x2       x=x_right
//
//  with x1 = - 5/13 bgel1 and x2 = + 8/13 bgel1
//
//  Moreover, the concentration u has the initial distribution
//
//    ---------
//    |       |
//    | u=u0  |           u=0             u=0
//   +--------+--+---------------------+----------+
//            x=x0
//
//  and boundary conditions u' = 0 at all boundaries. Between
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
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
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
#include <base/dof/constrainBoundary.hpp>
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

#include "../../manu/diffusion/FieldIntegralPart.hpp"

//--------------------------------------------------------------------------
const double tolerance = 1.e-10;

static const unsigned DIM=2;

typedef base::Vector<DIM>::Type VecDim;


//------------------------------------------------------------------------------
// Set the initial value to 'u0' in the left domain (x < x1)
template<typename DOF>
void initialConcentration( const VecDim& x, DOF* doFPtr,
                           const double x1, const double x2,
                           const double u0 )
{
    const double value =  (x[0] <= x1 ? u0 : 0.0 );
    doFPtr -> setValue( 0, value );
    doFPtr -> pushHistory();
}

//------------------------------------------------------------------------------
// Advection velocity according to Darcy's law for the region x1 < x < x2
template<typename DOF>
void advectionVelocity( const VecDim& x, DOF* doFPtr,
                        const double x1, const double x2,
                        const double p1, const double p2,
                        const double perm )
{
    const bool domain2 = ( x[0] > x1 ) and ( x[0] < x2 );

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
                                      const VecDim& xi,
                                      const FIELD& field,
                                      const double rate1, const double rate2, 
                                      const double x1,    const double x2 )
{
    VecDim x = base::Geometry<ELEMENT>()( geomEp, xi );

    const bool inDomain2 = ( (x[0] > x1) and (x[0] < x2) );
    
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
                     const double x1, const double x2 )
{
    typename FIELD::DoFPtrConstIter cIter = concentration.doFsBegin();
    typename FIELD::DoFPtrConstIter cEnd  = concentration.doFsEnd();
    typename FIELD::DoFPtrConstIter outIter = product.doFsBegin();
    for ( ; cIter != cEnd; ++cIter, ++outIter ) {

        const std::pair<std::size_t,VecDim> aux = doFLocation[ (*cIter) -> getID() ];

        typename MESH::Element* geomEp = mesh.elementPtr( aux.first );
        VecDim x  = base::Geometry<typename MESH::Element>()( geomEp, aux.second );

        const bool inDomain2 = (x[0] > x1) and (x[0] < x2);

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
                const double x1,
                const double x2,
                const unsigned writeEvery )
{
    // evaluate the concentration and the product at the FE nodes
    typedef base::Vector<1>::Type VecDoF;
    std::vector<VecDoF> nodalValuesC, nodalValuesR;
    base::post::evaluateFieldAtNodes( mesh, concentration, nodalValuesC );
    base::post::evaluateFieldAtNodes( mesh, product,       nodalValuesR );

    // write a VTK file
    if ( step % writeEvery == 0 ) {

        const std::string vtkFile =
            baseName + "." + base::io::leadingZeros( step, 5 ) + ".vtk";
        std::ofstream vtk(  vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        vtkWriter.writePointData( nodalValuesC.begin(), nodalValuesC.end(), "c" );
        vtkWriter.writePointData( nodalValuesR.begin(), nodalValuesR.end(), "r" );
        vtk.close();

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
                                          
        
        std::cout << step*stepSize << "  "
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
                          const VecDim& xi,
                          const double x1, const double x2,
                          const double D1, const double D2,
                          const double D3 )
{
    const VecDim x = base::Geometry<ELEMENT>()( geomEp, xi );

    const double mid = 0.5 * (x1 + x2);

    const bool domain2 = (x[0] >= x1 and x[0] <= mid);
    const bool domain3 = (x[0] > mid and x[0] <= x2 );

    return (domain2 ? D2 : (domain3 ? D3 : D1) );
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    const unsigned    geomDeg           = 1;
    const unsigned    fieldDeg          = 1;
    const base::Shape shape             = base::TRI;
    const unsigned    doFSizeC          = 1;
    const unsigned    kernelDegEstimate = fieldDeg + fieldDeg;
    const unsigned    tiOrder           = 1;

    if ( argc != 3 ) {
        std::cout << "Usage:  " << argv[0] << " meshFile  input \n\n";
        return -1;
    }

    // first command line argument is the input data file
    const std::string meshFile   = boost::lexical_cast<std::string>( argv[1] );
    const std::string inputFile  = boost::lexical_cast<std::string>( argv[2] );
    
    // read from input file
    double bgel1, bcan1, stepSize, u0, D1, D2, D3, rate1, rate2, p1, p2, perm;
    unsigned numSteps, writeEvery;
    std::string outFile;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "bgel1",       bgel1 );
        prop.registerPropertiesVar( "bcan1",       bcan1 );

        prop.registerPropertiesVar( "numSteps",    numSteps );
        prop.registerPropertiesVar( "stepSize",    stepSize );

        prop.registerPropertiesVar( "u0",          u0 );

        prop.registerPropertiesVar( "D1",          D1 );
        prop.registerPropertiesVar( "D2",          D2 );
        prop.registerPropertiesVar( "D3",          D3 );

        prop.registerPropertiesVar( "rate1",       rate1 );
        prop.registerPropertiesVar( "rate2",       rate2 );

        prop.registerPropertiesVar( "p1",          p1 );
        prop.registerPropertiesVar( "p2",          p2 );
        prop.registerPropertiesVar( "perm",        perm );

        prop.registerPropertiesVar( "outFile",     outFile );
        prop.registerPropertiesVar( "writeEvery",  writeEvery );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck(inp) , "Input error" );
        inp.close( );
    }

    // x_1 coordinates of the domain seperating lines
    const double x1 = -5./13. * bgel1;
    const double x2 =  8./13. * bgel1;

    // x_0 rightmost point in the post (approx.)
    const double x0 = x1 - 1443./460.*bcan1;

    // with advection
    const bool withAdvection =
        ( std::abs( p1 - p2 ) < tolerance ) or
        ( std::abs( perm )    < tolerance );

    //--------------------------------------------------------------------------
    // create a mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature and surface quadrature
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // time integration
    typedef base::time::BDF<tiOrder> MSM;
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
                                                           _1, _2, x0, x2, u0 ) );

    base::dof::setField( mesh, advection,     boost::bind( &advectionVelocity<DoFV>,
                                                           _1, _2, x1, x2, p1, p2, perm ) );

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

    // write initial state
    writeData( mesh, concentration, product, quadrature, outFile, 0, stepSize, x1, x2,
               writeEvery );

    typedef heat::Laplace<FTB::Tuple> Laplace;
    Laplace laplace( 1.0 );
    boost::function< double( const Mesh::Element *,
                             const Mesh::Element::GeomFun::VecDim & ) >
        diffusionFun = boost::bind( &diffusionConstant<Mesh::Element>, _1, _2,
                                    x1, x2, D1, D2, D3);
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
                        rate1, rate2, stepSize, x1, x2 );
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        // compute the RHS due to a mass production term
        base::asmb::bodyForceComputation2<FTB>(
            quadrature, solver, fieldBinder,
            boost::bind( &massProduction<Mesh::Element,Field>,
                         _1, _2, boost::ref( concentration ),
                         rate1, rate2, x1, x2 ) );
        
        // compute stiffness matrix
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                     fieldBinder, laplace );

        // compute inertia terms, d/dt, due to time integration
        const bool incremental = false;
        base::time::computeInertiaTerms<FTB,MSM>( quadrature, solver,
                                                  fieldBinder, stepSize, step,
                                                  1.0, incremental );

        // // compute history of residual forces due to time integration
        base::time::computeResidualForceHistory<FTB,MSM>( laplace, 
                                                          quadrature, solver,
                                                          fieldBinder, step );
        
        if ( withAdvection ) {
            // Compute system matrix from Convection
            base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver, 
                                                         fieldBinder, convection );


            // compute history of residual forces due to time integration
            base::time::computeResidualForceHistory<FTB,MSM>( convection, 
                                                              quadrature, solver,
                                                              fieldBinder, step );
        }

        // Finalise assembly
        solver.finishAssembly();

        // Solve
        if ( withAdvection ) // unsymmetric system matrix
            solver.superLUSolve();
        else                 //   symmetric system matrix
            solver.choleskySolve();

        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, concentration );

        // pass to history
        base::dof::pushHistory( concentration );

        // write output
        writeData( mesh, concentration, product, quadrature,
                   outFile, step+1, stepSize, x1, x2, writeEvery );

    }

    return 0;
}
