#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/CreateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <base/dof/scaleConstraints.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/setField.hpp>
#include <base/fe/Basis.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

#include <base/asmb/BodyForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/Lame.hpp>

#include <solid/HyperElastic.hpp>
#include <fluid/Stokes.hpp>
#include <heat/Laplace.hpp>

#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

#include <base/post/ErrorNorm.hpp>

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD1, typename FIELD2>
void writeVTK( const MESH& mesh,
               const FIELD1& displacement,
               const FIELD2& pressure,
               const std::string& meshFile,
               const unsigned numStep )
{
    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );
    // create file name with step number
    const std::string vtkFile = baseName + "." +
        base::io::leadingZeros( numStep ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );

    base::io::vtk::writePointData( vtkWriter, mesh, displacement, "disp" );
    base::io::vtk::writePointData( vtkWriter, mesh, pressure,     "pressure" );
    vtk.close();
}


//==============================================================================
const double A = 0.;
const double B = 0.;
const double alpha = 1.;
const double beta  = 2.;

template<unsigned DIM>
double f( const typename base::Vector<DIM>::Type& x )
{
    return std::sin( alpha * x[0] + A ) * std::sin( beta * x[1] + B );
}

template<unsigned DIM>
typename base::Vector<DIM>::Type
gradF( const typename base::Vector<DIM>::Type& x )
{
    typename base::Vector<DIM>::Type grad = base::Vector<DIM>::Type::Zero();
    grad[0] = alpha * std::cos( alpha*x[0] + A ) * std::sin( beta*x[1] + B );
    grad[1] = beta  * std::sin( alpha*x[0] + A ) * std::cos( beta*x[1] + B );
    return grad;
}

//------------------------------------------------------------------------------
const double E     = 1.0;
const double nu    = 0.0;
const double kappa = 1.0;

const double C = kappa *
    (mat::Lame::lambda(E, nu) + 2. * mat::Lame::mu(E,nu) )
    * (alpha*alpha + beta*beta);

double g(    const double t ){ return      std::exp( -C * t ); }
double gDot( const double t ){ return -C * std::exp( -C * t ); }

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM>::Type
u( const typename base::Vector<DIM>::Type& x, const double time )
{
    return gradF<DIM>(x) * g(time);
}

template<unsigned DIM>
base::Vector<1>::Type
p( const typename base::Vector<DIM>::Type& x, const double time )
{
    base::Vector<1>::Type result;
    result[0] = (1./kappa) * f<DIM>(x) * gDot(time);
    return result;
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void setU( const typename base::Vector<DIM>::Type& x, DOF* doFPtr,
           const double time ) 
{
    const typename base::Vector<DIM>::Type U = u<DIM>(x,time);
    doFPtr -> setValue( 0, U[0] );
    doFPtr -> setValue( 1, U[1] );

    doFPtr -> pushHistory();
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void dirichletU( const typename base::Vector<DIM>::Type& x, DOF* doFPtr,
                 const double time ) 
{
    const typename base::Vector<DIM>::Type U = u<DIM>(x,time);
    doFPtr -> constrainValue( 0, U[0] );
    doFPtr -> constrainValue( 1, U[1] );
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void setP( const typename base::Vector<DIM>::Type& x, DOF* doFPtr,
           const double time ) 
{
    const double P = ( p<DIM>(x, time))[0];
    doFPtr -> setValue( 0, P );

    doFPtr -> pushHistory();
}

//------------------------------------------------------------------------------
template<unsigned DIM>
base::Vector<1>::Type
neumannBC( const typename base::Vector<DIM>::Type& x,
           const typename base::Vector<DIM>::Type& normal,
           const double time )
{
    const typename base::Vector<DIM>::Type gF = gradF<DIM>(x);
    base::Vector<1>::Type result;
    result[0] = gDot( time ) * gF.dot( normal );
    return result;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM>::Type
bodyForce( const typename base::Vector<DIM>::Type& x, const double time )
{
    return base::constantVector<DIM>( 0. );
}

//------------------------------------------------------------------------------
template<typename QUAD, typename MESH, typename DISP, typename PRESS>
void error( const QUAD& quadrature,   const MESH&  mesh,
            const DISP& displacement, const PRESS& pressure,
            const double time )
{
    // compute L2-error and tell it to the user
    std::cout << time << "  "
              << base::post::errorComputation<0>(
                  quadrature, mesh, displacement,
                  boost::bind( &u<QUAD::dim>, _1, time ) )
              << "  "
              << base::post::errorComputation<0>(
                  quadrature, mesh, pressure,
                  boost::bind( &p<QUAD::dim>, _1, time ) )
              << '\n';
}

//------------------------------------------------------------------------------
template<typename MESH, typename DISP, typename PRESS>
void probe( const MESH& mesh, const DISP& displacement, const PRESS& pressure,
            const double time)
{
    // pick an element
    const std::size_t numElements = std::distance( mesh.elementsBegin(), mesh.elementsEnd() );
    const std::size_t elemNum = (numElements + std::sqrt(numElements)) / 2;

    // geometry
    const typename MESH::Element* geomEp = mesh.elementPtr( elemNum );
    const typename base::Vector<MESH::Element::dim>::Type xi =
        base::constantVector<MESH::Element::dim>( 0. );

    const typename base::Geometry<typename MESH::Element>::result_type x =
        base::Geometry<typename MESH::Element>()( geomEp, xi );

    std::cout << time << " ";

    // displacement
    const typename DISP::Element* dispEp = displacement.elementPtr( elemNum );
    const typename base::Vector<DISP::DegreeOfFreedom::size>::Type uh =
        base::post::evaluateField( geomEp, dispEp, xi );
    const typename base::Vector<DISP::DegreeOfFreedom::size>::Type ux =
        u<MESH::Element::dim>( x, time );
    
    std::cout << uh[0] << " " << uh[1] << "  "
              << ux[0] << " " << ux[1] << "  ";

    // pressure
    const typename PRESS::Element* pressEp = pressure.elementPtr( elemNum );
    const typename base::Vector<PRESS::DegreeOfFreedom::size>::Type ph =
        base::post::evaluateField( geomEp, pressEp, xi );
    const typename base::Vector<PRESS::DegreeOfFreedom::size>::Type px =
        p<MESH::Element::dim>( x, time );
    std::cout << ph[0] << "  " << px[0] << '\n';
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // basic attributes of the computation
    const unsigned    geomDeg   = 1;
    const unsigned    fieldDegU = 2;
    const unsigned    fieldDegP = 1;
    const unsigned    tiOrder   = 2;   // order of time integrator
    const base::Shape shape     = base::QUAD;

    typedef  base::time::BDF<tiOrder> MSM;
    //typedef base::time::AdamsMoulton<tiOrder> MSM;

    // time stepping method determines the history size
    const unsigned nHist = MSM::numSteps;


    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  input.dat \n";
        return 0;
    }

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    //double E, nu, kappa,
    double  deltaT;
    unsigned numSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        //prop.registerPropertiesVar( "E",                E );
        //prop.registerPropertiesVar( "nu",               nu );
        //prop.registerPropertiesVar( "kappa",            kappa );
        prop.registerPropertiesVar( "deltaT",           deltaT );
        prop.registerPropertiesVar( "numSteps",         numSteps );

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

    //--------------------------------------------------------------------------
    // define a mesh
    typedef base::Unstructured<shape,geomDeg>   Mesh;
    const unsigned dim = Mesh::Node::dim;
    
    // create a mesh and read from input
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // quadrature objects for volume and surface
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // Create a displacement field
    const unsigned    doFSizeU = dim;
    typedef base::fe::Basis<shape,fieldDegU>           FEBasisU;
    typedef base::Field<FEBasisU,doFSizeU,nHist>       Displacement;
    typedef Displacement::DegreeOfFreedom              DoFU;
    Displacement displacement;

    // Create a pressure field
    const unsigned    doFSizeP = 1;
    typedef base::fe::Basis<shape,fieldDegP>           FEBasisP;
    typedef base::Field<FEBasisP,doFSizeP,nHist>       Pressure;
    typedef Pressure::DegreeOfFreedom                  DoFP;
    Pressure pressure;

    // generate DoFs from mesh
    base::dof::generate<FEBasisU>( mesh, displacement );
    base::dof::generate<FEBasisP>( mesh, pressure );

    // Creates a list of <Element,faceNo> pairs along the boundary
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Create a boundary mesh from this list
    typedef base::mesh::CreateBoundaryMesh<Mesh::Element> CreateBoundaryMesh;
    typedef CreateBoundaryMesh::BoundaryMesh BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        CreateBoundaryMesh::apply( meshBoundary.begin(),
                                   meshBoundary.end(),
                                   mesh, boundaryMesh );
    }

    // material object
    typedef mat::hypel::StVenant Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    typedef base::asmb::FieldBinder<Mesh,Displacement,Pressure> Field;
    Field field( mesh, displacement, pressure );
    
    typedef Field::TupleBinder<1,1>::Type TopLeft;
    typedef Field::TupleBinder<1,2>::Type TopRight;
    typedef Field::TupleBinder<2,1>::Type BotLeft;
    typedef Field::TupleBinder<2,2>::Type BotRight;
    
    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Pressure> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, pressure );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SFTB;
    
    // kernel objects
    typedef solid::HyperElastic<Material,TopLeft::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    typedef fluid::PressureGradient<TopRight::Tuple> GradP;
    GradP gradP;

    typedef fluid::VelocityDivergence<BotLeft::Tuple> DivU;
    DivU divU( true );

    typedef heat::Laplace<BotRight::Tuple> Laplace;
    Laplace laplace( kappa );

    // constrain the boundary
    base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                            meshBoundary.end(),
                                            mesh, displacement, 
                                            boost::bind( &dirichletU<dim,DoFU>, _1, _2, 0.0 ) );

    // Fix one pressure dof (Neumann problem of the pressure!)
    Pressure::DoFPtrIter pIter = pressure.doFsBegin();
    (*pIter) -> constrainValue( 0, 0.0 );

    // Number the degrees of freedom
    const std::size_t numDoFsU =
        base::dof::numberDoFsConsecutively( displacement.doFsBegin(), displacement.doFsEnd() );
    std::cout << "# Number of displacement dofs " << numDoFsU << std::endl;
    const std::size_t numDoFsP =
        base::dof::numberDoFsConsecutively( pressure.doFsBegin(), pressure.doFsEnd(), numDoFsU );
    std::cout << "# Number of pressure     dofs " << numDoFsP << std::endl;

    // set initial conditions
    base::dof::setField( mesh, displacement,
                         boost::bind( &setU<dim,DoFU>, _1, _2, 0.0 ) );

    base::dof::setField( mesh, pressure,
                         boost::bind( &setP<dim,DoFP>, _1, _2, 0.0 ) );


    // write VTK file
    writeVTK( mesh, displacement, pressure, meshFile, 0 );

    error( quadrature, mesh, displacement, pressure, 0. );

    //probe( mesh, displacement, pressure, 0. );
    
    for ( unsigned n = 0; n < numSteps; n++ ) {

        //std::cout << "Step " << n << " ";;

        const double time = (n+1) * deltaT;

        //std::cout << time << "  ";

        // clear all the constraints
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoFU::clearConstraints, _1 ) );

        // initialise current displacement to zero (linear elasticity)
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoFU::clearValue, _1 ) );

        // constrain the boundary
        base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                                meshBoundary.end(),
                                                mesh, displacement, 
                                                boost::bind( &dirichletU<dim,DoFU>, _1, _2, time ) );


    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDoFsU + numDoFsP );

        //------------------------------------------------------------------
        base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature, solver, 
                                                         field, hyperElastic );

        base::asmb::stiffnessMatrixComputation<TopRight>( quadrature, solver,
                                                          field, gradP );

        base::asmb::stiffnessMatrixComputation<BotRight>( quadrature, solver,
                                                          field, laplace);


        base::time::computeReactionTerms<BotLeft,MSM>( divU, quadrature, solver,
                                                       field, deltaT, n );

        
        //typedef base::time::ReactionTerms<Quadrature,Solver,MSM,BotLeft> RT;
        //RT::Kernel kernel = boost::bind( &DivU::tangentStiffness, &divU,
        //                                 _1, _2, _3, _4 );
        //RT rt( kernel, quadrature, solver, deltaT, n, false );
        //std::for_each( field.elementsBegin(), field.elementsEnd(), rt );

#if 0
        // Compute RHS terms from history of forces
        {

            base::time::ResidualForceHistory<HyperElastic,Quadrature,Solver,MSM>
                rfh1( hyperElastic, quadrature, solver, n );
            std::for_each( fieldUUP.elementsBegin(), fieldUUP.elementsEnd(), rfh1 );

            base::time::ResidualForceHistory<GradP,Quadrature,Solver,MSM>
                rfh2( gradP, quadrature, solver, n );
            std::for_each( fieldUPU.elementsBegin(), fieldUPU.elementsEnd(), rfh2 );

            base::time::ResidualForceHistory<Laplace,Quadrature,Solver,MSM>
                rfh3( laplace, quadrature, solver, n );
            std::for_each( fieldPPU.elementsBegin(), fieldPPU.elementsEnd(), rfh3 );
        }
#endif

        // Body force
        base::asmb::bodyForceComputation<TopLeft>( quadrature, solver, field,
                                                   boost::bind( &bodyForce<dim>, _1, time ) );

        // Neumann boundary condition
        base::asmb::neumannForceComputation<SFTB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                                   boost::bind( &neumannBC<dim>, _1, _2, time ) );


        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.superLUSolve();
            
        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, displacement );
        base::dof::setDoFsFromSolver( solver, pressure );
        
        // push history
        std::for_each( displacement.doFsBegin(), displacement.doFsEnd(),
                       boost::bind( &DoFU::pushHistory, _1 ) );
        std::for_each( pressure.doFsBegin(), pressure.doFsEnd(),
                       boost::bind( &DoFP::pushHistory, _1 ) );


        // write VTK file
        writeVTK( mesh, displacement, pressure, meshFile, n+1 );
        
        // Finished time steps
        //--------------------------------------------------------------------------

        // compute L2-error and tell it to the user
        error( quadrature, mesh, displacement, pressure, time );

        //probe( mesh, displacement, pressure, time );

    }
    
    return 0;
}
