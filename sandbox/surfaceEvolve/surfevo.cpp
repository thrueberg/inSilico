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
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>

#include <base/solver/Eigen3.hpp>
#include <base/post/evaluateAtNodes.hpp>
#include <base/post/Monitor.hpp>

#include <base/kernel/Mass.hpp>
#include <base/kernel/Laplace.hpp>
#include <surf/Moments.hpp>

#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

//------------------------------------------------------------------------------
template<typename ELEMENT>
typename base::Vector<ELEMENT::Node::dim>::Type
prescribedCurvature( const ELEMENT* gep, 
                     const typename ELEMENT::GeomFun::VecDim& xi,
                     const double value )
{
    typename base::Vector<ELEMENT::Node::dim>::Type normal;
    base::SurfaceNormal<ELEMENT>()( gep, xi, normal );
    return 2. * value * normal;
}

    
//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class SurfaceMoments
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    //@}

    typedef void result_type;
    
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim   LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim  GlobalVecDim;

    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //--------------------------------------------------------------------------
    static void enclosedVolume( const FieldTuple&  fieldTuple,
                                const LocalVecDim& xi,
                                const double       weight,
                                double& volume )
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        

        // Get surface normal and measure
        GlobalVecDim normal;
        const double detJ = 
            base::SurfaceNormal<GeomElement>()( geomEp, xi, normal );

        // get global coordinate
        const typename GeomElement::Node::VecDim x =
            base::Geometry<GeomElement>()( geomEp, xi );

        //
        volume += detJ * weight * normal[0] *x[0];
    }

    //--------------------------------------------------------------------------
    static void centroid( const FieldTuple&  fieldTuple,
                          const LocalVecDim& xi,
                          const double       weight,
                          GlobalVecDim&      centroid )
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();

        // Get surface normal and measure
        GlobalVecDim normal;
        const double detJ = 
            base::SurfaceNormal<GeomElement>()( geomEp, xi, normal );

        // get global coordinate
        const typename GeomElement::Node::VecDim x =
            base::Geometry<GeomElement>()( geomEp, xi );

        //
        for ( unsigned d = 0; d < globalDim; d++ )
            centroid[d] += 0.5 * x[d] * x[d] * normal[d] * detJ * weight;
    }

};

//------------------------------------------------------------------------------
template<typename FTB, typename QUADRATURE, typename FIELDBINDER>
double enclosedMeshVolume( const QUADRATURE&  quadrature,
                           const FIELDBINDER& fieldBinder )
{
    // compute volume
    double result = 0.;
    
    base::asmb::simplyIntegrate<FTB>(
        quadrature, result, fieldBinder,
        boost::bind( &SurfaceMoments<typename FTB::Tuple>::enclosedVolume,
                     _1, _2, _3, _4 ) );

    return result;
}

//------------------------------------------------------------------------------
template<typename FTB, typename QUADRATURE, typename FIELDBINDER>
typename base::Vector<FIELDBINDER::Mesh::Node::dim>::Type
enclosedMeshCentroid( const QUADRATURE&  quadrature,
                      const FIELDBINDER& fieldBinder )
{
    static const unsigned dim = FIELDBINDER::Mesh::Node::dim;
    typename base::Vector<dim>::Type centroid = base::constantVector<dim>(0.);
            
    base::asmb::simplyIntegrate<FTB>(
        quadrature, centroid, fieldBinder,
        boost::bind( &SurfaceMoments<typename FTB::Tuple>::centroid,
                     _1, _2, _3, _4 ) );

    return centroid;
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void fixSomePart( const typename base::Vector<DIM,double>::Type& x,
                 DOF* doFPtr ) 
{
    bool lower = true;
    bool upper = true;
    for ( unsigned d = 0; d < DIM; d++ ) {
        if ( x[d] > 0. ) lower = false;
        if ( x[d] < 0. ) upper = false;
    }
    
    const bool fix = false; //lower or upper;    
    
    if ( ( doFPtr -> isActive(0) ) and fix ) {

        for ( unsigned d = 0; d < DIM; d++ )
            doFPtr -> constrainValue( d, x[d] );

    }
        
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM,double>::Type& x,
                  DOF* doFPtr ) 
{
    if ( doFPtr -> isActive(0) ) {

        for ( unsigned d = 0; d < DIM; d++ )
            doFPtr -> constrainValue( d, x[d] );

    }
        
}

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void currentLoc( const typename base::Vector<DIM,double>::Type& x,
                 DOF* doFPtr ) 
{
    for ( unsigned d = 0; d < DIM; d++ )
        doFPtr -> setValue( d, x[d] );

    doFPtr -> pushHistory();
}

//------------------------------------------------------------------------------
template<typename MESH>
void writeVTKFile( const std::string baseName,
                   const unsigned step,
                   const MESH&  mesh )
{
    // VTK Legacy
    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );
    vtk.close();
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    if ( argc != 6 ) {
        std::cout << "Usage:  " << argv[0] << " mesh.smf nSteps dt visc kappa\n\n";
        return -1;
    }

    const std::string meshFile = boost::lexical_cast<std::string>( argv[1] );
    const unsigned    nSteps   = boost::lexical_cast<unsigned   >( argv[2] );
    const double      dt       = boost::lexical_cast<double     >( argv[3] );
    const double      visc     = boost::lexical_cast<double     >( argv[4] );
    const double      kappa    = boost::lexical_cast<double     >( argv[5] );

    const std::string baseName = base::io::baseName( meshFile, ".smf" );
    
    //--------------------------------------------------------------------------
    const unsigned    dim      = 2;
    const unsigned    geomDeg  = 1;   // degree of geometry approximation
    const unsigned    tiOrder  = 1;   // order of time integrator
    const base::Shape shape    =
        base::SimplexShape<dim-1>::value;

    // choose a time stepping method
    typedef  base::time::BDF<tiOrder> MSM;
    //typedef base::time::AdamsMoulton<tiOrder> MSM;

    // time stepping method determines the history size
    const unsigned nHist = MSM::numSteps;
    
    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg,dim>     Mesh;

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
    const unsigned    doFSize  = dim;

    // DOF handling 
    typedef base::Field<FEBasis,doFSize,nHist>        Location;
    Location location;
    base::dof::generate<FEBasis>( mesh, location );

    // set initial condition to the identity 
    base::dof::setField( mesh, location,
                         boost::bind( &currentLoc<dim,
                                      Location::DegreeOfFreedom>, _1, _2 ) );

    // set initial condition to the identity 
    base::dof::setField( mesh, location,
                         boost::bind( &fixSomePart<dim,
                                      Location::DegreeOfFreedom>, _1, _2 ) );
    
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // empty mesh boundary implies a closed surface
    const bool closedSurface =
        (std::distance( meshBoundary.begin(), meshBoundary.end() )
         == 0);

    // constrain the boundary degrees of freedom
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, location, 
                                           boost::bind( &dirichletBC<dim,
                                                        Location::DegreeOfFreedom>,
                                                        _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( location.doFsBegin(),
                                            location.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs
              << ", number of time steps " << nSteps
              << std::endl;

    typedef base::asmb::FieldBinder<Mesh,Location> FieldBinder;
    FieldBinder fieldBinder( mesh, location );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // Object for the surface laplacian
    typedef base::kernel::Laplace<FTB::Tuple> LaplaceBeltrami;
    LaplaceBeltrami laplace( visc );

    typedef base::kernel::Mass<FTB::Tuple> Mass;
    Mass          mass( 1.0 );

    // write initial state
    writeVTKFile( baseName, 0, mesh );

    const double initialVolume =
        ( closedSurface ? enclosedMeshVolume<FTB>( quadrature, fieldBinder ) :
          1.0 );

    //--------------------------------------------------------------------------
    // Time loop
    for ( unsigned step = 0; step < nSteps; step ++ ) {

        //std::cout << "Time step " << step << std::endl;
        const double time = step * dt;
        std::cout << time << " ";
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        // Compute system matrix from Laplacian
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                     fieldBinder,
                                                     laplace );

        
        // compute inertia terms, d/dt, due to time integration
        base::time::computeReactionTerms<FTB,MSM>( mass, quadrature, solver,
                                                   fieldBinder, dt, step );

        base::asmb::bodyForceComputation2<FTB>( quadrature, solver, fieldBinder,
                                                boost::bind( &prescribedCurvature<Mesh::Element>,
                                                             _1, _2, kappa ) );



        // compute history of residual forces due to time integration
        //base::time::computeResidualForceHistory<FieldTupleBinder,MSM>( staticHeat,
        //                                                               quadrature, solver,
        //                                                               fieldBinder, step );

        // Finalise assembly
        solver.finishAssembly();

        // Solve a possibly non-symmetric system
        solver.superLUSolve();

        // distribute results back to dofs
        {
            base::dof::setDoFsFromSolver( solver, location );

        }

        // Pass 'location' datum to coordiantes of the mesh
        {
            std::vector<Mesh::Node::VecDim> newCoordinates;
            base::post::evaluateFieldAtNodes( mesh, location, newCoordinates );

            std::vector<Mesh::Node::VecDim>::iterator newIter = newCoordinates.begin();
            Mesh::NodePtrIter begin = mesh.nodesBegin();
            Mesh::NodePtrIter   end = mesh.nodesEnd();
            for ( ; begin != end; ++begin, ++newIter ) {

                (*begin) -> setX( &( (*newIter)[0] ) );
            }
        }
        
        // Scaling
        if ( false ) {

            // compute volume
            const double volume = enclosedMeshVolume<FTB>( quadrature, fieldBinder );
            
            // compute geometric centroid
            base::Vector<dim>::Type centroid =
                enclosedMeshCentroid<FTB>( quadrature, fieldBinder );

            centroid /= volume;

            // rescale coordinates to maintain the volume
            const double scaling = std::pow( (initialVolume/volume),
                                             1./static_cast<double>( dim ) );

            //
            Mesh::NodePtrIter begin = mesh.nodesBegin();
            Mesh::NodePtrIter   end = mesh.nodesEnd();
            for ( ; begin != end; ++begin ) {
                base::Vector<dim>::Type oldX;
                (*begin) -> getX( &(oldX[0]) );

                const base::Vector<dim>::Type newX =
                    centroid + scaling * (oldX - centroid);
                
                (*begin) -> setX( &(newX[0]) );
            }
            
            // compute volume
            const double volume2 = enclosedMeshVolume<FTB>( quadrature, fieldBinder );
            
            // bla
            std::cout << "V = " << volume << " --> " << volume2
                      << ", X = " << centroid.transpose();
            
            // pass back the rescaled coordinates to the degrees of freedom
            base::dof::setField( mesh, location,
                                 boost::bind( &currentLoc<dim,
                                              Location::DegreeOfFreedom>, _1, _2 ) );

        }

        std::cout << std::endl;

        // push history
        std::for_each( location.doFsBegin(), location.doFsEnd(),
                       boost::bind( &Location::DegreeOfFreedom::pushHistory, _1 ) );

        
        writeVTKFile( baseName, step+1, mesh );



        
    } // end time loop

    return 0;
}
