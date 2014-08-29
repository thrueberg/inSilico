#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/io/PropertiesParser.hpp>
#include <base/io/smf/Reader.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/Lame.hpp>

#include <base/BoundaryValueProblem.hpp>
#include <solid/CompressibleDriver.hpp>

#include "../generateMesh.hpp"
#include "extractInternalMeshSurface.hpp"

#include "NeumannForce2.hpp"

#include "AdjointBodyForce.hpp"
#include "InternalTraction.hpp"

#include "ElasticityDriver.hpp"
#include "CostFunction.hpp"

#include "collapseMesh.hpp"

//------------------------------------------------------------------------------
namespace tfm{

    const double coordTol = 1.e-3;

    // Fix (parts of) the boundary
    template<unsigned DIM, typename DOF>
    void fixBoundary( const typename base::Vector<DIM>::Type& x,
                      DOF* doFPtr )
    {
        const bool lr =
            (std::abs( x[0] - 0. ) < coordTol ) or
            (std::abs( x[0] - 1. ) < coordTol );

        if ( not lr ) return;
        
        for ( unsigned d = 0; d < DIM; d++ ) { 
            if ( doFPtr -> isActive(d) )
                doFPtr -> constrainValue( d, 0. );
        }
    }

    // Filter for element faces belonging to the internal surface
    template<unsigned DIM>
    bool internalSurface( const typename base::Vector<DIM>::Type& x )
    {
        //if ( std::abs( x[0] - 0.5 ) < coordTol ) return true;

        if ( ( std::abs( x[0] - 0.25 ) < coordTol ) or
             ( std::abs( x[0] - 0.75 ) < coordTol ) ) {
            if ( ( x[1] >= 0.25 ) and ( x[1] <= 0.75 ) ) return true;
        }
        
        if ( ( std::abs( x[1] - 0.25 ) < coordTol ) or
             ( std::abs( x[1] - 0.75 ) < coordTol ) ) {
            if ( ( x[0] >= 0.25 ) and ( x[0] <= 0.75 ) ) return true;
        }

        return false;
    }

    // Measure filter
    template<unsigned DIM>
    bool measureHere( const typename base::Vector<DIM>::Type& x )
    {
        if ( x[1] < x[0] - 0.5 ) return true;
        if ( x[1] > x[0] + 0.5 ) return true;
        return false;
    }


        
    int tfm( int argc, char * argv[] );

}


//------------------------------------------------------------------------------
int tfm::tfm( int argc, char * argv[] )
{
    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  input.dat \n";
        return 0;
    }

    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 1;
    const base::Shape shape    = base::HyperCubeShape<SPACEDIM>::value;

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    double E, nu, traction, delta, stepSize;
    unsigned numElements, maxIter;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "traction",         traction );
        prop.registerPropertiesVar( "numElements",      numElements );
        prop.registerPropertiesVar( "delta",            delta );
        prop.registerPropertiesVar( "maxIter",          maxIter );
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

    //--------------------------------------------------------------------------
    // define a mesh
    typedef base::Unstructured<shape,geomDeg>    Mesh;
    const unsigned dim = Mesh::Node::dim;
    typedef Mesh::Node::VecDim VecDim;

    // create a mesh and read from input
    Mesh mesh;
    VecDim a, b;
    {
        a = base::constantVector<dim>(0.);
        b = base::constantVector<dim>(1.);
        
        base::Vector<dim,unsigned>::Type N;
        for ( unsigned d = 0; d < dim; d++ ) N[d] = numElements;
        generateMesh<dim>( mesh, N, a, b );
    }

    // elastic material
    typedef mat::hypel::StVenant Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );
    
    // boundary value problems
    typedef base::BoundaryValueProblem<Mesh,dim,fieldDeg> BVP;
    BVP bvpDirect(    mesh );
    BVP bvpAdjoint(   mesh );
    BVP bvpReference( mesh );
 

    // driver for more convenient bvp solve
    typedef tfm::ElasticityDriver<BVP> Driver;
    Driver direct(    bvpDirect,    material );
    Driver adjoint(   bvpAdjoint,   material );
    Driver reference( bvpReference, material );

    // fix boundaries
    direct.dirichlet(    boost::bind( &tfm::fixBoundary<dim,BVP::DoF>, _1, _2 ) );
    adjoint.dirichlet(   boost::bind( &tfm::fixBoundary<dim,BVP::DoF>, _1, _2 ) );
    reference.dirichlet( boost::bind( &tfm::fixBoundary<dim,BVP::DoF>, _1, _2 ) );

    // extract internal faces
    std::vector< std::pair<std::size_t, unsigned> > faceList;
    tfm::extractInternalMeshSurface( mesh.elementsBegin(),
                                     mesh.elementsEnd(),
                                     boost::bind( &tfm::internalSurface<dim>, _1 ),
                                     faceList );

    if ( faceList.empty() ) {
        std::cout << "Cannot create surface mesh with this filter\n";
        return 0;
    }

    // make a mesh from this
    typedef base::mesh::BoundaryMeshBinder<Mesh>::Type SurfaceMesh;
    SurfaceMesh surfaceMesh;
    {
        SurfaceMesh surfaceMeshTmp;
        base::mesh::generateBoundaryMesh( faceList.begin(), faceList.end(),
                                          mesh, surfaceMeshTmp );
    
        tfm::collapseMesh( surfaceMeshTmp, surfaceMesh, 1.e-6 );

        // surface mesh --> copy element pointer
        SurfaceMesh::ElementPtrConstIter eIter = surfaceMeshTmp.elementsBegin();
        SurfaceMesh::ElementPtrConstIter eEnd   = surfaceMeshTmp.elementsEnd();
        SurfaceMesh::ElementPtrIter      eNew   = surfaceMesh.elementsBegin();
        for ( ; eIter != eEnd; ++eIter, ++eNew ) {
            const std::size_t domainID = (*eIter) -> getDomainID();
            (*eNew) -> setDomainElementPointer( mesh.elementPtr( domainID ) );

            // copy parametric coordinates one-by-one
            std::copy( (*eIter) -> parametricBegin(),
                       (*eIter) -> parametricEnd(),
                       (*eNew)  -> parametricBegin() );

        }
    }

    // base::mesh::generateBoundaryMesh( faceList.begin(), faceList.end(),
    //                                   mesh, surfaceMesh );

    // handler for internal tractions
    typedef tfm::InternalTraction<BVP::Field,SurfaceMesh> TH;
    TH tractionHandler( bvpAdjoint.getField(), surfaceMesh, delta );
    tractionHandler.writeVTKFile( "bla", 0 );
    
    // set to unit vector
    base::Vector<dim>::Type t = base::constantVector<dim>( 0. );
    //t[0] = 1.;
    //tractionHandler.setToConstant( t );
    //tractionHandler.randomize();
    tractionHandler.setToBla();
    
    // create reference solution
    {
        reference.neumann( boost::bind( &TH::apply, &tractionHandler, _1, _2 ) );
        reference.solve( surfaceMesh );
        bvpReference.writeVTKFile( "reference" );

        tractionHandler.writeVTKFile( "ref", 0 );
    }

    // reset tractions to zero
    //t[0] = 0.4; t[1] = 0.4;
    //tractionHandler.setToConstant( t );
    tractionHandler.randomize();

    // body force
    typedef tfm::AdjointBodyForce<Mesh,BVP::Field> ABF;
    ABF abf( mesh, bvpDirect.getField(), bvpReference.getField(),
             boost::bind( &measureHere<dim>, _1 ) );

    // cost function
    typedef tfm::CostFunction<Mesh,BVP::Field,TH> CostFunction;
    CostFunction costFunction( mesh, bvpDirect.getField(), bvpReference.getField(),
                               tractionHandler, delta,
                               boost::bind( &measureHere<dim>, _1 ) );

    //
    const double J = costFunction.evaluate();
    boost::array<double,3> costs = {{ J, J, J }};
    double prevStepSize = stepSize;



    //--------------------------------------------------------------------------
    // Optimisation loop
    for ( unsigned iter = 0; iter < maxIter; iter++ ) {

        std::cout << iter << " " << stepSize << "  ";

        // solve direct problem
        direct.neumann( boost::bind( &TH::apply, &tractionHandler, _1, _2 ) );
        direct.solve(   surfaceMesh );
        bvpDirect.writeVTKFile( "direct." + base::io::leadingZeros( iter ) );

        // solve adjoint problem
        adjoint.force( boost::bind( &ABF::apply, &abf, _1, _2 ) );
        adjoint.solve( surfaceMesh );
        bvpAdjoint.writeVTKFile( "adjoint." + base::io::leadingZeros( iter ) );

        const double Jnew = costFunction.evaluate();

        std::cout << Jnew
                  << std::endl;

        // update costs
        {
            costs[2] = costs[1];
            costs[1] = costs[0];
            costs[0] = Jnew;
        }

        // quadratic equation
#if 0
        if ( iter > 5 ) {
            const double alpha =
                ( (costs[0]-costs[1])*stepSize/prevStepSize +
                  (costs[1]-costs[2])*prevStepSize/stepSize) / (stepSize+prevStepSize);

            const double beta =
                ( (costs[0]-costs[1])/stepSize -
                  (costs[1]-costs[2])/prevStepSize) /  (stepSize+prevStepSize);

            
            if ( beta > 0 ) { // positive Hessian!

                const double tmp = stepSize;

                const double candidate = - alpha / 2. / beta;

                if ( candidate < 100. )               
                    stepSize = candidate;

                prevStepSize = tmp;
                
            }
            else
                prevStepSize = stepSize;
        }
#endif        

        // update the traction
        tractionHandler.update( stepSize );
        tractionHandler.writeVTKFile( "surface", iter );
    }
    
    
    return 0;
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    return tfm::tfm( argc, argv );
}
