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
#include <base/auxi/compareNumbers.hpp>

#include "../../generateMesh.hpp"
#include "../extractInternalMeshSurface.hpp"
#include "../NeumannForce2.hpp"
#include "../AdjointBodyForce.hpp"
#include "../InternalTraction.hpp"
#include "../ElasticityDriver.hpp"
#include "../CostFunction.hpp"
#include "../collapseMesh.hpp"

//------------------------------------------------------------------------------
namespace tfm{

    const double coordTol = 1.e-3;

    // Filter for element faces belonging to the internal surface
    template<typename MESH>
    class InternalSurface
    {
    public:

        static const unsigned numVerticesPerFace =
            base::NumNFaces<base::FaceShape<MESH::Element::shape>::value,
                            base::VERTEX>::value;

        InternalSurface( std::ifstream& surfIn )
        {
            while( not surfIn.eof() )
            {
                std::vector< std::size_t > face( numVerticesPerFace );
                    
                for ( unsigned v = 0; v < numVerticesPerFace; v++ ) {
                    surfIn >> face[v]; 
                }

                faces_.push_back( face );
            }
        }
        
        bool apply( const std::vector<std::size_t>& face )
        {
            // for ( unsigned v = 0; v < numVerticesPerFace; v++ ) {
            //     std::cout << face[v] << "  ";
            // }
            // std::cout << std::endl;
            
            typename std::vector<std::vector< std::size_t > >::iterator
                it = std::find( faces_.begin(), faces_.end(), face );

            if ( it != faces_.end() ) return true;
            
            return false;
        }

    private:
        std::vector< std::vector<std::size_t> > faces_;
    };

    // Measure filter
    template<unsigned DIM>
    bool measureHere( const typename base::Vector<DIM>::Type& x )
    {
        if ( (x[0] > 45.) and (x[0] < 75.) ) {
            if ( (x[1] > 15.) and (x[1] < 50.) ) {
                return false;
            }
        }

        return true;
    }

    template<typename MESH,typename FIELD>
    void setMeasuredDisplacements( const MESH& mesh,
                                   FIELD& field,
                                   std::ifstream& din )
    {
        while ( not din.eof() ) {
            
            std::size_t id;
            double ux, uy;

            din >> id >> ux >> uy;

            if ( tfm::measureHere<MESH::Node::dim>( mesh.nodePtr( id ) -> getX() ) ) {

                field.doFPtr( id ) -> setValue( 0, ux );
                field.doFPtr( id ) -> setValue( 1, uy );

                VERIFY_MSG( FIELD::DegreeOfFreedom::size == 2,
                            "Only 2D allowed" );
            }

        }
    }



    // Fix (parts of) the boundary
    template<unsigned DIM, typename FIELD>
    void fixBoundary( const typename base::Vector<DIM>::Type& x,
                      typename FIELD::DegreeOfFreedom* doFPtr,
                      const FIELD& field )
    {

        for ( unsigned d = 0; d < FIELD::DegreeOfFreedom::size; d++  ) {
            if ( doFPtr -> isActive( d ) ) {

                const double u = field.doFPtr( doFPtr -> getID() ) -> getValue( d );
                doFPtr -> constrainValue( d, u );
                
            }
        }

        return;
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
    const base::Shape shape    = base::SimplexShape<SPACEDIM>::value;

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    double E, nu, traction, delta, stepSize;
    unsigned maxIter;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "E",                E );
        prop.registerPropertiesVar( "nu",               nu );
        prop.registerPropertiesVar( "traction",         traction );
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
    {
        std::ifstream smfIn( "mesh.smf" );
        base::io::smf::readMesh( smfIn, mesh );
    }

    // elastic material
    typedef mat::hypel::StVenant Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );
    
    // boundary value problems
    typedef base::BoundaryValueProblem<Mesh,dim,fieldDeg> BVP;
    BVP bvpDirect(    mesh );
    BVP bvpAdjoint(   mesh );
    BVP bvpMeasured(  mesh );


    {
        std::ifstream din( "displacements_prueba2D" );
        tfm::setMeasuredDisplacements( mesh, bvpMeasured.getField(), din );   

        bvpMeasured.writeVTKFile( "reference" );
    }
 

    // driver for more convenient bvp solve
    typedef tfm::ElasticityDriver<BVP> Driver;
    Driver direct(    bvpDirect,    material );
    Driver adjoint(   bvpAdjoint,   material );

    // fix boundaries
    direct.dirichlet(    boost::bind( &tfm::fixBoundary<dim,BVP::Field>, _1, _2,
                                      boost::ref( bvpMeasured.getField() ) ) );
    adjoint.dirichlet(   boost::bind( &tfm::fixBoundary<dim,BVP::Field>, _1, _2,
                                      boost::ref( bvpMeasured.getField() ) ) );

    // extract internal faces
    std::ifstream surfIn( "surface" );
    tfm::InternalSurface<Mesh> isurf( surfIn );
    
    std::vector< std::pair<std::size_t, unsigned> > faceList;
    tfm::extractInternalMeshSurface2( mesh.elementsBegin(),
                                      mesh.elementsEnd(),
                                      boost::bind( &tfm::InternalSurface<Mesh>::apply,
                                                   &isurf, _1 ),
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
    base::Vector<dim>::Type t = base::constantVector<dim>( 0. );
    //t[0] = 1.;
    tractionHandler.setToConstant( t );

    {
        BVP    bvpBlub( mesh );
        Driver blub(   bvpBlub,   material );

        // fix boundaries
        blub.dirichlet(    boost::bind( &tfm::fixBoundary<dim,BVP::Field>, _1, _2,
                                        boost::ref( bvpMeasured.getField() ) ) );

        //blub.neumann( boost::bind( &TH::apply, &tractionHandler, _1, _2 ) );
        blub.solve( surfaceMesh );
        bvpBlub.writeVTKFile( "blub" );
    }

    
    // set to unit vector
    //tractionHandler.randomize();
    //tractionHandler.setToBla();


    // reset tractions to zero
    //t[0] = 0.4; t[1] = 0.4;
    //tractionHandler.setToConstant( t );
    tractionHandler.randomize();

    // body force
    typedef tfm::AdjointBodyForce<Mesh,BVP::Field> ABF;
    ABF abf( mesh, bvpDirect.getField(),
             bvpMeasured.getField(), 
             boost::bind( &measureHere<dim>, _1 ) );

    // cost function
    typedef tfm::CostFunction<Mesh,BVP::Field,TH> CostFunction;
    CostFunction costFunction( mesh, bvpDirect.getField(),
                               bvpMeasured.getField(), 
                               tractionHandler, delta,
                               boost::bind( &measureHere<dim>, _1 ) );

    //
    const double J = costFunction.evaluate();

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

        std::cout << Jnew << std::endl;



        // update the traction
        tractionHandler.update( stepSize );
        tractionHandler.writeVTKFile( "surface", iter );

        //if ( iter == 10 ) stepSize *= 2;
    }
    
    
    return 0;
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    return tfm::tfm( argc, argv );
}
