//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smfRandom.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <ctime>
// boost includes
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
// base includes
#include <base/verify.hpp>
// tools includes
#include <tools/converter/smf2xx/SMFHeader.hpp>
#include <tools/converter/smf2xx/Conversion.hpp>

#include <base/linearAlgebra.hpp>
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/Size.hpp>
#include <base/fe/Policies.hpp>
#include <base/fe/Basis.hpp>
// base/io includes
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace smfRandom{

            //------------------------------------------------------------------
            // fix the boundary ?
            static bool fixBoundary = true;
            // how far to push the nodes
            static double maxDist   = 0.1;
            
            // random number generator, seeded with time 
            static boost::random::mt19937 rng( static_cast<unsigned>( std::time(0)) );

            // return a random number from [-1,1]
            static double randomNum()
            {
                // use boost to get a random number
                boost::random::uniform_int_distribution<> randInt;
                const int x = randInt( rng );

                // random double from [0,1]
                const double y =
                    static_cast<double>( x ) / static_cast<double>( randInt.max() );

                // map from [0,1] to [-1,1]
                return 2. * y - 1.;
            }
            
            // implemenation below
            template<base::Shape SHAPE,unsigned DEGREE>
            struct Converter;
        }
    }
}


    
//------------------------------------------------------------------------------
/** Read smf file, apply random coordinate shift and write new smf file.
 *  A mesh in SMF format is read as usual and, before outputting to the same
 *  format, the coordinates in the mesh are displacement in a random fashion.
 *  The nodes are pushed in a random direction by a random distance. Here, the
 *  spatial dimension of the problem is deduced from the dimension of the shape
 *  of the elements.
 *  \note This approach will fail in case of a surface mesh
 *
 *  In detail, the method does the following
 *  - find the minimal mesh size as the upper bound of random displacement
 *  - determine the maximal random displacement \f$ d \f$ as the minimal mesh
 *    size times a number from [0,1] provided by the user
 *  - find the nodes which are on the boundary, they are not moved
 *  - go through all nodes, get their given coordinate \f$ x \f$, construct a
 *    random direction vector, normalise it \f$ e \f$ and do
 *    \f[
 *          x \to  x + d \times e
 *    \f]
 * 
 *  \tparam SHAPE  Type of element shape
 *  \tparam DEGREE Polynomial degree of the element
 */
template<base::Shape SHAPE,unsigned DEGREE>
struct tools::converter::smfRandom::Converter
{
    static void apply( std::istream& smfIn,
                       std::ostream& smfOut )
    {
        // Attributes of the mesh
        static const base::Shape shape   = SHAPE;
        static const unsigned degree     = DEGREE;
        static const unsigned    dim     = 3;

        // the manifold dimension
        static const unsigned elemDim    = base::ShapeDim<shape>::value;
        
        // Mesh type and object
        typedef base::Unstructured<shape,degree,dim>  Mesh;
        Mesh mesh;

        // SMF input
        base::io::smf::readMesh( smfIn, mesh );

        //----------------------------------------------------------------------
        // determine mesh size
        double hMin = std::numeric_limits<double>::max();
        {
            typename Mesh::ElementPtrConstIter eIter = mesh.elementsBegin();
            typename Mesh::ElementPtrConstIter eEnd  = mesh.elementsEnd();
            for ( ; eIter != eEnd; ++eIter ) {

                // size of this element
                const double h =
                    base::mesh::Size<typename Mesh::Element>::apply( *eIter );

                // store minimal value
                if ( h < hMin ) hMin = h;
            }
        }

        // maximal random displacement
        const double randDisp = hMin * maxDist;

        //----------------------------------------------------------------------
        // boundary nodes
        std::set<std::size_t> boundaryNodes;
        if ( fixBoundary ) {
            
            // Creates a list of <Element,faceNo> pairs
            base::mesh::MeshBoundary meshBoundary;
            meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

            // Type of face for extraction
            static const base::NFace surface = base::ShapeSurface<shape>::value;

            // Helper to get face extraction done
            typedef base::fe::Basis<shape,degree> FEBasis;
            typedef typename FEBasis::FiniteElement FiniteElement;
            typedef base::fe::FaceExtraction<FiniteElement,
                                             surface> FaceExtraction;

            // go through all boundary faces
            typename base::mesh::MeshBoundary::BoundConstIter bIter =
                meshBoundary.begin();
            for ( ; bIter != meshBoundary.end(); ++bIter ) {

                // element and local face number
                const std::size_t elemNo = bIter -> first;
                const unsigned    faceNo = bIter -> second;

                // extract the local node numbers from the face
                std::vector<unsigned> faceIndices;
                FaceExtraction::apply( faceNo, faceIndices );

                // access to element and its nodes
                typename Mesh::Element* geomEp = mesh.elementPtr( elemNo );
                typename Mesh::Element::NodePtrConstIter nFirst =
                    geomEp -> nodesBegin();

                // go through all nodes of the face
                for ( std::size_t f = 0; f < faceIndices.size(); f++ ) {
                    // local node number
                    const unsigned localNum = faceIndices[f];
                    // iterator to corresponding node object
                    typename Mesh::Element::NodePtrConstIter nIter = nFirst;
                    std::advance( nIter, localNum );
                    // store global node ID
                    boundaryNodes.insert( (*nIter) -> getID() );
                }
            }
        } // if boundary shall be fixed

        //----------------------------------------------------------------------
        // go through nodes and modify coordinates
        typename Mesh::NodePtrIter nIter = mesh.nodesBegin();
        typename Mesh::NodePtrIter nEnd  = mesh.nodesEnd();
        for (; nIter != nEnd; ++nIter ) {
            // global node ID
            const std::size_t nodeID = (*nIter) -> getID();
            // check if to be fixed
            const typename std::set<std::size_t>::iterator find =
                boundaryNodes.find( nodeID );

            if ( find == boundaryNodes.end() ) {
                // get given node coordinate
                typename Mesh::Node::VecDim xOld;
                (*nIter) -> getX( &(xOld[0]) );

                // generate a direction vector
                typename Mesh::Node::VecDim dir =
                    base::constantVector<dim>( 0. );
                for ( unsigned d = 0; d < elemDim; d++ ) {
                        dir[d] = randomNum();
                }

                // normalise
                const double dirNorm = dir.norm();
                if ( dirNorm > 1.e-10 ) 
                    dir /= dir.norm();

                // generate a distance
                const double dist = randomNum() * randDisp;
                
                // displace the given coordinate
                const typename Mesh::Node::VecDim xNew =
                    xOld + dist * dir;
                (*nIter) -> setX( &(xNew[0]) );
            }
        }

        // SMF output
        base::io::smf::writeMesh( mesh, smfOut );
    }

};


//------------------------------------------------------------------------------
/** Read smf formatted file, create a temporary mesh and write transformed mesh
 */
int main( int argc, char * argv[] )
{
    namespace smfRandom = tools::converter::smfRandom;
    namespace smf2xx    = tools::converter::smf2xx;
    
    // Sanity check of the number of input arguments
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0]
                  << " file.smf \n\n";
        return 0;
    }

    std::cout << "--------------------------------------------\n"
              << "**  Random motion of the mesh coordinates **\n"
              << "--------------------------------------------\n";

    // Name of smf input file, its basename and the data output file name
    const std::string smfFileIn  = boost::lexical_cast<std::string>( argv[1] );
    const std::string base       = smfFileIn.substr(0, smfFileIn.find( ".smf") );
    const std::string smfFileOut  = base + ".rand.smf";

    // Element attributes
    base::Shape elementShape;
    unsigned    elementNumPoints;
    
    {
        // extract data from header
        std::ifstream smf( smfFileIn.c_str() );
        smf2xx::readSMFHeader( smf, elementShape, elementNumPoints );
        smf.close();
    }

    // User input of random act
    {
        // Not fixing the boundaries does not make sense, I think
        // std::cout << "Fix the boundaries (0/1) ? ";
        // std::cin  >> smfRandom::fixBoundary;

        std::cout << "Maximal displacement (relative to minimal mesh size) ? ";
        std::cin  >> smfRandom::maxDist;

        VERIFY_MSG( (smfRandom::maxDist >= 0.) and (smfRandom::maxDist <  1.),
                    "Distance out of interval [0,1) !!" );
    }

    // Input and output file streams
    std::ifstream smfIn(   smfFileIn.c_str() );
    std::ofstream smfOut( smfFileOut.c_str() );

    // write to file for traceback
    smfOut << "# Generated by smfRandom \n";

    // Call generic conversion helper
    smf2xx::Conversion< smfRandom::Converter >::apply( elementShape,
                                                       elementNumPoints,
                                                       smfIn, smfOut );

    // Close the streams
    smfIn.close();
    smfOut.close();
    
    return 0;
}
