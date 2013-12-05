//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smfBBox.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <fstream>
#include <string>
// boost includes
#include <boost/lexical_cast.hpp>
// base includes
#include <base/verify.hpp>
#include <base/numbers.hpp>
#include <base/Unstructured.hpp>
#include <base/io/smf/Reader.hpp>
// tools includes
#include <tools/converter/smf2xx/SMFHeader.hpp>
#include <tools/converter/smf2xx/Conversion.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace analysis{
        namespace bbox{
        
            //------------------------------------------------------------------
            template<base::Shape SHAPE,unsigned DEGREE>
            struct Analysis
            {
                static void apply( std::istream& smfIn,
                                   std::ostream& out )
                {
                    // Attributes of the mesh
                    static const base::Shape shape   = SHAPE;
                    static const unsigned degree     = DEGREE;
                    static const unsigned    dim     = 3;

                    // Mesh type and object
                    typedef base::Unstructured<shape,degree,dim>  Mesh;
                    Mesh mesh;

                    // SMF input
                    base::io::smf::readMesh( smfIn, mesh );

                    // Corner coordinates
                    typedef typename Mesh::Node::VecDim VecDim;
                    VecDim min, max;
                    for ( unsigned d = 0; d < dim; d++ ) {
                        min[d] = +base::invalidReal();
                        max[d] = -base::invalidReal();
                    }
                    
                    typename Mesh::NodePtrConstIter nIter = mesh.nodesBegin();
                    typename Mesh::NodePtrConstIter nEnd  = mesh.nodesEnd();
                    for ( ; nIter != nEnd; ++nIter ) {

                        VecDim x;
                        (*nIter) -> getX( &(x[0]) );

                        for ( unsigned d = 0; d < dim; d++ ) {
                            if ( x[d] < min[d] ) min[d] = x[d];
                            if ( x[d] > max[d] ) max[d] = x[d];
                        }

                    }

                    // write message
                    const std::string coord( "xyz" );
                    std::cout << "Bounding box: \n";
                    for ( unsigned d = 0; d < dim; d++ ) {

                        std::cout << coord[d] << "-min = " << min[d]
                                  << ", "
                                  << coord[d] << "-max = " << max[d] << "\n";

                    }

                }
            };
            
        }
    }
}
//------------------------------------------------------------------------------
/** Read smf formatted file, compute the bounding box and write to stdout
 */
int main( int argc, char * argv[] )
{
    namespace smf2xx = tools::converter::smf2xx;
    namespace bbox   = tools::analysis::bbox;
    
    // Sanity check of the number of input arguments
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0]
                  << " fileIn.smf \n\n";
        return 0;
    }

    // Name of smf input file, its basename and the data output file name
    const std::string smfFileIn  = boost::lexical_cast<std::string>( argv[1] );

    // Element attributes
    base::Shape elementShape;
    unsigned    elementNumPoints;
    
    {
        // extract data from header
        std::ifstream smf( smfFileIn.c_str() );
        smf2xx::readSMFHeader( smf, elementShape, elementNumPoints );
        smf.close();
    }
    
    // Input and output file streams
    std::ifstream smfIn(   smfFileIn.c_str() );

    // Call generic conversion helper
    smf2xx::Conversion< bbox::Analysis >::apply( elementShape,
                                                 elementNumPoints,
                                                 smfIn, std::cout );

    // Close the streams
    smfIn.close();
    
    return 0;
}
