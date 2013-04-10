//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smf2gp.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <fstream>
#include <string>
// boost includes
#include <boost/lexical_cast.hpp>
// base/mesh includes
#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
// base/io includes
#include <base/io/smf/Reader.hpp>
#include <base/io/gp/Writer.hpp>
// tools/converter/smf2xx includes
#include <tools/converter/smf2xx/SMFHeader.hpp>
#include <tools/converter/smf2xx/Conversion.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace smf2gp{
    
            //------------------------------------------------------------------
            /** Read smf file and write a gnuplot data file
             *  \tparam SHAPE  Type of element shape
             *  \tparam DEGREE Polynomial degree of the element
             */
            template<base::Shape SHAPE,unsigned DEGREE>
            struct Converter
            {
                static void apply( std::istream& smf,
                                   std::ostream& gp )
                {
                    // Attributes of the mesh
                    static const base::Shape shape   = SHAPE;
                    static const unsigned degree     = DEGREE;
                    static const unsigned    dim     = 3;
        
                    // Typedefs for defining a mesh
                    typedef base::mesh::Node<dim>                 Node;
                    typedef base::LagrangeShapeFun<degree,shape>  SFun;
                    typedef base::mesh::Element<Node,SFun>        Element;
                    typedef base::mesh::Unstructured<Element>     Mesh;

                    // Mesh object
                    Mesh mesh;

                    // SMF input
                    base::io::smf::Reader<Mesh> smfReader;
                    smfReader( mesh, smf ); 

                    // Gnuplot data output
                    base::io::gp::Writer::apply( mesh.elementsBegin(),
                                                 mesh.elementsEnd(), gp );
                }
            };
            
        } // namespace smf2gp
    } // namespace converter
} // namespace tools

//------------------------------------------------------------------------------
/** Read smf formatted file, create a temporary mesh and write gnuplot data file
 */
int main( int argc, char * argv[] )
{
    namespace smf2gp    = tools::converter::smf2gp;
    namespace smf2xx    = tools::converter::smf2xx;

    // Sanity check of the number of input arguments
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0]
                  << " file.smf \n\n";
        return 0;
    }

    // Name of smf input file, its basename and the data output file name
    const std::string smfFile = boost::lexical_cast<std::string>( argv[1] );
    const std::string base    = smfFile.substr(0, smfFile.find( ".smf") );
    const std::string gpFile  = base + ".dat";

    // Element attributes
    base::Shape elementShape;
    unsigned    elementNumPoints;
    
    {
        // extract data from header
        std::ifstream smf( smfFile.c_str() );
        smf2xx::readSMFHeader( smf, elementShape, elementNumPoints );
        smf.close();
    }

    // Input and output file streams
    std::ifstream smf( smfFile.c_str() );
    std::ofstream gp(  gpFile.c_str() );

    // Call generic conversion helper
    smf2xx::Conversion< smf2gp::Converter >::apply( elementShape,
                                                    elementNumPoints,
                                                    smf, gp );

    // Close the streams
    smf.close();
    gp.close();
    
    return 0;
}
