//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   sgf2vtk.cpp
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
#include <base/mesh/Structured.hpp>
#include <base/sfun/BSpline.hpp>
// base/io includes
#include <base/io/sgf/Reader.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
// tools/converter/sgf2xx includes
#include <tools/converter/sgf2xx/SGFHeader.hpp>
#include <tools/converter/sgf2xx/Conversion.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace sgf2vtk{
    
            //------------------------------------------------------------------
            /** Read sgf file and write vtk file.
             *  \tparam DIM    Dimension of the grid
             *  \tparam DEGREE Degree of the geometry representation
             */
            template<unsigned DIM, unsigned DEGREE>
            struct Converter
            {
                static void apply( std::istream& sgf,
                                   std::ostream& vtk )
                {
                    // Typedefs for defining a mesh
                    typedef base::mesh::Node<DIM>                  Node;
                    typedef base::sfun::BSpline<DEGREE>            BSpline;
                    typedef base::sfun::TensorProduct<BSpline,DIM> GeomFun;
                    typedef base::mesh::Element<Node,GeomFun>      Element;
                    typedef base::mesh::Structured<Element>        Grid;

                    // Grid object
                    Grid grid;

                    // SGF input
                    base::io::sgf::readGrid( sgf, grid );
            
                    // VTK Legacy output
                    base::io::vtk::LegacyWriter vtkWriter( vtk );
                    vtkWriter.writeStructuredGrid( grid );
                }
            };
            
        } // namespace sgf2vtk
    } // namespace converter
} // namespace tools

//------------------------------------------------------------------------------
/** Read sgf formatted file, create a temporary grid and write a vtk file
 */
int main( int argc, char * argv[] )
{
    namespace sgf2vtk   = tools::converter::sgf2vtk;
    namespace sgf2xx    = tools::converter::sgf2xx;

    // Sanity check of the number of input arguments
    if ( (argc != 2) and (argc != 3) ) {
        std::cout << "Usage:  " << argv[0] << " file.sgf [degree=1] \n\n";
        return 0;
    }

    // Name of sgf input file, its basename and the vtk output file name
    const std::string sgfFile = boost::lexical_cast<std::string>( argv[1] );
    const std::string base    = sgfFile.substr(0, sgfFile.find( ".sgf") );
    const std::string vtkFile = base + ".vtk";

    // geomtry degree
    const unsigned degree =
        (argc == 2) ? 1 : boost::lexical_cast<unsigned>( argv[2] );

    // Grid dimension
    unsigned dim;
    {
        // extract data from header
        std::ifstream sgf( sgfFile.c_str() );
        dim = sgf2xx::readDimFromSGFHeader( sgf );
        sgf.close();
    }

    // Input and output file streams
    std::ifstream sgf( sgfFile.c_str() );
    std::ofstream vtk( vtkFile.c_str() );

    // Call generic conversion helper
    sgf2xx::Conversion< sgf2vtk::Converter >::apply( dim, degree, sgf, vtk );

    // Close the streams
    sgf.close();
    vtk.close();
    
    return 0;
}
