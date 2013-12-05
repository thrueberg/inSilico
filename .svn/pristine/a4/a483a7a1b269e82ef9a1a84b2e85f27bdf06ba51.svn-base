//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   sgfRefine.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <fstream>
#include <string>
// boost includes
#include <boost/lexical_cast.hpp>
// base  includes
#include <base/Structured.hpp>
// base/io includes
#include <base/io/sgf/Reader.hpp>
#include <base/io/sgf/Writer.hpp>
#include <base/io/Format.hpp>
// tools/converter/sgf2xx includes
#include <tools/converter/sgf2xx/SGFHeader.hpp>
#include <tools/converter/sgf2xx/Conversion.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace sgfRefine{

            static unsigned refinement = 1;
    
            //------------------------------------------------------------------
            /** Read sgf file and write refined sgf file
             *  \tparam DIM    Dimension of the grid
             *  \tparam DEGREE Polynomial degree of geometry representation
             */
            template<unsigned DIM, unsigned DEGREE>
            struct Refiner
            {
                static void apply( std::istream& sgfIn,
                                   std::ostream& sgfOut )
                {
                    // Defining a grid
                    typedef base::Structured<DIM,DEGREE>        Grid;

                    // Grid object
                    Grid grid;

                    // SGF input
                    base::io::sgf::readGrid( sgfIn, grid );

                    // SGF refined output
                    base::io::sgf::writeGrid( grid, sgfOut, refinement );
                }

            };
            
        } // namespace sgfRefine
    } // namespace converter
} // namespace tools

//------------------------------------------------------------------------------
/** Read sgf formatted file, create a temporary grid and write a vtk file
 */
int main( int argc, char * argv[] )
{
    namespace sgfRefine = tools::converter::sgfRefine;
    namespace sgf2xx    = tools::converter::sgf2xx;

    // Sanity check of the number of input arguments
    if ( (argc != 3) and ( argc != 4) ) {
        std::cout << "Usage:  " << argv[0]
                  << " file.sgf refinement [degree=1] \n\n";
        return 0;
    }

    // Name of sgf input file, its basename and the vtk output file name
    const std::string sgfInFile  = boost::lexical_cast<std::string>( argv[1] );
    const unsigned refinement    = boost::lexical_cast<unsigned>(    argv[2] );
    const std::string base       = sgfInFile.substr(0, sgfInFile.find( ".sgf") );

    const std::string refString = base::io::leadingZeros( refinement, 2 );
    const std::string sgfOutFile = base + "." + refString + ".sgf";

    // Degree of geometry representation
    const unsigned degree = (argc==3) ? 1 : boost::lexical_cast<unsigned>( argv[3] );

    // Important set the resolution
    sgfRefine::refinement = refinement;

    // Grid dimension
    unsigned dim;
    {
        // extract data from header
        std::ifstream sgf( sgfInFile.c_str() );
        dim = sgf2xx::readDimFromSGFHeader( sgf );
        sgf.close();
    }

    // Input and output file streams
    std::ifstream sgfIn(   sgfInFile.c_str() );
    std::ofstream sgfOut( sgfOutFile.c_str() );

    // Call generic conversion helper
    sgf2xx::Conversion< sgfRefine::Refiner >::apply( dim, degree, sgfIn, sgfOut );

    // Close the streams
    sgfIn.close();
    sgfOut.close();
    
    return 0;
}
