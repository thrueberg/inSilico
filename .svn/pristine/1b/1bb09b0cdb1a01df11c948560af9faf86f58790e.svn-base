//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SGFHeader.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef sgf2xx_sgfheader_hpp
#define sgf2xx_sgfheader_hpp

//------------------------------------------------------------------------------
// std includes
#include <istream>
// base includes
#include <base/verify.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace sgf2xx{
            //------------------------------------------------------------------
            //! Read SGF header, deduce the spatial dimension from the grid size
            unsigned readDimFromSGFHeader( std::istream& sgf )
            {
                // Skip leading comment lines
                const char commentChar = '#';
                while ( sgf.peek() == commentChar )
                    sgf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );


                // Read grid dimension
                boost::array<unsigned,3> gridDim;
                for ( unsigned d = 0; d < 3; d++ )
                    sgf >> gridDim[d];

                const unsigned dim =
                    ( gridDim[2] > 0 ? 3 : ( gridDim[1] > 0 ? 2 : 1 ) );
                VERIFY_MSG( gridDim[0] > 0, "Something strange happened" );

                return dim;
            }
            
        } // namespace sgf2xx
    } // namespace converter
} // namespace tools
#endif
