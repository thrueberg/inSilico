//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SMFHeader.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef smf2xx_smfheader_hpp
#define smf2xx_smfheader_hpp

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <string>
// boost includes
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
// base includes
#include <base/verify.hpp>
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace smf2xx{

            //------------------------------------------------------------------
            //! Find shape value in a given string
            base::Shape findInString( const std::string& shapeString )
            {
    
                if      ( shapeString.find( "point" ) != std::string::npos )
                    return base::POINT;
                else if ( shapeString.find( "line"  ) != std::string::npos )
                    return base::LINE;
                else if ( shapeString.find( "tria"  ) != std::string::npos )
                    return base::TRI;
                else if ( shapeString.find( "quad"  ) != std::string::npos )
                    return base::QUAD;
                else if ( shapeString.find( "tet"   ) != std::string::npos )
                    return base::TET;
                else if ( shapeString.find( "hex"   ) != std::string::npos )
                    return base::HEX;

                VERIFY_MSG(  false, "Could not detect a shape descriptor" );

                return base::POINT;
            }

            //------------------------------------------------------------------
            //! Read SMF header and find values of shape and number of points
            void readSMFHeader( std::istream& smf,
                                base::Shape&  elementShape,
                                unsigned&     elementNumPoints )
            {
                // Skip leading comment lines
                const char commentChar = '#';
                while ( smf.peek() == commentChar )
                    smf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

                // Check if everything has been found
                bool foundShape  = false;
                bool foundNumber = false;

                // Look for header with element description
                const char headerChar = '!';
                while ( smf.peek() == headerChar )
                {
                    // Eat header character
                    char dummy;
                    smf >> dummy;

                    // Read rest of the line
                    std::string line;
                    getline( smf, line );

                    boost::tokenizer<> tokens( line );
                    boost::tokenizer<>::iterator iter = tokens.begin();
                    const std::string descriptor = *iter;

                    if      ( descriptor.find( "elementShape" ) != std::string::npos ) {
                        // Read shape from file
                        foundShape = true;
                        ++iter;
                        const std::string value = *iter;

                        // Convert to base::Shape
                        elementShape = smf2xx::findInString( value );
            
                    }
                    else if ( descriptor.find( "elementNumPoints" ) != std::string::npos ) {
            
                        // Read number of element nodes from file
                        foundNumber = true;
                        ++iter;
                        const std::string value = *iter;

                        // Convert to unsigned
                        elementNumPoints = boost::lexical_cast<unsigned>( value );
                    }
            
                }

                VERIFY_MSG( (foundShape and foundNumber),
                            "Could not find a valid SMF header" );
            }
            
        } // namespace smf2xx
    } // namespace converter
} // namespace tools
#endif
