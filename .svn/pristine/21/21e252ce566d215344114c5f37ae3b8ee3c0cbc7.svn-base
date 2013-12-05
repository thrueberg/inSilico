//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smfOutput.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_converter_gmsh2smf_smfoutput_hpp
#define tools_converter_gmsh2smf_smfoutput_hpp

// std includes
#include <string>
#include <vector>
#include <iomanip>

// local includes
#include <tools/converter/gmsh2smf/elementTypes.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace gmsh2smf{

            //! Write coordinates to file with given name
            void writeCoordinates( const std::string& coordFile,
                                   const std::vector<gmsh2smf::Node>& nodes );


            //! Write connectivity and smf files
            void writeTopology(  const std::string& outFileBaseName,
                                 const std::map< unsigned, unsigned >& elementMap,
                                 const unsigned domNum,
                                 const std::string coordFile,
                                 const std::size_t numNodes,
                                 const std::vector<unsigned>& elementTypes,
                                 const std::vector<unsigned>& elementFirstTags,
                                 const std::vector< std::vector<std::size_t> >& connectivities );
        }
    }
}

//------------------------------------------------------------------------------
void tools::converter::gmsh2smf::writeCoordinates( const std::string& coordFile,
                                                   const std::vector<gmsh2smf::Node>& nodes )
{
    std::cout << "Found " << nodes.size() << " nodes\n";
    std::cout << "Writing coordinates to with " << coordFile << " ...\n";
    std::ofstream coord( coordFile.c_str() );
    for ( unsigned n = 0; n < nodes.size(); n ++ ) {
        coord << std::setprecision( std::numeric_limits<double>::digits10 )
              << nodes[n][0] << " " << nodes[n][1] << " "
              << nodes[n][2] << "\n";
    }
    coord.close();

    return;
}

//------------------------------------------------------------------------------
void tools::converter::gmsh2smf::writeTopology(  const std::string& outFileBaseName,
                                                 const std::map< unsigned, unsigned >& elementMap,
                                                 const unsigned domNum,
                                                 const std::string coordFile,
                                                 const std::size_t numNodes,
                                                 const std::vector<unsigned>& elementTypes,
                                                 const std::vector<unsigned>& elementFirstTags,
                                                 const std::vector< std::vector<std::size_t> >&
                                                 connectivities )
{

    //! go through all types and element of given domain
    for ( std::map<unsigned,unsigned>::const_iterator
              iter = elementMap.begin();
          iter != elementMap.end(); ++iter ) {

        //! Type of element
        const unsigned eType         = iter -> first;
        const unsigned numOfThisType = iter -> second;

        std::cout << "      " << numOfThisType << " Elements of type ("
                  << smfNameOfElementType( eType ) << ", "
                  << numNodesPerElement(   eType )
                  << ")\n";

        //! Add element type to file name if necessary
        const std::string extendedFileName =
            outFileBaseName +
            ( (elementMap.size() > 1 ) ? ("." + smfNameOfElementType(eType) ) : "" );

        //! Connectivity file name
        const std::string elementFileName = extendedFileName + ".conn";

        //! Output stream for connectivity
        std::ofstream conn( elementFileName.c_str() );

        //! Go through all elements
        for ( unsigned e = 0; e < elementTypes.size(); e ++ ) {
            
            //! If elment has the right domain and type number
            if ( ( elementFirstTags[e] == domNum ) and
                 ( elementTypes[e]     == eType ) ){

                //! Extract element
                const std::vector<std::size_t> element =
                    connectivities[e];
                
                //! Write connectivity
                std::copy( element.begin(), element.end(),
                           std::ostream_iterator<std::size_t>( conn, " ") );
                conn << "\n";
            }
            
        }
        conn.close();

        //! smf file (header only version)
        const std::string smfFileName = extendedFileName + ".smf";
        std::ofstream smf( smfFileName.c_str() );
        smf << "! elementShape " << smfNameOfElementType(eType) << "\n"
            << "! elementNumPoints " << numNodesPerElement(eType) << "\n"
            << "! externalNodes " << coordFile << "\n"
            << "! externalElements " << elementFileName << "\n"
            << "  " << numNodes << "  " << numOfThisType
            << "\n";
        smf.close();

    }
    return;
}

#endif
