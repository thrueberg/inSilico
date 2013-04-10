//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   gmsh2smf.cpp
//! @author Thomas Rueberg
//! @date   2013

// std includes
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <map>

// boost includes
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

// local includes
#include <tools/converter/gmsh2smf/elementTypes.hpp>
#include <tools/converter/gmsh2smf/mshInput.hpp>
#include <tools/converter/gmsh2smf/smfOutput.hpp>


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    namespace gmsh2smf = tools::converter::gmsh2smf;
    
    if ( argc != 2 ) {
        std::cerr << " Usage: " << argv[0] << "  file.msh \n";
        exit(-1);
    }

    //! Mandatory input argument: Name of msh file (ends in .msh) 
    const std::string mshFile = boost::lexical_cast<std::string>( argv[1] );

    if ( mshFile.find( ".msh" ) == std::string::npos) {
        std::cerr << "(EE) Failed to find file ending \".msh\" \n";
        exit(-1);
    }

    //! Base name of input file determines base names of output
    const std::string baseName = mshFile.substr( 0, mshFile.find( ".msh" ) );

    //! Input stream
    std::ifstream msh( mshFile.c_str() );

    //! Validate stream
    if ( not msh.is_open() ) {
        std::cerr << "(EE) Could not open file " << mshFile << "\n";
        exit(-1);
    }

    //! @name Data reprenting the a msh-file
    //@{
    std::vector<std::string> physicalNames;
    physicalNames.push_back( "(unnamed)" );
    std::vector<gmsh2smf::Node> nodes;
    std::vector<unsigned> elementTypes, elementFirstTags;
    std::vector<std::vector<unsigned> > connectivities;
    //@}

    //! Create a tag reader to manage the input
    gmsh2smf::TagReader tagReader;

    //! Read until a complement mesh has been found
    bool carryOn = true; unsigned numTags = 0;
    while ( carryOn ) {

        const gmsh2smf::TagReader::Tag tag = tagReader.readAndValidate( msh );

        if ( tag == gmsh2smf::TagReader::MeshFormat ) {
            double version;
            unsigned dummy;
            msh >> version >> dummy >> dummy;
            msh.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        }
        else if ( tag == gmsh2smf::TagReader::PhysicalNames ) {
            gmsh2smf::readPhysicalNames( msh, physicalNames );
        }
        else if ( tag == gmsh2smf::TagReader::Nodes ) {
            gmsh2smf::readNodes( msh, nodes );
        }
        else if ( tag == gmsh2smf::TagReader::Elements ) {
            gmsh2smf::readElements( msh, elementTypes, elementFirstTags, connectivities );
            assert( elementTypes.size() == elementFirstTags.size() );
            assert( elementTypes.size() == connectivities.size() );
        }

        //! Read closing tag and check if complete mesh has been found
        carryOn = not ( tagReader.readEndTag( msh ) );

        // Prevent infinite loop in malicious msh file
        if ( ( numTags++ > 5 ) and carryOn ) {
            std::cerr << "(EE) Stuck in reading loop, file corrupt\n";
            exit(-1);
        }
        
    }
    msh.close();
        
    //--------------------------------------------------------------------------
    //! Output section

    //! write coordinates
    const std::string coordFile = baseName + ".coords";
    gmsh2smf::writeCoordinates( coordFile, nodes );

    //! tell caller what has been found
    std::cout << "----------------------------------------------------\n";
    std::cout << "Found " << physicalNames.size() << " physical names\n";
    std::cout << "Found " << elementTypes.size() << " elements\n";

    // map between physical name (or default) and a map of element type
    std::vector< std::map<unsigned, unsigned> >  elementMap;
    for ( unsigned d = 0; d < physicalNames.size(); d ++ ) {

        // count number of elements per type apprearing in this domain
        std::map< unsigned, unsigned > types;
        for ( unsigned e = 0; e < elementTypes.size(); e++ ) {
            const unsigned eType = elementTypes[e];
                
            if ( elementFirstTags[e] == d ) {
                if ( types.find( eType ) == types.end() )
                    types[ eType ] = 1;
                else
                    types[ eType ] += 1;
                    
            }
        }
        elementMap.push_back( types );
    }

    //! Write topology per domain
    std::cout << "----------------------------------------------------\n";
    for ( unsigned d = 0; d < elementMap.size(); d ++ ) {
        if ( elementMap[d].size() ) {

            const std::string outFileBaseName =
                baseName + ( d == 0 ? "" : ("."+ physicalNames[d]) );

            std::cout << "Domain " << d << " with name "
                      << physicalNames[d] << " has \n";

            gmsh2smf::writeTopology( outFileBaseName, elementMap[d], d,
                                     coordFile, nodes.size(),
                                     elementTypes, elementFirstTags,
                                     connectivities );
                
        }
            
    }


        
    return 0;
}
