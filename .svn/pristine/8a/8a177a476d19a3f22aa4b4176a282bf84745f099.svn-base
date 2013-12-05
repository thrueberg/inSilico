//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   mshInput.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_converter_gmsh2smf_mshinput_hpp
#define tools_converter_gmsh2smf_mshinput_hpp

// std includes
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <iomanip>
// local includes
#include <tools/converter/gmsh2smf/elementTypes.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace gmsh2smf{

            //! Read section of physical names
            void readPhysicalNames( std::istream& inp,
                                    std::vector<std::string>& physicalNames );

            //! Read section of nodes
            void readNodes( std::istream& inp,
                            std::vector<Node>& nodes );

            //! Read sectio of elements
            void readElements( std::istream& inp,
                               std::vector<unsigned>&                  elementTypes,
                               std::vector<unsigned>&                  elementFirstTags, 
                               std::vector<std::vector<std::size_t> >& connectivities );

            //! Manage the tags and their state of being read
            class TagReader;

        }
    }
}

//------------------------------------------------------------------------------
void tools::converter::gmsh2smf::readPhysicalNames( std::istream& inp,
                                                    std::vector<std::string>& physicalNames )
{
    unsigned numPhysicalNames;
    inp >> numPhysicalNames;
    inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    
    for ( unsigned p = 0; p < numPhysicalNames; p ++ ) {

        if ( inp.peek() == '$' ) {
            std::cerr << "(EE) Line begins with $ character \n";
            exit(-1);
        }
        
        unsigned dim, id;
        std::string name;
        inp >> dim >> id >> name;
        inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

        // remove quotation marks
        name.erase( 0, 1 );
        name.erase( name.length()-1, 1 );

        // store name
        physicalNames.push_back( name );

    }
}

//------------------------------------------------------------------------------
//! Function to read numNodes Nodes from the input stream inp
void  tools::converter::gmsh2smf::readNodes( std::istream& inp,
                                             std::vector<Node>& nodes )
{
    unsigned numNodes;
    inp >> numNodes;
    inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
                
    for ( unsigned n = 0; n < numNodes; n ++ ) {

        if ( inp.peek() == '$' ) {
            std::cerr << "(EE) Line begins with $ character \n";
            exit(-1);
        }

        unsigned id;
        double x, y, z;
        inp >> id;
        inp >> std::setprecision( std::numeric_limits<double>::digits10 )
            >> x >> y >> z;
        inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

        // sanity check
        if ( (n+1) != id ) {
            std::cerr << "(EE) Nodes are not numbered consecutively\n";
            exit(-1);
        }

        // create and store node
        Node node = {{ x, y, z }};
        nodes.push_back( node );
    }
    
}

//------------------------------------------------------------------------------
//! Function to read element types, tags and connectivities
void tools::converter::gmsh2smf::readElements( std::istream& inp,
                                               std::vector<unsigned>& elementTypes,
                                               std::vector<unsigned>& elementFirstTags, 
                                               std::vector<std::vector<std::size_t> >&
                                               connectivities )
{
    unsigned numElements;
    inp >> numElements;
    inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
                   
    for ( unsigned e = 0; e < numElements; e ++ ) {

        if ( inp.peek() == '$' ) {
            std::cerr << "(EE) Line begins with $ character \n";
            exit(-1);
        }

        unsigned  id, type, numTags, firstTag, dummy;
        inp >> id >> type >> numTags >> firstTag;

        // ignore all other tags
        for ( unsigned i = 0; i < numTags-1; i ++) inp >> dummy;

        // number of nodes depends on element type
        const unsigned numNodes = gmsh2smf::numNodesPerElement( type );

        // read connectivity
        std::vector<std::size_t> connec;
        for ( unsigned n = 0; n < numNodes; n++ ) {
            std::size_t nodeId;
            inp >> nodeId;
            connec.push_back( nodeId-1 );

        }
        
        gmsh2smf::reorderConnectivity( type, connec );


        inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

        // store data
        elementTypes.push_back(     type );
        elementFirstTags.push_back( firstTag );
        connectivities.push_back(    connec );
    }

    return;
}

//------------------------------------------------------------------------------
class tools::converter::gmsh2smf::TagReader
{
public:
    //! Possible tags to read
    enum Tag{ INVALID = -1, MeshFormat, PhysicalNames, Nodes, Elements };

private:
    //! State flags
    enum State { OPEN, READ, CLOSED };

public:
    //! Initialise arrays
    TagReader()
    {
        tagNames_[0] = "MeshFormat";
        tagNames_[1] = "PhysicalNames";
        tagNames_[2] = "Nodes";
        tagNames_[3] = "Elements";

        tagStates_.assign( OPEN );
    }

private:
    //! Read line and extract the tag name
    std::string readTagline_( std::istream& inp ) const
    {
        if ( not (inp.peek() == '$') ) {
            std::cerr << "(EE) Line does not contain a tag \n ";
            exit(-1);
        }

        std::string word;
        inp >> word;
        inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        word.erase( 0, 1 ); // delete first character

        return word;
    }

public:
    //! Read a tag, validate it, and return the corresponding Tag enum
    Tag readAndValidate( std::istream& inp )
    {
        const std::string word = this -> readTagline_( inp );

        // check against existing tags
        int tagNum = -1;
        std::size_t pos;
        for ( unsigned t = 0; t < tagNames_.size(); t++ ) {
            const std::string tagFromList = tagNames_[t];

            pos = word.find( tagFromList );
            if ( pos != std::string::npos ) {
                tagNum = t;
                break;
            }
        }

        // Check if known tag has been found
        if ( tagNum < 0 ) {
            std::cerr << "(EE) Found tag \"" << word << "\" which is not in list\n";
            exit(-1);
        }

        // Check that tag is labelled as open
        if ( tagStates_.at( tagNum ) != OPEN ) {
            std::cerr << "(EE) Tag \""<< word << "\" is not OPEN \n";
            exit(-1);
        }

        // Check that tag has been found at the beginning of string
        if( pos != 0 ) {
            std::cerr << "(EE) Found tag " << word << " and expected  string \""
                      << tagNames_[tagNum] << "\" at position 0 \n";
            exit( -1 );
        }

        tagStates_.at(tagNum) = READ;
        
        return static_cast<Tag>( tagNum );
    }

    //! Read the end-tag corresponding to the last read tag
    bool readEndTag( std::istream& inp )
    {
        const std::string word = this -> readTagline_( inp );

        const std::size_t numReadTags =
            std::count_if( tagStates_.begin(), tagStates_.end(),
                           boost::bind( std::equal_to<State>(), _1, READ ) );

        if ( numReadTags != 1 ) {
            std::cerr << "(EE) Found " << numReadTags << " tags to be in "
                      << "unfinished state\n";
            exit(-1);
        }

        const  boost::array<State,4>::iterator iter = std::find( tagStates_.begin(),
                                                                 tagStates_.end(), READ );
        const std::size_t readTagNum = std::distance( tagStates_.begin(), iter );

        const std::string expectedWord = "End" + tagNames_[ readTagNum ];

        if ( expectedWord != word ) {
            std::cerr << "(EE) Extected an end tag \"" << expectedWord
                      << "\" but found instead \"" << word << "\" \n";
        }

        tagStates_.at( readTagNum ) = CLOSED;

        //! Check if complete mesh has been read (note PhysicalNames are optional)
        bool foundCompleteMesh = (
            (tagStates_[0] == CLOSED) and  // read mesh format
            (tagStates_[2] == CLOSED) and  // read nodes
            (tagStates_[3] == CLOSED) );   // read elements

        return foundCompleteMesh;
    }

private:
    boost::array< std::string, 4 > tagNames_;
    boost::array< State, 4 >       tagStates_;
};

#endif
