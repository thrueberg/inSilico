//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   codeExtractor.cpp
//! @author Thomas Rueberg
//! @date   2013

// std includes
#include <string>
#include <fstream>
#include <cassert>

// boost includes
#include <boost/lexical_cast.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/algorithm/string.hpp>

//------------------------------------------------------------------------------
namespace codeExtractor{
    
    // make boost::xpressive more concise
    namespace xpr = boost::xpressive;

    //--------------------------------------------------------------------------
    // Type of tag found in line of code
    enum TagType{
        CODE,        //!< no comment line
        BEGIN,       //!< begin-tag of code block to be extracted
        END,         //!< end-tag   of code block to be extracted
        ONELINE,     //!< tag for a one-line snippet
        COMMENTONLY  //!< tag that a plain comment line is found (no tags)
    };

    //--------------------------------------------------------------------------
    /** Look for code extraction tags in comment lines.
     *  Given a line of code, it will be checked if this line contains special
     *  code extraction tags.  The following options exist:
     *
     *     string     |   Meaning
     *  ------------- | ---------------------------------------------------------
     *  \c //         | standard comment line, nothing to be done (COMMENTONLY)
     *  \c //[tagID]  | the next line of code shall be extracted  (ONELINE)
     *  \c //[tagID]{ | this marks the begin of a code block      (BEGIN)
     *  \c //[tagID]} | this marks the end of a code block        (END)
     *     else       | The line is assumed to be a code line
     *
     *  Note that all other comment markers are not supported, e.g. the
     *  pairs of C-tyle comments (slash and asterisk).
     *
     *  \param[in]  line        Line of code
     *  \param[out] tagIDString String carrying the tag ID
     *  \returns                Type of tag
     */
    TagType findTagInLine( const std::string& line,
                           std::string& tagIDString )
    {
        // result value
        TagType tagType = CODE;

        // Make a trimmed copy
        // The motivation is to make sure that the comment symbols come as
        // the first characters in the line. Then, it is safe to work with
        // xpr::bos as 'begin-of-sequence' expression. This in turn is needed
        // in order to avoid faulty behaviour with post-fix comments.
        const std::string lineCopy = boost::trim_left_copy( line );

        // place holder for the tag ID
        xpr::mark_tag tagID(1);

        // regular expression of begin-block comment: //[tagID]{
        xpr::sregex beginBlock
            = xpr::bos >> xpr::as_xpr('/') >> '/'
                       >> '[' >> (tagID = +xpr::_w) >> ']' >> '{';

        // regular expression of begin-oneLiner comment: //[tagID]
        xpr::sregex oneLiner
            = xpr::bos >> xpr::as_xpr('/') >> '/'
                       >> '[' >> (tagID = +xpr::_w) >> ']'
                       >>  ( ~(xpr::set= '{', '}' ) | xpr::eos );

        // regular expression of end-block comment: //[tagID]}
        xpr::sregex endBlock
            = xpr::bos >> xpr::as_xpr('/') >> '/'
                       >> '[' >> (tagID = +xpr::_w) >> ']' >> '}';

        // regular expression for just a plain comment: //
        xpr::sregex justComment
            = xpr::bos >> xpr::as_xpr( '/' ) >> '/'
                       >> ~( xpr::as_xpr( '[' ) );

        // string match container
        xpr::smatch matches;

        // perform matching
        if ( xpr::regex_search( lineCopy, matches, beginBlock ) ){
            tagIDString = matches[ tagID ];
            tagType = BEGIN;
        }
        else if ( xpr::regex_search( lineCopy, matches, oneLiner ) ){
            tagIDString = matches[ tagID ];
            tagType = ONELINE;
        }
        else if ( xpr::regex_search( lineCopy, matches, endBlock ) ){
            tagIDString = matches[ tagID ];
            tagType = END;
        }
        else if ( xpr::regex_search( lineCopy, matches, justComment ) ){
            tagIDString = matches[ tagID ];
            tagType = COMMENTONLY;
        }

        return tagType;
    }

    //--------------------------------------------------------------------------
    /** Lines which are detected as code lines (see findTagInLine), can still
     *  contain post-fix comments. Since these are not desirable in the
     *  current context, they will be removed if detected.
     *  \param[in,out] line  Line of code to be freed of a post-fix comment
     */
    void trimPostFix( std::string& line )
    {
        // regular expression for a postfix comment //
        xpr::sregex comment = xpr::as_xpr( '/' )  >> '/';
        
        xpr::smatch matches;

        // if match is found remove the part of the line including and
        // following the comment symbol
        if ( xpr::regex_search( line, matches, comment ) ) {
            line = line.substr( 0, matches.position(0) );
        }
    }

} // end namespace codeExtractor

//------------------------------------------------------------------------------
/** This program is intended to be used for documenting the tutorial files.
 *  By applying it to a code file (.cpp), it creates file containing extracted
 *  code: one file with the complete code without comments and files according
 *  to tagged code snippets.
 *  These files can be included in the tutorial documentation as external code
 *  files.
 *
 *  A special convention is used here as described in
 *  codeExtractor::findTagInLine in order to mark code snippets. Any deviation
 *  from this rule will most likely not yield the desired result, in the worst
 *  case crash even. Moreover, consider the rules:
 *
 *  -  every begin marker must be followed an end marker
 *  -  every end marker must carry the same ID as its preceeding begin marker
 *  -  code snippets must not be nested
 *  -  neither a one-line nor a begin marker must be placed in the last line
 *     of the code
 *
 *  Any violation of these rules will probably lead to malicious behaviour.
 *
 *  \param[in] argc  Number of command line arguments (including program name)
 *  \param[in] argv  Value of the command line arguments
 *                   -  the first is the program name
 *                   -  the second should be the code file (.cpp)
 *                   -  the third is an optional output path
 *                   
 */
int main( int argc, char* argv[] )
{
    // make sure the right number if arguments is provided
    if ( (argc != 2) and (argc != 3) ) {
        std::cout << "Usage: " << argv[0] << "  file.cpp [outpath] \n\n"
                  << "  Creates file_plain.cpp (with the uncommented code) and\n"
                  << "  the files file_{id}.cpp (with snippets for every id) \n"
                  << "  by parsing the c++ code in file.cpp. The files are \n"
                  << "  placed in outpath wich defaults to the PWD.\n";
        return 0;
    }
    
    // first input command: the name of the code file
    const std::string cppFile = boost::lexical_cast<std::string>( argv[1] );

    // user-provided output path
    const bool withOutputPath = ( argc == 3 );
    const std::string outputPathName = ( withOutputPath ? 
                                         boost::lexical_cast<std::string>( argv[2] ) :
                                         "." );

    // make a director if the path requires that
    if ( withOutputPath ) {
        const std::string mkdirCmd = "mkdir -p " + outputPathName;
        const int check = system( mkdirCmd.c_str() );
        assert( check == 0 );
    }

    // extract its basename
    const std::string baseName = cppFile.substr( 0, cppFile.find( ".cpp" ) );
    
    // name of the plain-code file
    const std::string plainName = outputPathName + "/" + baseName + "_plain.cpp";

    // open input file 
    std::ifstream cpp( cppFile.c_str() );

    // sanity check
    assert( cpp.is_open() );

    // open output stream for plain code
    std::ofstream plain( plainName.c_str() );

    // some sort of worst-scenario catch
    const unsigned maxCodeBlockSize = 1000; // educated guess ;-)

    // as long as lines are found, read them
    std::string line;
    while ( std::getline(cpp, line) ) {

        // read tag and ID from code line
        std::string tagID;
        codeExtractor::TagType tagInLine =
            codeExtractor::findTagInLine( line, tagID );

        // go through the cases of possible tags
        if ( tagInLine == codeExtractor::CODE ) {

            // remove possible post-fix comments
            codeExtractor::trimPostFix( line );

            // plain code
            plain << line << '\n';
        }
        else if ( tagInLine == codeExtractor::ONELINE ) {

            // create file name and open the stream
            const std::string taggedCodeFile =
                outputPathName + "/" +
                baseName + "_" + tagID + ".cpp";
            
            std::ofstream tg( taggedCodeFile.c_str() );

            // just a line of code
            std::getline( cpp, line );
            // remove possible post-fix comments
            codeExtractor::trimPostFix( line );
            // write to plain code
            plain << line << '\n';
            // remove all leading whitespaces
            boost::trim_left( line );
            // write to tagged file
            tg    << line << '\n';
            tg.close();
            
        }
        else if ( tagInLine == codeExtractor::BEGIN ) {
            
            // create file name and open the stream
            const std::string taggedCodeFile =
                outputPathName + "/" +
                baseName + "_" + tagID + ".cpp";
            
            std::ofstream tg( taggedCodeFile.c_str() );

            // count number of leading whitespace for the first code line
            unsigned numWhiteSpaces = 0;
            unsigned codeLineCounter = 0;

            // found begin of a code to be extracted
            // read following lines until end tag
            while ( true ) {

                // get line of code and analyse it
                std::getline( cpp, line );
                std::string cmpTag;
                tagInLine = codeExtractor::findTagInLine( line, cmpTag );

                // possible tags in this code line
                if ( tagInLine == codeExtractor::END ) {
                    // verify the equality of the tag IDs
                    assert( cmpTag==tagID );
                    break;
                }
                else if ( tagInLine == codeExtractor::CODE ) {

                    // remove possible post-fix comments
                    codeExtractor::trimPostFix( line );

                    // write to plain code
                    plain << line << '\n';

                    // count number of leading whitespaces for first code line
                    if ( codeLineCounter == 0 ) {
                        // make trimmed copy
                        const std::string copy
                            = boost::trim_left_copy( line );
                        // difference in size = number of whitespaces
                        numWhiteSpaces = line.length() - copy.length();

                        // increment code line counter
                        codeLineCounter++;
                    }

                    // trim from left
                    if (line.length() > numWhiteSpaces ) {
                        // sanity check
                        const std::string copy
                            = boost::trim_left_copy( line );
                        // difference in size = number of whitespaces
                        const unsigned currentNumWhiteSpaces = line.length() - copy.length();

                        // make sure that not too much space is trimmed
                        // note this can look ugly, but the only correct solution
                        // would require to have a look ahead for all lines in the
                        // block, find the minimum number of leading white spaces
                        // and trim only by that number
                        if ( currentNumWhiteSpaces > numWhiteSpaces ) 
                            line = line.substr( numWhiteSpaces );
                        else
                            line = line.substr( currentNumWhiteSpaces );
                    }

                    // write to tagged code file (if not empty)
                    if ( boost::trim_left_copy( line ).length() )
                        tg    << line << '\n';

                }
                else if ( tagInLine == codeExtractor::COMMENTONLY ) {
                    // ordinary comment lines are just ignored here
                    ; // do nothing
                }
                else  // this case must not be reached
                    assert( false );

                // If the lines of this block exceed a certain number,
                // better do something about it
                if ( codeLineCounter > maxCodeBlockSize ) {
                    std::cerr << "(EE) Exceeded " << maxCodeBlockSize
                              << " lines of ocde in the block with ID "
                              << tagID << "\n"
                              << "(EE) Aborting \n";
                    exit(-1);
                        
                }

            }

            // end tag must be reached
            tg.close();
        }
        else if ( tagInLine == codeExtractor::COMMENTONLY ) {
            // ignore the plain code lines
            ; // just plain comments, do nothing

        }
        else {
            // this case is illegal
            assert( false );
        }

        
    } // done reading code lines
    
    cpp.close();
    plain.close();
    
    return 0;
}
