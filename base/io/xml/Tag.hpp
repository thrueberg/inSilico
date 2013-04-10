//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Tag.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_xml_tag_hpp
#define base_io_xml_tag_hpp

//------------------------------------------------------------------------------
// std includes
#include <string>
#include <vector>
#include <utility>
#include <sstream>
// boost includes
#include <boost/lexical_cast.hpp>
// base includes
#include <base/verify.hpp>
// base/io/XML includes
#include <base/io/xml/CData.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace xml{

            class Tag;
        }
    }
}

//------------------------------------------------------------------------------
/** Representation of an XML tag with attributes and children or CData
 *
 *  This object represents and writes itself as an XML tag of the form
 *  \code{.xml}
 *        <Name  Attr1="Val1"  Attr2="Val2" ... />
 *  \endcode
 *  in case of simple tags,
 *  \code{.xml}
 *        <Name  Attr1="Val1"  Attr2="Val2" ... >
 *          CDATA
 *        </Name>
 *  \endcode
 *  in case of a tag embracing CDATA (which are given as strings onlye), and
 *  \code{.xml}
 *        <Name  Attr1="Val1"  Attr2="Val2" ... >
 *           ... (children tags)
 *        </Name>
 *  \endcode
 *  where the tag has a list of other tags as children. The children themselves
 *  are again of one of these three forms.
 *  In the writing, indentation is carried out according to the children level. 
 */
class base::io::xml::Tag
{
public:
    //! Constructor with name
    Tag( const std::string & name )
        : name_( name ), cData_( NULL )
    { }

    //--------------------------------------------------------------------------
    /** Store an attribute of the tag with name and value.
     *  The given value will be converted to a string.
     *  \param[in] attrName  Name of the attribute
     *  \param[in] attrValue Value of the attribute
     *  \tparam Type of attribute (internally converted to a string)
     */
    template<typename T>
    void addAttribute( const std::string& attrName,
                       const T&           attrValue )
    {
        const std::string valueString =
            boost::lexical_cast<std::string>( attrValue );

        attributes_.push_back( std::make_pair( attrName, valueString ) );
    }

    //--------------------------------------------------------------------------
    /** Store a tag as child.
     *  \param[in] tag  Child tag
     */
    void addTag( Tag & tag )
    {
        children_.push_back( tag );
    }

    //--------------------------------------------------------------------------
    /** Store a reference to a object hidden in the CData class
     */
    void setCData( CData * cData )
    {
        cData_ = cData;
    }

    //--------------------------------------------------------------------------
    /** Write itself and write either CData or child tags
     *  \param[in] out         Stream for output
     *  \param[in] indentLevel Level of indentation for legibility
     */
    void write( std::ostream & out, const unsigned indentLevel = 0 ) const
    {
        std::string indent( 2 * indentLevel, ' ' );
        
        out << indent << "<" << name_ << " ";

        // write attributes;
        for ( unsigned a = 0; a < attributes_.size(); a ++ ) {
            out << attributes_[a].first  << "=\""
                << attributes_[a].second << "\"";

            if ( a < (attributes_.size() - 1) ) out << " ";
        }

        // check children
        const bool noChildren = children_.empty();

        // check cdata_
        const bool noCData    = (cData_ == NULL);

        // closing symbol
        if ( noChildren and noCData ) {
            out << "/> \n";
        }
        else out << "> \n";

        // sanity check
        VERIFY_MSG( (noCData or noChildren),
                    "Cannot have data AND children" );

        //! Write CDATA if available (not indented)
        if ( not noCData ) {
            cData_ -> write( out );
            out << "\n";
        }

        //! Write children if available
        if ( not noChildren ) {
            for ( unsigned c = 0; c < children_.size(); c ++ )
                children_[c].write( out, indentLevel+1 );
        }

        //! write closing tag if their had been data or children
        if ( not( noChildren and noCData ) ) {
            out << indent << "</" << name_ << ">\n";
        }
        
    }

private:
    //--------------------------------------------------------------------------
    /** The name of this tag, e.g.
     *  \code{.xml} <Name ..> ... </Name> \endcode
     *  would be the XML output
     */
    const std::string                                   name_;
    //! A list of attributes with names and values, the latter as strings
    std::vector< std::pair<std::string,std::string> >   attributes_;
    //! The children tags
    std::vector<Tag>                                    children_;
    //! Pointer to a CData object 
    CData*                                              cData_;
};


#endif
