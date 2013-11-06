//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   xdmf/Writer.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_xdmf_writer_hpp
#define base_io_xdmf_writer_hpp

//------------------------------------------------------------------------------
// std includes
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
// boost includes
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
// base/io includes
#include <base/io/OStreamIterator.hpp>
// base/io/raw includes
#include <base/io/raw/ascii.hpp>
// base/io/xml includes
#include <base/io/xml/Tag.hpp>
#include <base/io/xml/CData.hpp>
// base/io/xdmf includes
#include <base/io/xdmf/ElementTraits.hpp>
#include <base/io/xdmf/DoFTraits.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace xdmf{
            class Writer;
        }
    }
}

//------------------------------------------------------------------------------
/**
 *  see  http://www.xdmf.org/index.php/XDMF_Model_and_Format
 */
class base::io::xdmf::Writer
    : public boost::noncopyable
{
public:
    typedef base::io::xml::Tag                       Tag;
    typedef base::io::xml::CData                     CData;
    typedef base::io::xml::SimpleCData<std::string>  StringCData;

    //! Constructor initialises the main tags
    Writer()
        : xdmf_( "Xdmf" ), domain_( "Domain" ), grid_( "Grid" )
    { }

    //--------------------------------------------------------------------------
    //! Destructor
    ~Writer()
    {
        // For reason unknown to me the second, commented version produces
        // memory leaks !?!
        
        for ( unsigned s = 0; s < cDataPtrs_.size(); s ++ )
            delete cDataPtrs_[s];
        
        //std::for_each( cDataPtrs_.begin(), cDataPtrs_.end(),
        //               boost::bind( &operator delete, _1 ) );
    }

    //! Write coordinates to file and register xdmf-tags
    template<typename MESH> 
    void writeGeomtry( const MESH & mesh, const std::string & geomFile );

    //! Write connectivity to file and register xdmf-tags
    template<typename MESH> 
    void writeTopolgy( const MESH & mesh, const std::string & topoFile );
    
    //! Write attribute
    template<typename VALITER> 
    void writeAttribute( VALITER first, VALITER last,
                         const std::string& attrFile,
                         const std::string& name );


    void finalise( std::ostream & out ) 
    {
        //! Write header
        out << "<?xml version=\"1.0\" ?>\n"
            << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";

        xdmf_.addAttribute( "Version",  "2.0" );
        xdmf_.addAttribute( "xmlns:xi", "http://www.w3.org/2001/XInclude" );

        //Possibly create a nice ID for the grid here 
        grid_.addAttribute( "Name", "grid" );
        domain_.addTag( grid_ );
        xdmf_.addTag( domain_ );

        xdmf_.write( out );
        
    }
    
private:
    //! The main tag (everything else is a child)
    Tag  xdmf_;
    //! Domain tag
    Tag  domain_;
    //! Grid tag
    Tag grid_;
        
    //! Associate with a tag-name (string) the CData to be written with it
    std::vector<CData*> cDataPtrs_;
};

//------------------------------------------------------------------------------
//! Write the geometry of the mesh
template<typename MESH>
void base::io::xdmf::Writer::writeGeomtry( const MESH & mesh,
                                           const std::string & geomFile )
{
    std::ofstream geom( geomFile.c_str() );

    //! Go through all nodes and write their coordinates to a file
    typename MESH::NodePtrConstIter node = mesh.nodesBegin();
    typename MESH::NodePtrConstIter end  = mesh.nodesEnd();
    const std::size_t numNodes = std::distance( node, end );

    typedef typename MESH::Node Node;
    std::copy( node, end, 
               base::io::OStreamIterator<Node*>(
                   geom,
                   &base::io::raw::template writeNodeCoordinates<Node> ) );

    geom.close();


    //! Create geometry tag
    Tag geometry( "Geometry" );
    {
        geometry.addAttribute( "GeometryType", "XYZ" );
        geometry.addAttribute( "Dimensions",   numNodes );

        //! Create data item tag for reference
        Tag dataItem( "DataItem" );
        {
            dataItem.addAttribute( "Reference", "XML" );
            
            const std::string cdataValue( "/Xdmf/Domain/DataItem[@Name=\"Coordinates\"]" );
            CData * cdata = new StringCData( cdataValue );

            cDataPtrs_.push_back( cdata );
            
            dataItem.setCData( cdata );
        }
        geometry.addTag( dataItem );
    }
    grid_.addTag( geometry );


    //! Create a data item tag for the external coordinate file
    Tag dataItem( "DataItem" );
    {
        dataItem.addAttribute( "Name",       "Coordinates" );
        dataItem.addAttribute( "Format",     "XML" );
        dataItem.addAttribute( "NumberType", "Float" );
        const std::string dimName
            = boost::lexical_cast<std::string>( numNodes ) + " 3";
        dataItem.addAttribute( "Dimensions", dimName );

        //! Child tag with file name
        Tag child( "xi:include" );
        {
            child.addAttribute( "parse", "text" );
            child.addAttribute( "href",  geomFile );
        }
        dataItem.addTag( child );
    }
    domain_.addTag( dataItem );
}

//------------------------------------------------------------------------------
template<typename MESH> 
void base::io::xdmf::Writer::writeTopolgy( const MESH & mesh,
                                           const std::string & topoFile )
{
    std::ofstream topo( topoFile.c_str() );

    //! Element type
    typedef typename MESH::Element Element;

    //! Get basic details from ElementTraits
    typedef base::io::xdmf::ElementName<Element::shape,
                                        Element::numNodes> EN;
    typedef base::io::xdmf::ElementNumOutputNodes<Element::shape,
                                                  Element::numNodes> ENON;

    //! Go through all elements
    typename MESH::ElementPtrConstIter element = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter end     = mesh.elementsEnd();

    const std::size_t numElements = std::distance( element, end );

    //! Write QUAD and HEX as serendipity elements
    typedef base::io::xdmf::ElementNumOutputNodes<MESH::Element::shape,
                                                  MESH::Element::numNodes> ENON;

    //! Copy element connectivity to ostream iterator
    std::copy( element, end, 
               base::io::OStreamIterator<Element*>(
                   topo,
                   &base::io::raw::template writeElementConnectivity<Element,
                                                                     ENON::value> ) );

    topo.close();
    
    //! Create topology tag
    Tag topology( "Topology" );
    {
        topology.addAttribute( "TopologyType",     EN::value() );
        topology.addAttribute( "NumberOfElements", numElements );
        //! Create data item tag for reference
        Tag dataItem( "DataItem" );
        {
            dataItem.addAttribute( "Reference", "XML" );
            const std::string cdataValue( "/Xdmf/Domain/DataItem[@Name=\"Connectivity\"]" );

            CData * cdata = new StringCData( cdataValue );
            cDataPtrs_.push_back( cdata );
            
            dataItem.setCData( cdata );
        }
        topology.addTag( dataItem );
    }
    grid_.addTag( topology );

    //! Create a data item tag for the external coordinate file
    Tag dataItem( "DataItem" );
    {
        dataItem.addAttribute( "Name",       "Connectivity" );
        dataItem.addAttribute( "Format",     "XML" );
        dataItem.addAttribute( "NumberType", "Int" );
        const std::string dimName =
            boost::lexical_cast<std::string>(  numElements ) + " " +
            boost::lexical_cast<std::string>( +ENON::value );
        dataItem.addAttribute( "Dimensions", dimName );

        //! Child tag with file name
        Tag child( "xi:include" );
        {
            child.addAttribute( "parse", "text" );
            child.addAttribute( "href",  topoFile );
        }
        dataItem.addTag( child );
    }
    domain_.addTag( dataItem );
    
}

//------------------------------------------------------------------------------
template<typename VALITER> 
void base::io::xdmf::Writer::writeAttribute( VALITER first,
                                             VALITER last,
                                             const std::string& attrFile,
                                             const std::string& name )
{
    std::ofstream attr( attrFile.c_str() );

    const std::size_t numDoFs = std::distance( first, last );

    typedef typename std::iterator_traits<VALITER>::value_type DoFValue;

    std::copy( first, last, 
               base::io::OStreamIterator<DoFValue>(
                   attr,
                   &base::io::raw::template writeDoFValues<DoFValue> ) );
    
    attr.close();

    const unsigned size = first -> size();

    //! Create attribute tag
    Tag attribute( "Attribute" );
    {
        const std::string centerName = "Node";
        const std::string dofTypeName = ( size == 1 ? "Scalar" : "Vector" );
        
        attribute.addAttribute( "Name",          name );
        attribute.addAttribute( "Dimensions",    numDoFs );
        attribute.addAttribute( "AttributeType", dofTypeName );
        attribute.addAttribute( "Center",        centerName );
        
        //! Create data item tag for reference
        Tag dataItem( "DataItem" );
        {
            dataItem.addAttribute( "Reference", "XML" );
            const std::string cdataValue = 
                "/Xdmf/Domain/DataItem[@Name=\"" + name + "\"]";
        
            CData * cdata = new StringCData( cdataValue );
            cDataPtrs_.push_back( cdata );
            
            dataItem.setCData( cdata );
        }
        attribute.addTag( dataItem );
    }
    grid_.addTag( attribute );

    //! Create a data item tag for the external coordinate file
    Tag dataItem( "DataItem" );
    {
        dataItem.addAttribute( "Name",       name );
        dataItem.addAttribute( "Format",     "XML" );
        dataItem.addAttribute( "NumberType", "Float" );
        const std::string dimName =
            boost::lexical_cast<std::string>(  numDoFs ) + " " +
            boost::lexical_cast<std::string>(  size );
        dataItem.addAttribute( "Dimensions", dimName );

        //! Child tag with file name
        Tag child( "xi:include" );
        {
            child.addAttribute( "parse", "text" );
            child.addAttribute( "href",  attrFile );
        }
        dataItem.addTag( child );
    }
    domain_.addTag( dataItem );

}

#endif
