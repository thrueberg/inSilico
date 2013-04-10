//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   CData.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_xml_cdata_hpp
#define base_io_xml_cdata_hpp

//------------------------------------------------------------------------------
// std  includes
#include <ostream>
// boost includes
#include <boost/utility.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace xml{

            class CData;

            template<typename T>
            class SimpleCData;
        }
    }
}


//------------------------------------------------------------------------------
/** Abstract base class for the typed CDATA
 *  Using this object, the Tag can simple store a pointer to it without being
 *  aware of the underlying complexity.
 */
class base::io::xml::CData : public boost::noncopyable
{
public:
    virtual ~CData() { }
    virtual void write( std::ostream & out ) const = 0;
};

//------------------------------------------------------------------------------
/** Specialised CData object for simple types which have an 'operator<<'.
 *  If you need anything more complex (e.g. an array of numbers), you should
 *  provide your own object which inherits from  CData and provides the correct
 *  write function
 *  \tparam T  Type of datum with operator<<
 */
template<typename T>
class base::io::xml::SimpleCData
    : public base::io::xml::CData
{
public:
    SimpleCData( const T & t )
    : t_( t ) { }
                
    void write( std::ostream & out ) const
    {
        out << t_;
    }
                
private:
    const T t_;
};
                

#endif
