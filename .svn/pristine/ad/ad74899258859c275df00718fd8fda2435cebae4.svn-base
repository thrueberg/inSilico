//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   io/List.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_io_list_hpp
#define base_io_list_hpp

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <vector>
#include <limits>
// base includes
#include <base/verify.hpp>
#include <base/io/Format.hpp>
#include <base/io/PropertiesParser.hpp>

//------------------------------------------------------------------------------
// declarations
namespace base {
    namespace io {

        template<typename T> class List;
        
    }
}

//------------------------------------------------------------------------------
/**  Read and store a list of values of type 'T'.
 *   In combination with the base::io::PropertiesParser it is sometimes useful
 *   to read a variable number of values, for instance a list of file names.
 *   This object stores such a list an provides in- and output routines.
 *   The following snippet of a text file illustrates the functionality
 *   \code{.txt}
 *   variable     value
 *   # other variables ...
 *   variableList {
 *                  value1  # a comment
 *                  value2
 *
 *                  # another comment
 *                  value3  value4
 *                 }
 *
 *   variable2   anotherValue
 *   \endcode
 *   From this text file the variables 'variable' and 'variable2' are read
 *   as usual by the PropertiesParser. The variable 'variableList' is
 *   represented by this class and the values 'value1', 'value2', 'value3'
 *   and 'value4' will be stored here.
 *
 *   \note All values of the variable list are assumed to be of the same type.
 *         Any type with an istream and ostream operator is allowed.
 *         The symbols '{' and '}' are mandatory to identify the begin and
 *         end of the list. Whitespaces between these symbols and the values
 *         are necessary if the type of the variables is a string.
 *         Separation is done by a newline or any other whitespace character.
 *
 *  \tparam T  Type of variables to be stored in this list
 */
template<typename T>
class base::io::List
{
public:

    //! Constructor initialises the begin- and end-symbol
    List()
        : beginSymbol_( '{' ),
          endSymbol_(   '}' ) { }
    
    //! @name Accessors
    //@{
    typedef typename std::vector<T>::const_iterator ListIter;
    ListIter begin() const { return data_.begin(); }
    ListIter   end() const { return data_.end(); }
    const T& operator[]( const size_t index ) const { return data_[index]; }
    //@}

    //! Size query
    std::size_t size() const { return data_.size(); }

    
    //--------------------------------------------------------------------------
    /** Main function for reading the variables from a stream.
     *  The list is found within '{' and '}', comments and whitespaces are
     *  discarded.
     *  \param[in]  in  Stream in which a list is found
     *  \return         Reference to the input stream
     */
    std::istream& read( std::istream& in )
    {
        // ead whitespace characters
        in >> std::ws;
        // read next non-whitespace character
        char beginList;
        in.get( beginList );
        // assert this character
        VERIFY_MSG( beginList == beginSymbol_,
                    x2s("List must start with ") + x2s( beginSymbol_ ) );

        // until the endSymbol has been found
        while ( in.peek() != endSymbol_ ) {
            // skip any encountered comment line
            base::io::detail_::skip_comment( in );
            // read in value for the list
            T newValue;
            in >> newValue;
            // store new value
            data_.push_back( newValue );
            // eat white spaces
            in >> std::ws;
        }

        return in;
    }

    //--------------------------------------------------------------------------
    //! Write every value of the list to the stream, seperated by newline
    std::ostream& write( std::ostream& out ) const
    {
        const std::size_t listSize = this -> size();
        
        for ( std::size_t el = 0; el < listSize; el++ ) {
            const char endLine = ( el == listSize-1 ? ' ' : '\n' );
            
            out << data_[el] << endLine;
        }

        return out;
    }


    //--------------------------------------------------------------------------
    //! @name Friend i/o-operators to be used by PropertiesParser
    //@{
    friend
    std::istream& operator>>( std::istream& in, List<T>& list )
    {
        return list.read( in );
    }

    friend
    std::ostream& operator<<( std::ostream& out, const List<T>& list )
    {
        return list.write( out );
    }
    //@}


private:
    std::vector<T> data_;     //!< Dynamic storage of data
    const char beginSymbol_;  //!< Begin of list symbol
    const char endSymbol_;    //!< End of list symbol
};


#endif
