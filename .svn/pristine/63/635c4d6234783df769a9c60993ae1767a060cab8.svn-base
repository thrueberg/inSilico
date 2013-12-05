//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   OStreamIterator.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_io_ostreamiterator_hpp
#define base_io_ostreamiterator_hpp

//------------------------------------------------------------------------------
// std includes
#include <ostream>
#include <string>
#include <iterator>
// boost includes
#include <boost/iterator/iterator_facade.hpp>
#include <boost/function.hpp>
// base includes
#include <base/types.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{

        template<typename DATUM, typename WRITER>
        class OStreamIterator;

        namespace detail_{

            template<typename DATUM>
            struct WriteFunctorType
            {
                typedef boost::function<void( const DATUM&,
                                              std::ostream & )> Type;
            };
            
        }
       

    }
}

//------------------------------------------------------------------------------
/** Custom ostream iterator.
 *  Whereas the class ostream_iterator provided by STL inherently relies on
 *  the operator<< defined for the datum to be written, this iterator uses
 *  a custom functor provided by the caller. Therefore, the output can be
 *  formatted in that functor even for objects of external libraries which
 *  already come with an operator<< (e.g boost's and Eigen3's vector types).
 *  Also, there can only be one such operator<< in the global namespace which
 *  is quite restrictive. This custom iterator makes use of boost's
 *  iterator_facade in order to reduce the implemtation to a minimum of
 *  functions.
 *  \tparam DATUM   Type of datum to be written to an ostream
 *  \tparam WRITER  Type of functor which does the actual (formatted) writing
 *                  (defaults, for convenience, to a boost::function)
 */
template<typename DATUM,
         typename WRITER =
         typename base::io::detail_::WriteFunctorType<DATUM>::Type>
class base::io::OStreamIterator
    : public boost::iterator_facade< OStreamIterator<DATUM,WRITER>,
                                     DATUM, 
                                     std::output_iterator_tag,
                                     OStreamIterator<DATUM,WRITER> &>
{
    //! Grant access 'from above'
    friend class boost::iterator_core_access;

public:

    //! @name Template parameter
    //@{
    typedef DATUM  Datum;
    typedef WRITER Writer;
    //@}

    //! @name Convenience typedefs
    //@{
    typedef OStreamIterator<Datum,Writer>   SelfType;
    typedef std::basic_ostream<char>        OStreamType;
    //@}

private:
    //! @name Sanity checks (assure the right signature of the functor)
    //@{
    static base::TypeEquality<typename Writer::first_argument_type,
                              Datum>       sanity1;
    static base::TypeEquality<typename Writer::second_argument_type,
                              OStreamType> sanity2;
    //@}


public:
    //! Constructor with stream and Writing functor
    OStreamIterator( OStreamType& os, const Writer& writer )
        : oStreamPtr_( &os ),
          writer_(     writer )
    {}

    //! Copy constructor
    OStreamIterator( const OStreamIterator & rhs )
    : oStreamPtr_( rhs.oStreamPtr_ ),
      writer_(     rhs.writer_ )
    {}

    //--------------------------------------------------------------------------
    /** Assignment operator does the output.
     *  Using an object of this type, e.g., iter, there will be the statement
     *  '*iter = value' at some point. As the dereferencing does not do anything
     *  this can be read as 'iter = value' and the following operator is
     *  called which itself calls the internal write function.
     *  \param[in] d Reference to datum to be written
     */
    SelfType & operator=( const Datum & d )
    {
        this -> write( d );
        return *this;
    }

private:
    //! @name Implemenation required by boost's iterator facade
    //@{
    SelfType& dereference() const
    {
        return const_cast<SelfType&>( *this );
    }
    
    bool equal( const SelfType& rhs ) const
    {
        return ( oStreamPtr_ == rhs.oStreamPtr_ ); 
    }
    
    void increment(){} 
    //@}    
    
private:

    //--------------------------------------------------------------------------
    /** Main action of the iterator.
     *  Use member functor #writer_ in order to write the given datum to
     *  the stream.
     *  \param[in] d  Reference to the datum to be written
     */
    void write( const Datum& d )
    {
        if( oStreamPtr_ != NULL ) {
        
            //! use functor here
            writer_( d,  *oStreamPtr_ );

            //! sanity check
            if ( not oStreamPtr_->good() )
                oStreamPtr_ = NULL;
        }
    }

private:
    OStreamType*      oStreamPtr_; //!< pointer to stream
    const Writer&     writer_;     //!< Functor for formatted writing
};

#endif 
