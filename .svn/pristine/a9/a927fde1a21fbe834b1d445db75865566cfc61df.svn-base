//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   DegreeOfFreedom.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_dof_degreeoffreedom_hpp
#define base_dof_degreeoffreedom_hpp

//------------------------------------------------------------------------------
// std   includes
#include <algorithm>
// boost includes
#include <boost/utility.hpp>
#include <boost/array.hpp>
// base  includes
#include <base/numbers.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace dof{
        
        template<unsigned SIZE, unsigned NHIST = 0>
        class DegreeOfFreedom;
    }
}

//------------------------------------------------------------------------------
/** Storage of a local degree of freedom which can be a vector quantity.
 *  Such a 'dof' is represented by an array of values and an array of indices.
 *  \tparam SIZE  Number of values representing the local dof
 *  \tparam NHIST Number of previous solution values to be stored
 */
template<unsigned SIZE, unsigned NHIST>
class base::dof::DegreeOfFreedom
    : public boost::noncopyable
{
public:
    //! @name Template parameter
    //@{
    static const unsigned size  = SIZE;  //!< dimension of the degree of freedom
    static const unsigned nHist = NHIST; //!< number of history terms stored
    //@}

    //! Index storage
    typedef boost::array<unsigned,size>       IndexArray;
    
    //! Value storage
    typedef boost::array<number,size>         ValueArray;

    //! Activity storage
    typedef boost::array<bool,size>           ActivityArray;

    //! History storage
    typedef boost::array<ValueArray,nHist+1>  HistoryStorage;

    //--------------------------------------------------------------------------
    /** Empty constructor invalidates the member data
     */
    DegreeOfFreedom()
    {
        indices_.assign( base::invalidInt );
        for ( unsigned h = 0; h < nHist+1; h++ )
            values_[h].assign(  0. );
        activity_.assign( true );
    }

    //--------------------------------------------------------------------------
    //! @name Index Mutators
    //@{
    template<typename INPITER>
    void setIndices( INPITER iter )
    {
        for ( unsigned d = 0; d < size; d ++ ) indices_[d] = *iter++;
    }

    void setIndex( const unsigned which, const unsigned value )
    {
        indices_[ which ] = value;
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Index Accessors
    //@{
    template<typename OUTITER>
    void getIndices( OUTITER iter ) const
    {
        std::copy( indices_.begin(), indices_.end(), iter );
    }

    unsigned getIndex( const unsigned which ) const
    {
        return indices_[ which ];
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Values
    //@{
    void setValue( const unsigned which, const number value )
    {
        values_[0][ which ] = value;
    }

    number getValue( const unsigned which ) const
    {
        return values_[0][which];
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name History
    //@{
    template<unsigned H>
    number getHistoryValue( const unsigned which ) const
    {
        STATIC_ASSERT_MSG( (H <= nHist), "Invalid history access" );
        return values_[H][which];
    }

    void pushHistory()
    {
        for ( unsigned h = nHist; h > 0; h-- )
            for ( unsigned s = 0; s < size; s++ )
                values_[h][s] = values_[h-1][s];
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Activity
    //@{
    template<typename OUTITER>
    void getActivity( OUTITER iter ) const
    {
        std::copy( activity_.begin(), activity_.end(), iter );
    }

    bool isActive( const unsigned which ) const { return activity_[ which ]; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Constraints
    //@{
    void constrainValue( const unsigned which, const number value )
    {
        activity_[which] = false;
        
        if ( constraints_.empty() )
            constraints_.resize( size, base::invalidReal() );
        
        constraints_[ which ] = value;
    }

    number getConstraint( const unsigned which ) const
    {
        return constraints_[which];
    }

    template<typename OUTITER>
    void getConstraints( OUTITER iter ) const
    {
        for ( unsigned d = 0; d < size; d ++ ) {
            if ( activity_[d] )
                *iter = base::invalidReal();
            else
                *iter = constraints_[d];
            ++iter;
        }
    }
    //@}
    
private:
    IndexArray     indices_;  //!< Storage of the dof indices
    //ValueArray    values_;   //!< Storage of the dof values
    HistoryStorage values_;
    ActivityArray  activity_; //!< Flags if DoF-entries are active

    std::vector<number> constraints_;
};

#endif
