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
#include <vector>
// boost includes
#include <boost/utility.hpp>
#include <boost/array.hpp>
// base  includes
#include <base/numbers.hpp>
#include <base/dof/Constraint.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace dof{
        
        template<unsigned SIZE, unsigned NHIST = 0>
        class DegreeOfFreedom;

        //! Give the possible status of a DoF a name
        enum DoFStatus
        {
            ACTIVE,
            CONSTRAINED,
            INACTIVE
        };

        //! Helper for output of a DoF status
        template<typename DOF>
        struct Status
        {
            static double apply( const DOF* doFPtr, const unsigned component )
            {
                if ( doFPtr -> isConstrained( component ) ) return 0.;
                if ( doFPtr -> isActive(      component ) ) return 1.;
                return -1.;
            }
        };

       
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
    typedef boost::array<std::size_t,size>    IndexArray;
    
    //! Value storage
    typedef boost::array<number,size>         ValueArray;

    //! Activity storage
    typedef boost::array<DoFStatus,size>      StatusArray;

    //! History storage
    typedef boost::array<ValueArray,nHist+1>  HistoryStorage;

    //! Self type
    typedef base::dof::DegreeOfFreedom<size,nHist> SelfType;

    //! Type of constraint
    typedef base::dof::Constraint<SelfType>       Constraint;

    //--------------------------------------------------------------------------
    /** Empty constructor invalidates the member data
     */
    DegreeOfFreedom()
        : id_( base::invalidInt )
    {
        indices_.assign( base::invalidInt );
        for ( unsigned h = 0; h < nHist+1; h++ )
            for ( unsigned s = 0; s < size; s++ )
                values_[h][s] = 0.;
        
        status_.assign( ACTIVE );

        constraints_.assign( NULL );
    }

    ~DegreeOfFreedom()
    {
        for ( unsigned s = 0; s < constraints_.size(); s++ )
            delete constraints_[s];
    }

    //--------------------------------------------------------------------------
    //! @name ID methods
    //@{
    void setID( const std::size_t id ) { id_ = id; }
    std::size_t getID() const { return id_; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Index Mutators
    //@{
    template<typename INPITER>
    void setIndices( INPITER iter )
    {
        for ( unsigned d = 0; d < size; d ++ ) indices_[d] = *iter++;
    }

    void setIndex( const unsigned which, const std::size_t value )
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

    std::size_t getIndex( const unsigned which ) const
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

    void clearValue( )
    {
        values_[0].assign( 0. );
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

    template<unsigned H>
    void setHistoryValue( const unsigned which, const number value )
    {
        values_[H][which] = value;
    }

    void pushHistory()
    {
        for ( unsigned h = nHist; h > 0; h-- )
            for ( unsigned s = 0; s < size; s++ )
                values_[h][s] = values_[h-1][s];
    }
    //@}

    //--------------------------------------------------------------------------
    //! Scale all stored values uniformly by a factor
    void scaleAllValues( const double factor )
    {
        for ( unsigned h = 0; h <= nHist; h++ )
            for ( unsigned s = 0; s < size; s++ )
                values_[h][s] *= factor;
        return;
    }

    //--------------------------------------------------------------------------
    //! @name Status
    //@{
    template<typename OUTITER>
    void getStatus( OUTITER iter ) const
    {
        std::copy( status_.begin(), status_.end(), iter );
    }

    bool isActive( const unsigned which ) const
    {
        return ( status_[ which ] == ACTIVE );
    }

    bool isConstrained( const unsigned which ) const
    {
        return ( status_[ which ] == CONSTRAINED );
    }

    void activateAll()   { status_.assign( ACTIVE );   }
    void deactivateAll() { status_.assign( INACTIVE ); }
    //@}

    //--------------------------------------------------------------------------
    //! @name Constraints
    //@{
    void constrainValue( const unsigned which, const number value )
    {
        // if dof had been active, inactivate and create a constraint
        if ( status_[which] == ACTIVE ) {
            status_[which] = CONSTRAINED;
            constraints_[ which ] = new Constraint( value );
        }
        // otherwise pass the value to the constraint
        else if ( status_[which] == CONSTRAINED )
            constraints_[ which ] -> setValue( value );
    }

    //--------------------------------------------------------------------------
    /** Get the values the constrained DoF components are prescribed.
     *  In case of no constraint give an invalid number. For the case of an
     *  incremental analysis, return the difference between the prescribed
     *  value and the current value of the DoF components.
     *  \tparam OUTITER Type of iterator to copy into.
     *  \param[in,out]  iter        Iterator to some storage
     *  \param[in]      incremental True for an incremental analysis
     */
    template<typename OUTITER>
    void getPrescribedValues( OUTITER iter,
                              const bool incremental ) const
    {
        for ( unsigned d = 0; d < size; d ++ ) {
            if ( status_[d] == CONSTRAINED ) {
                if ( incremental )
                    *iter = (constraints_[d] -> getValue()) - values_[0][d];
                else
                    *iter = (constraints_[d] -> getValue());
            }
            else
                *iter = base::invalidReal();

            ++iter;
        }
    }

    //! Destroy the constraints
    void clearConstraints()
    {
        for ( unsigned s = 0; s < constraints_.size(); s++ ) {
            delete constraints_[s];
            constraints_[s] = NULL;
        }
        
        status_.assign( ACTIVE );
    }

    //! Generate a new constraint
    void makeConstraint( const unsigned which )
    {
        if ( not (status_[which] == CONSTRAINED) ) {
            status_[      which ] = CONSTRAINED;
            constraints_[ which ] = new Constraint();
        }
    }

    //! Return specific constraint object
    Constraint* getConstraint( const unsigned which )
    {
        return constraints_[which];
    }

    //! Multiply constraint RHS by a scalar
    void scaleConstraint( const double factor )
    {
        for ( unsigned d = 0; d < size; d++ ) {
            if ( status_[d] == CONSTRAINED ) {
                const number old = constraints_[d] -> getValue();
                constraints_[d] -> setValue( old * factor );
            }
        }
    }
    
    //@}

    
private:
    std::size_t    id_;       //!< ID for convenience
    
    IndexArray     indices_;  //!< Storage of the dof indices
    HistoryStorage values_;   //!< Result values and history storage
    StatusArray    status_;   //!< Flags if DoF-entries are active

    //! Array of pointers to possible constraints
    boost::array<Constraint*,size> constraints_;
};

#endif
