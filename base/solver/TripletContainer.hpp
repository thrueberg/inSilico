//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   TripletContainer.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_solver_tripletcontainer_hpp
#define base_solver_tripletcontainer_hpp

// std includes
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <utility>
// boost includes
#include <boost/unordered_set.hpp>
#include <boost/utility.hpp>
// Eigen includes
#include <Eigen/Sparse>
#include <Eigen/Core>
// base includes
#include <base/verify.hpp>
#include <base/auxi/EqualPointers.hpp>
#include <base/auxi/parallel.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/asmb/collectFromDoFs.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace solver{
        namespace detail_{

            template<typename SCALAR,
                     typename INDEX = typename Eigen::SparseMatrix<SCALAR>::Index>
            class Triplet;

            //! Overloaded hash-value function for the unordered_set of boost
            template<typename TRIP>
            std::size_t hash_value(const TRIP& t)
            {
                // Use boost's hasher for a pair of integers
                boost::hash< std::pair<typename TRIP::Index,
                                       typename TRIP::Index> > hasher;
                return hasher( std::make_pair(t.row(), t.col() ) );
            }

        }

        class TripletContainer;
    }
}

//------------------------------------------------------------------------------
/** Add more flexibility to a triplet storage.
 *  The original Triplet provided by Eigen is rather rigid as only allows
 *  constant access functions and no comparison operation.
 *  This object copies the code from the Eigen::Triplet and adds two features:
 *  1. addTo which allows to add onto an existing value
 *  2. less-than and equality comparison for sorting algorithms
 *     (the equality is necessary for the boost::unordered containers)
 *  \tparam SCALAR Type of number
 *  \tparam INDEX  Type of index
 */
template<typename SCALAR, typename INDEX>
class base::solver::detail_::Triplet
{
public:
    //! @name Template parameter
    //@{
    typedef SCALAR  Scalar;
    typedef INDEX   Index;
    //@}

    //! Default constructor
    Triplet() : row_( 0 ), col_( 0 ), value_( 0 ) {}
    
    //! Constructor with data
    Triplet(const Index& i, const Index& j, const Scalar& v = Scalar(0) )
        : row_( i ), col_( j ), value_( v )
    {}

    //! @name Private data accessor
    //@{
    const Index& row() const { return row_; }
    const Index& col() const { return col_; }
    const Scalar& value() const { return value_; }
    //@}

    //! Modifier: add to the stored value
    void addTo( const Scalar& value )
    {
#pragma omp atomic
        value_ += value;
    }

    //! @name Comparison operators
    //@{

    //! Required by std::set
    bool operator<( const Triplet<Scalar,Index>& other ) const
    {
        return
            ( row_  < other.row() ) or
            ( not( other.row() < row_ ) and
              ( col_ < other.col() ) );
    }

    //! Required by boost::unordered_set
    bool operator==( const Triplet<Scalar,Index>& other ) const
    {
        return ( row_ == other.row() ) and ( col_ == other.col() );
    }
    //@}
        
private:
    Index  row_, col_; //!< Row and column indices
    Scalar value_;     //!< Matrix entry
};

//------------------------------------------------------------------------------
/** Container for the temporary matrix storage.
 *  This object has two fundamental modes
 *  1.  Pre-structured (storing a non-zero matrix pattern)
 *  2.  Dynamic
 *  The first case is useful to use parallel assembly since the addition of
 *  values to the pre-defined non-zero structure is thread-safe. For simplicity
 *  and downward compatibility, there is dynamic version in which the triplet
 *  storage grows on demand.
 */
class base::solver::TripletContainer : boost::noncopyable
{
public:
    // Type of underlying storage for (i,j,value) tuples
    typedef base::solver::detail_::Triplet<base::number> Triplet;

    //! Constructure
    TripletContainer() : preStructured_( false )
    { }

    //--------------------------------------------------------------------------
    /** Registering of active DoF IDs for test and trial fields in order to
     *  pre-determine the non-zero pattern of the system matrix.
     *  Going through all the tuples of a given field-binder, this object
     *  creates storage for all requested triplets. This storage is sorted
     *  and copied into a vector.
     *  \tparam FIELDTUPLEBINDER Type to determine which is test and trial field
     *  \tparam FIELDBINDER      Type of mesh and field binder
     *  \param[in] fieldBinder Binder of mesh and fields
     */
    template<typename FIELDTUPLEBINDER, typename FIELDBINDER>
    void registerFields( const FIELDBINDER& fieldBinder )
    {
        // Convenience typedefs
        typedef typename FIELDTUPLEBINDER::Tuple Tuple;
        typedef typename Tuple::TestElement      TestElement;
        typedef typename Tuple::TrialElement     TrialElement;

        // flag for symmetry
        bool isBubnov = false;
        
        // Go through all elements
        typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
        typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
        for ( std::size_t i = 0; iter != end; ++iter, i++ ) {
            
            // make the tuple of field element pointers
            Tuple tuple = FIELDTUPLEBINDER::makeTuple( *iter );

            // extract test and trial elements from tuple
            TestElement*  testEp  = tuple.testElementPtr();
            TrialElement* trialEp = tuple.trialElementPtr();

            // if pointers are identical, Galerkin-Bubnov scheme (check only once)
            if ( i == 0 ) isBubnov =
                              base::auxi::EqualPointers<TestElement,
                                                        TrialElement>::apply( testEp,
                                                                              trialEp );

            // dof statuses 
            std::vector<base::dof::DoFStatus> rowDoFStatus, colDoFStatus;

            // dof IDs 
            std::vector<std::size_t> rowDoFIDs, colDoFIDs;

            // dof values (just placeholder)
            std::vector<base::number> rowDoFValues, colDoFValues;

            // dof constraints 
            typedef std::pair<unsigned, std::vector< std::pair<base::number,std::size_t> > >
                WeightedDoFIDs;
            std::vector<WeightedDoFIDs> rowConstraints, colConstraints;

            // Collect dof entities from element
            bool doSomething = 
                base::asmb::collectFromDoFs( testEp, rowDoFStatus,
                                             rowDoFIDs, rowDoFValues,
                                             rowConstraints,
                                             false );

            // if no row dof is ACTIVE or CONSTRAINED, go to next element
            if ( not doSomething ) continue;
            
            // In case of identical test and trial spaces, just copy
            if ( isBubnov ) {
                colDoFStatus   = rowDoFStatus;
                colDoFIDs      = rowDoFIDs;
                colDoFValues   = rowDoFValues;
                colConstraints = rowConstraints;
            }
            else // otherwise, collect for trial space
                doSomething = 
                    base::asmb::collectFromDoFs( trialEp, colDoFStatus,
                                                 colDoFIDs, colDoFValues,
                                                 colConstraints,
                                                 false );

            // if no col dof is ACTIVE or CONSTRAINED, go to next element
            if ( not doSomething ) continue;


            // Collect all ID numbers (ACTIVE and CONSTRAINED)
            std::vector<std::size_t> effRowDoFIDs, effColDoFIDs;
            {
                // collect active row dof IDs from this element
                for ( unsigned r = 0; r < rowDoFIDs.size(); r++ )
                    if ( rowDoFStatus[r] == base::dof::ACTIVE )
                        effRowDoFIDs.push_back( rowDoFIDs[r] );

                // collect master dof IDs from linear constraints
                for ( unsigned r = 0; r < rowConstraints.size(); r++ ) {
                    for ( unsigned r2 = 0; r2 < (rowConstraints[r]).second.size(); r2++ ) {
                        effRowDoFIDs.push_back( (rowConstraints[r].second[r2].second) );
                    }
                }

                // if test- and trial-spaces are equal, just copy the IDs
                if ( isBubnov ) effColDoFIDs = effRowDoFIDs;
                // otherwise, do the same for the column (i.e trial) space
                else {
                    for ( unsigned c = 0; c < colDoFIDs.size(); c++ )
                        if ( colDoFStatus[c] == base::dof::ACTIVE )
                            effColDoFIDs.push_back( colDoFIDs[c] );

                    for ( unsigned c = 0; c < colConstraints.size(); c++ ) {
                        for ( unsigned c2 = 0; c2 < (colConstraints[c]).second.size(); c2++ ) {
                            effColDoFIDs.push_back( (colConstraints[c].second[c2].second) );
                        }
                    }
                }
            }

            // go through all indices
            for ( std::size_t r = 0; r < effRowDoFIDs.size(); r++ ) {
                for ( std::size_t c = 0; c < effColDoFIDs.size(); c++ ) {

                    Triplet triplet( static_cast<typename Triplet::Index>(effRowDoFIDs[r]),
                                     static_cast<typename Triplet::Index>(effColDoFIDs[c]), 0. );
                    // insertion only works if (r,c) has not been insertet yet
                    tmpTriplets_.insert( triplet );
                } // for all columns
            } // for all rows

            
        } // for all element tuples

        // append to vector of triplets
        {
            // size of existing storage
            const std::size_t currentNum = triplets_.size();
            // resize the storage
            triplets_.resize( currentNum +
                              std::distance( tmpTriplets_.begin(), tmpTriplets_.end() ) );

            // iterator to the end of the existing storage
            typename std::vector<Triplet>::iterator insertIter = triplets_.begin();
            std::advance( insertIter, currentNum );

            // copy from temporary storage into new storage
            std::copy( tmpTriplets_.begin(), tmpTriplets_.end(), insertIter );

            // destroy current storage
            tmpTriplets_.clear();

            // if container was not empty, it needs to be sorted
            if ( currentNum > 0 ) {
                std::sort( triplets_.begin(), triplets_.end() );
            }
        }

        // set flag that pre-structured version is used
        preStructured_ = true;
        return;
    }

    //--------------------------------------------------------------------------
    // Insert a value at given row and column position
    void insert( const unsigned i, const unsigned j, const base::number value )
    {
        // new triplet
        Triplet triplet( i, j, value );

        // For non-pre-structured case, use dynamic memory
        if ( not preStructured_ ){

            // important check
#ifdef _OPENMP
#if NTHREADS > 1
            VERIFY_MSG( false, "Multiple threads are not allowed for this method" );
#endif
#endif

            // try inserting into triplet storage (only if new)
            const std::pair<std::set<Triplet>::iterator,bool>
                check = tmpTriplets_.insert( triplet );

            // in case insertion values, replace triplet by sum of entries
            if ( not check.second ) {
                triplet.addTo( check.first -> value() );
                
                std::set<Triplet>::iterator hint = check.first;
                hint++;
                tmpTriplets_.erase( check.first );
                tmpTriplets_.insert( hint, triplet );
            }
            
        }
        else{
            // use pre-determined non-zero pattern

            // find the existing (i,j) pair (pre-condition: sorted array)
            std::vector<Triplet>::iterator iter =
                std::lower_bound( triplets_.begin(), triplets_.end(), triplet );

            // make sure that this is the right one (lower_bound!!)
            const bool check =
                (iter -> row() == static_cast<Triplet::Index>(i) ) and
                (iter -> col() == static_cast<Triplet::Index>(j) );
            VERIFY_MSG( check, "TripletContainer had not been properly set up" );

            // add to entry
            iter -> addTo( value );
        }

    }

    //--------------------------------------------------------------------------
    //! In case of no pre-structure, create vector from set
    void prepare()
    {
        if ( not preStructured_ ) {

            // Eigen requires random-access iterators
            triplets_.resize( std::distance( tmpTriplets_.begin(), tmpTriplets_.end() ) );
            std::copy( tmpTriplets_.begin(), tmpTriplets_.end(), triplets_.begin() );

            // destroy current storage
            tmpTriplets_.clear();
        }
        return;
    }

    //--------------------------------------------------------------------------
    //! Begin iterator to Triplet storage
    std::vector<Triplet>::const_iterator begin() const
    {
        return triplets_.begin();
    }

    //! End iterator to Triplet storage
    std::vector<Triplet>::const_iterator end() const
    {
        return triplets_.end();
    }

    //--------------------------------------------------------------------------
    //! Clear local storage
    void destroy()
    {
        triplets_.clear();
        std::vector<Triplet>().swap( triplets_ );
    }

    //! Debug routine: print triplet
    std::ostream& write( std::ostream& out ) const
    {
        for ( std::size_t t = 0; t < triplets_.size(); t++ ) {
            out << triplets_[t].row() << " "
                << triplets_[t].col() << " "
                << triplets_[t].value() << std::endl;
        }

        return out;
    }
    
private:
    bool                 preStructured_; //!< Mode of storage
    std::set<Triplet>    tmpTriplets_;   //!< Temporary dynamic storage
    std::vector<Triplet> triplets_;      //!< Final triplet storage (sorted)
};



#endif
