#ifndef boundingbox_h
#define boundingbox_h

#include <base/linearAlgebra.hpp>
#include <base/verify.hpp>
#include <base/auxi/compareNumbers.hpp>

//------------------------------------------------------------------------------
/** Represenation of an axis-aligned bounding box with some test functions.
 *
 *  \tparam DIM Spatial dimension of the box
 */
template<unsigned DIM>
class BoundingBox
{
public:
    static const unsigned dim = DIM;
    
    typedef typename base::Vector<dim,double>::Type VecDim;

    BoundingBox( const VecDim lower,
                 const VecDim upper )
        : lower_( lower ),
          upper_( upper )
    {
        for ( unsigned d = 0; d < dim; d++ )
            VERIFY_MSG( lower_[d] < upper_[d],
                        "upper has to be larger than lower" );
    }

    bool isInside( const VecDim& x ) const
    {
        for ( unsigned d = 0; d < dim; d++ ) {
            if ( x[d] > upper_[d] ) return false;
            if ( x[d] < lower_[d] ) return false;
        }
        return true;
    }

    VecDim snap( const VecDim& x ) const
    {
        VecDim y = x;
        for ( unsigned d = 0; d < dim; d++ ) {
            if ( x[d] > upper_[d] ) y[d] = upper_[d];
            if ( x[d] < lower_[d] ) y[d] = lower_[d];
        }
        return y;
    }

    bool isOnLowerBoundary( const VecDim&  x,
                            const unsigned dir,
                            const double   tol ) const
    {
        return base::auxi::almostEqualNumbers( x[dir], lower_[dir], tol );
    }

    bool isOnUpperBoundary( const VecDim&  x,
                            const unsigned dir,
                            const double   tol ) const
    {
        return base::auxi::almostEqualNumbers( x[dir], upper_[dir], tol );
    }

    bool isOnAnyBoundary( const VecDim& x, const double tol ) const
    {
        for ( unsigned d = 0; d < dim; d++ ) {
            if ( base::auxi::almostEqualNumbers( x[d], upper_[d], tol ) ) return true;
            if ( base::auxi::almostEqualNumbers( x[d], lower_[d], tol ) ) return true;
        }
        return false;
    }

    VecDim localCoordinate( const VecDim& x ) const
    {
        VecDim xi;
        for ( unsigned d = 0; d < dim; d++ )
            xi[d] = (x[d] - lower_[d]) / (upper_[d] - lower_[d]);

        return xi;
    }

private:
    const VecDim lower_;
    const VecDim upper_;
};


#endif
