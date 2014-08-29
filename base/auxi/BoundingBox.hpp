//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   BoundingBox.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_auxi_boundingbox_h
#define base_auxi_boundingbox_h

#include <base/linearAlgebra.hpp>
#include <base/verify.hpp>
#include <base/auxi/compareNumbers.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace auxi{

        template<unsigned DIM> class BoundingBox;
    }
}

//------------------------------------------------------------------------------
/** Represenation of an axis-aligned bounding box.
 *
 *  \tparam DIM Spatial dimension of the box
 */
template<unsigned DIM>
class base::auxi::BoundingBox
{
public:
    //! Template parameter: spatial dimension of the BBox
    static const unsigned dim = DIM;

    //! Convenience typedef for coordinates
    typedef typename base::Vector<dim,double>::Type VecDim;

    //! Construct from given lower and upper points
    BoundingBox( const VecDim lower,
                 const VecDim upper )
        : lower_( lower ),
          upper_( upper )
    {
        for ( unsigned d = 0; d < dim; d++ )
            VERIFY_MSG( lower_[d] < upper_[d],
                        "upper has to be larger than lower" );
    }

    //! Predicate if point is inside the BBox's closure
    bool isInside( const VecDim& x ) const
    {
        for ( unsigned d = 0; d < dim; d++ ) {
            if ( x[d] > upper_[d] ) return false;
            if ( x[d] < lower_[d] ) return false;
        }
        return true;
    }

    //! Project point BBox, i.e. if outside to the BBox surface
    VecDim snap( const VecDim& x ) const
    {
        VecDim y = x;
        for ( unsigned d = 0; d < dim; d++ ) {
            if ( x[d] > upper_[d] ) y[d] = upper_[d];
            if ( x[d] < lower_[d] ) y[d] = lower_[d];
        }
        return y;
    }

    //! Check if point is on the lower boundary in specified direction
    bool isOnLowerBoundary( const VecDim&  x,
                            const unsigned dir,
                            const double   tol ) const
    {
        return base::auxi::almostEqualNumbers( x[dir], lower_[dir], tol );
    }

    //! Check if point is on the upper boundary in specified direction
    bool isOnUpperBoundary( const VecDim&  x,
                            const unsigned dir,
                            const double   tol ) const
    {
        return base::auxi::almostEqualNumbers( x[dir], upper_[dir], tol );
    }

    //! Check with tolerance if given point is on the BBox boundary
    bool isOnAnyBoundary( const VecDim& x, const double tol ) const
    {
        for ( unsigned d = 0; d < dim; d++ ) {
            if ( base::auxi::almostEqualNumbers( x[d], upper_[d], tol ) ) return true;
            if ( base::auxi::almostEqualNumbers( x[d], lower_[d], tol ) ) return true;
        }
        return false;
    }

    //! Map given point x to BBox parameter space from [0,1]
    VecDim localCoordinate( const VecDim& x ) const
    {
        VecDim xi;
        for ( unsigned d = 0; d < dim; d++ )
            xi[d] = (x[d] - lower_[d]) / (upper_[d] - lower_[d]);

        return xi;
    }

    //! Return the centre of the box
    VecDim centre() const { return 0.5 * (lower_ + upper_); }

    //! Return centre point of a surface
    VecDim surfaceCentre( const unsigned dir, const bool lower ) const
    {
        VecDim cs = this -> centre();
        cs[ dir ] = (lower ? lower_[dir] : upper_[dir] );
        return cs;
    }

    //! Normal vector to surface
    VecDim surfaceNormal( const unsigned dir, const bool lower ) const
    {
        VecDim n = base::constantVector<dim>( 0. );
        n[ dir ] = (lower ? -1.0 : 1.0 );
        return n;
    }

private:
    const VecDim lower_; //!< Coordinates of lower point
    const VecDim upper_; //!< Coordinates of upper point
};


#endif
