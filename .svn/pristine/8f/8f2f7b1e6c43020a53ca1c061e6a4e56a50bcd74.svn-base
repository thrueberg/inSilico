//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   MultiIndex.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_multiindex_hpp
#define base_multiindex_hpp

//------------------------------------------------------------------------------
// std includes
#include <cassert>
// eigen3 includes
#include <Eigen/Core>
#include <Eigen/Dense>

//------------------------------------------------------------------------------
namespace base{

    template<unsigned SIZE>
    struct MultiIndex;

}


//------------------------------------------------------------------------------
template<unsigned SIZE>
struct base::MultiIndex
{
    //! Size of the multi-index
    static const unsigned size = SIZE;
    
    //--------------------------------------------------------------------------
    /** Use Eigen3's array object for a multi-index representation
     *  \tparam size Dimension of the multi-index
     */
    typedef Eigen::Array<int,size,1> Type;

    //--------------------------------------------------------------------------
    /**
     */
    static Type constant( const int value = 0 )
    {
        return Type::Constant( value );
    }

    //--------------------------------------------------------------------------
    /** Compute the 'length' \f$ l \f$ of a multi-index \f$ m \f$ of size
     *  \f$ s \f$  defined as
     *  \f[
     *         l = \prod_{i=1}^s m[i]
     *  \f]
     *  \param[in] mi Multi-index whose length is returned
     *  \returns length of the index
     */
    static std::size_t length( const Type & mi )
    {
        std::size_t result = 1;
        for ( unsigned s = 0; s < size; s++ ) {
            assert( mi[s] > 0 );
            result *= mi[s];
        }
        return result;
    }

    //--------------------------------------------------------------------------
    /** Unwrap a multi-index \f$ m \f$ which is defined with respect to the
     *  dimensions \f$ d \f$ by computing
     *  \f[
     *       \sum_{i=1}^s ( m[i] \prod_{j=0}^{s-1} d[j] )
     *  \f]
     *  using the convention \f$ d[0] = 1 \f$. This gives the linear index
     *  corresponding to a multi-index in a multi-dimensional array with the
     *  given dimensions.
     *  \param[in] mi          Index to linearise
     *  \param[in] dimensions  Dimensions of the multi-array
     *  \returns Linearised index
     */
    static std::size_t unwrap( const Type & mi, const Type & dimensions )
    {
        std::size_t result = 0;
        std::size_t stride = 1;
        for ( unsigned s = 0; s < size; s++ ) {
            //! Pre-conditions for the arguments
            assert( mi[s]         >= 0 );
            assert( dimensions[s] >  0 );
            assert( mi[s] < dimensions[s] );
                
            result += mi[s] * stride;
            stride *= dimensions[s];
        }
        return result;
    }

    //--------------------------------------------------------------------------
    /** Undo the action of 'unrwap' above. Given a linear index of some
     *  position in a multi-dimensional array, construct the corresponding
     *  multi-index.
     *  For size = 3 this would be the multi-index \f$ m \f$ corresponding
     *  to the linear index \f$ l \f$ in the array of dimensions \f$ d \f$
     *  \f[
     *      m[3] = l / (d[1] d[2]), \quad
     *      m[2] = (l - (d[1] d[2]) m[3]) / d[1], \quad 
     *      m[1] = l - (d[1] d[2] m[3] + d[1] m[2])
     *  \f]
     *  (Integer division!)
     *  \param[in] unwrapped  The linearised index
     *  \param[in] dimensions The dimensions of the multi-index
     *  \returns   Constructed multi-index
     */
    static Type wrap( const std::size_t unwrapped, const Type & dimensions )
    {
        //! first compute the strides
        typename MultiIndex<size>::Type strides;
        unsigned factor = 1;
        
        for ( unsigned s = 0; s < size; s++ ) {
            assert( dimensions[s] > 0 );
            strides[s] = factor;
            factor *= dimensions[s];
        }

        //! Re-construct the multi-index
        typename MultiIndex<size>::Type result;
        std::size_t aux = unwrapped;
        for ( int s = size-1; s >=0; s-- ) {
            result[s] = static_cast<int>( aux / strides[s] );
            aux -= (result[s] * strides[s]);
        }

        return result;
        
    }
    
};

#endif

    
