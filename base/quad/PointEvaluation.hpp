//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   PointEvaluation.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_quad_pointevaluation_hpp
#define base_quad_pointevaluation_hpp

//------------------------------------------------------------------------------
// std   includes
#include <utility>
// boost includes
#include <boost/array.hpp>
#include <boost/utility.hpp>
// base  includes
#include <base/verify.hpp>
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace quad{
        
        class PointEvaluation;

    }
}

//------------------------------------------------------------------------------
/** \brief Point evaluation as zero-dimensional quadrature.
 *  \tparam DEGREE Dummy degree
 */
class base::quad::PointEvaluation
    : public boost::noncopyable
{
public:
    //! The local coordinate dimension
    static const unsigned dim = 0;

    //! Number of points
    static const unsigned numPoints = 1;
    
    //! Type of coordinate-vector
    typedef base::Vector<dim>::Type VecDim;
    
    //! Type of iterator for external access
    typedef boost::array<std::pair<double, VecDim>, numPoints>::const_iterator
    Iter;

    PointEvaluation()
    {
        VecDim dummy;
        weightsAndPoints_[0] = std::make_pair( 1.0, dummy );
    }

    //! Begin of array iterator
    Iter begin() const  { return weightsAndPoints_.begin(); }
    //! End of array iterator
    Iter end()   const  { return weightsAndPoints_.end();   }

private:
    //! Pair of weights and points representing the NPTS-point quadrature rule
    boost::array<std::pair<double,VecDim>, numPoints> weightsAndPoints_;
};

#endif
