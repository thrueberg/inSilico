//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ScaledBasis.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_scaledbasis_hpp
#define base_cut_scaledbasis_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <algorithm>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename SFUNBASIS> class ScaledBasis;
    }
}

//------------------------------------------------------------------------------
/** Wrapper around a shape function basis which scales the function values.
 *  Mimicking diagonal pre-conditioning in a situation with cut-cells, every
 *  shape function is scaled by the inverse of the size of its active support.
 *  Scaling values are stored in this object and applied when function or
 *  gradient evaluations are called.
 *  \tparam SFUNBASIS Shape function basis to inherit from and scale
 */
template<typename SFUNBASIS>
class base::cut::ScaledBasis : public SFUNBASIS
{
public:
    //! Template parameter: inherited shape function basis
    typedef SFUNBASIS SFunBasis;

    //! Allocate storage for scaling factors and set all to one
    ScaledBasis()
    {
        scalars_.resize( SFunBasis::Base::numFun );
        std::fill( scalars_.begin(), scalars_.end(), 1.0 );
    }

    //! Set all scalars to one
    void reset()
    {
        std::fill( scalars_.begin(), scalars_.end(), 1.0 );
    }
   
    //! Set scalar of the basis
    void setScalar( const unsigned which, const double value )
    {
        scalars_[ which ] = value;
    }

    /** Evaluate function in physical space.
     *  Evaluate shape function basis and scale the values accordingly
     */
    template<typename GEOMELEMENT>
    void evaluate( const GEOMELEMENT* geomElemPtr,
                   const typename SFunBasis::Base::VecDim& xi,
                   typename SFunBasis::Base::FunArray &result ) const
    {
        SFunBasis::evaluate( geomElemPtr, xi, result );
        for ( unsigned s = 0; s < result.size(); s++ )
            result[s] *= scalars_[s];
    }

    /** Evaluate gradient in physical space.
     *  Evalute gradiens from basis and scale the values accordingly
     */
    template<typename GEOMELEM>
    double evaluateGradient( const GEOMELEM* geomElemPtr,
                             const typename SFunBasis::Base::VecDim& xi,
                             std::vector<typename GEOMELEM::Node::VecDim>& result ) const
    {
        const double detJ = SFunBasis::evaluateGradient( geomElemPtr, xi, result );

        for ( unsigned s = 0; s < result.size(); s++ )
            result[s]  *= scalars_[s];

        // return det J
        return detJ;
    }

private:
    //! Storage of scaling values
    std::vector<double> scalars_;
};

#endif
