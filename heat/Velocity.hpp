//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Velocity.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef heat_velocity_hpp
#define heat_velocity_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>

//------------------------------------------------------------------------------
namespace heat{

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    mat::Vector velocityHistory( const GEOMELEMENT*  geomEp,
                                 const FIELDELEMENT* fieldEp,
                                 const typename
                                 FIELDELEMENT::FEFun::VecDim& xi )
    {
        // Evaluate the displacement gradient
        const typename base::VectorType<FIELDELEMENT::DegreeOfFreedom::size,
                                        double>::Type
            v = base::post::evaluateFieldHistory<HIST>( geomEp, fieldEp, xi );

        // return the entry
        mat::Vector v3 = mat::Vector::Zero();
        v3.head( FIELDELEMENT::DegreeOfFreedom::size ) = v;
        return v3;
    }
    
    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    mat::Vector velocity( const GEOMELEMENT*  geomEp,
                          const FIELDELEMENT* fieldEp,
                          const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return  velocityHistory<0>( geomEp, fieldEp, xi );
    }
    
}
#endif
