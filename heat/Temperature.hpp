//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Temperature.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef heat_temperature_hpp
#define heat_temperature_hpp

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
    double temperatureHistory( const GEOMELEMENT*  geomEp,
                               const FIELDELEMENT* fieldEp,
                               const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        // Evaluate the displacement gradient
        const typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                        double>::Type
            u = base::post::evaluateFieldHistory<HIST>( geomEp, fieldEp, xi );

        // return the entry
        return u[0];
    }


    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    double temperature( const GEOMELEMENT*  geomEp,
                        const FIELDELEMENT* fieldEp,
                        const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return temperatureHistory<0>( geomEp, fieldEp, xi );
    }


    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    mat::Vector temperatureGradientHistory( const GEOMELEMENT*  geomEp,
                                            const FIELDELEMENT* fieldEp,
                                            const typename
                                            FIELDELEMENT::FEFun::VecDim& xi )
    {
        mat::Vector gradU = mat::Vector::Constant( 0. );

        // Evaluate the displacement gradient
        const typename base::Matrix<GEOMELEMENT::Node::dim,
                                        FIELDELEMENT::DegreeOfFreedom::size,
                                        double>::Type
            aux = base::post::evaluateFieldGradientHistory<HIST>( geomEp,
                                                                  fieldEp, xi );
        
        // Add displacement gradient to tensor
        gradU.head( aux.size() ) += aux;;

        return gradU;
    }

    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    mat::Vector temperatureGradient( const GEOMELEMENT*  geomEp,
                                     const FIELDELEMENT* fieldEp,
                                     const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return temperatureGradientHistory<0>( geomEp, fieldEp, xi );
    }

}
#endif
