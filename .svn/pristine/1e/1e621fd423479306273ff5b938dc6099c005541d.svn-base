//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   evaluations.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef fluid_evaluations_hpp
#define fluid_evaluations_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace fluid{

    //--------------------------------------------------------------------------
    template<unsigned HIST,typename GEOMELEMENT, typename FIELDELEMENT>
    typename base::VectorType<FIELDELEMENT::DegreeOfFreedom::size,double>::Type
    velocityHistory( const GEOMELEMENT*  geomEp,
                     const FIELDELEMENT* fieldEp,
                     const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return base::post::evaluateFieldHistory<HIST>( geomEp, fieldEp, xi );
    }

    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    typename base::VectorType<FIELDELEMENT::DegreeOfFreedom::size,double>::Type
    velocity( const GEOMELEMENT*  geomEp,
              const FIELDELEMENT* fieldEp,
              const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return base::post::evaluateField( geomEp, fieldEp, xi );
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    typename base::MatrixType<GEOMELEMENT::Node::dim,
                              FIELDELEMENT::DegreeOfFreedom::size,
                              double>::Type
    velocityGradientHistory( const GEOMELEMENT*  geomEp,
                             const FIELDELEMENT* fieldEp,
                             const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        // return the entry
        return base::post::evaluateFieldGradientHistory<HIST>( geomEp, fieldEp, xi );
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    double pressureHistory( const GEOMELEMENT*  geomEp,
                            const FIELDELEMENT* fieldEp,
                            const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        // return the entry
        return base::post::evaluateFieldHistory<HIST>( geomEp, fieldEp, xi )[0];
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    double velocityDivergenceHistory( const GEOMELEMENT*  geomEp,
                                      const FIELDELEMENT* fieldEp,
                                      const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        // Evaluate the displacement gradient
        const typename base::MatrixType<GEOMELEMENT::Node::dim,
                                        FIELDELEMENT::DegreeOfFreedom::size,
                                        double>::Type
            gradU = base::post::evaluateFieldGradientHistory<HIST>( geomEp,
                                                                    fieldEp, xi );
        
        double divU = 0.;
        for ( unsigned d = 0; d < FIELDELEMENT::DegreeOfFreedom::size; d++ )
            divU += gradU(d,d);

        return divU;
    }

    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    double velocityDivergence( const GEOMELEMENT*  geomEp,
                               const FIELDELEMENT* fieldEp,
                               const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return velocityDivergenceHistory<0>( geomEp, fieldEp, xi );
    }


}
#endif
