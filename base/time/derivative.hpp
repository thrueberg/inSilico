//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   derivative.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_time_derivative_hpp
#define base_time_derivative_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
// base includes
#include <base/linearAlgebra.hpp>
// base/time includes
#include <base/time/BDF.hpp>
// base/post includes
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace time{

        namespace detail_{

            //! Evaluate the field variable in history times some weight
            template<typename GEOMELEMENT, typename FIELDELEMENT, int CTR>
            struct WeightedFieldValue
            {
                static void apply( const GEOMELEMENT*  geomElemPtr,
                                   const FIELDELEMENT* fieldElemPtr,
                                   const typename FIELDELEMENT::FEFun::VecDim& xi,
                                   const std::vector<double> weights,
                                   typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                                         base::number>::Type& result )
                {
                    // which history to access
                    static const unsigned past =
                        FIELDELEMENT::DegreeOfFreedom::nHist - CTR;// + 1;

                    
                    // number of weights
                    const unsigned numWeights = static_cast<unsigned>( weights.size() );
                    // quit if zero
                    if ( numWeights == 0 ) return;

                    // call history evaluation
                    typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                          base::number>::Type tmp
                        = base::post::evaluateFieldHistory<past>( geomElemPtr, fieldElemPtr, xi );
                    
                    // add weighted history to result
                    result += weights[0] * tmp;

                    // quit if this was the last weight
                    if ( numWeights == 1 ) return;

                    // truncate weights by discarding the first one
                    std::vector<double>::const_iterator first = weights.begin();
                    first++;
                    const std::vector<double> remainingWeights =
                        std::vector<double>( first, weights.end() );

                    // recursive call with remaning weights
                    WeightedFieldValue<GEOMELEMENT,FIELDELEMENT,CTR-1>::apply( geomElemPtr,
                                                                               fieldElemPtr,
                                                                               xi,
                                                                               remainingWeights,
                                                                               result );
                    return;
                }
            };

            //! Recursion stop
            template<typename GEOMELEMENT, typename FIELDELEMENT>
            struct WeightedFieldValue<GEOMELEMENT,FIELDELEMENT,-1>
            {
                static void apply( const GEOMELEMENT*  geomElemPtr,
                                   const FIELDELEMENT* fieldElemPtr,
                                   const typename FIELDELEMENT::FEFun::VecDim& xi,
                                   const std::vector<double> weights,
                                   typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                                         base::number>::Type& result )
                {
                    return; // Stop recursion
                }
            };
            
        }
        
        //----------------------------------------------------------------------
        /** Evaluate the first time derivative of a field.
         *  In order to evaluate the time derivative of a field variable, a
         *  BDF approximation is used
         *  \f[
         *      \dot{u}^{n+1} \approx \sum_{s=0}^A a_s u^{n+1-s}
         *  \f]
         *  where the method order \f$ Q = A-1 \f$ is determined from the number
         *  of history terms stored in the degrees of freedom.
         *  The weighted sum in this formula is achieved by delegation to a
         *  recursive template. The startup phase for \f$ n+1 < Q \f$ is taken
         *  into account via the base::time::MultiStep interface.
         *  \tparam GEOMELEMENT  Type of geometry element
         *  \tparam FIELDELEMENT Type of field element
         *  \param[in] geomElemPtr  Pointer to geometry element
         *  \param[in] fieldElemPtr Pointer to field element
         *  \param[in] xi           Local evaluation coordinate
         *  \param[in] stepSize     Size of the time step
         *  \param[in] step         Number of step (needed for the startup)
         *  \return                 Approximate time derivative of the field
         */
        template<typename GEOMELEMENT, typename FIELDELEMENT>
        typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                  base::number>::Type
        evaluateTimeDerivative( const GEOMELEMENT*  geomElemPtr,
                                const FIELDELEMENT* fieldElemPtr,
                                const typename FIELDELEMENT::FEFun::VecDim& xi,
                                const double stepSize, 
                                const unsigned step )
        {
            // evaluate with maximal available order, determined by storage size
            static const unsigned evalOrder =
                FIELDELEMENT::DegreeOfFreedom::nHist;

            // use BDF for approximation of the time derivative 
            std::vector<double> bdfWeights;
            base::time::BDF<evalOrder>::derivativeWeights( step, bdfWeights );
            
            // divide by step size
            for ( unsigned s = 0; s < bdfWeights.size(); s++ )
                bdfWeights[s] /= stepSize;

            // initialise result with zero
            typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                  base::number>::Type result =
                base::constantVector<FIELDELEMENT::DegreeOfFreedom::size>( 0. );


            // make recursive call
            detail_::WeightedFieldValue<GEOMELEMENT,FIELDELEMENT,
                                        evalOrder>::apply( geomElemPtr,
                                                           fieldElemPtr,
                                                           xi, bdfWeights,
                                                           result );

            return result;
        }
        
    }
}
#endif
