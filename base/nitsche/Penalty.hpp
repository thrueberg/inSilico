//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Penalty.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_nitsche_penalty_hpp
#define base_nitsche_penalty_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/mesh includes
#include <base/mesh/Size.hpp>
// base/asmb includes
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
// base/auxi includes
#include <base/auxi/FunEvaluationPolicy.hpp>
#include <base/auxi/EqualPointers.hpp>
// base/post includes
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace nitsche{

        template<typename SURFFIELDTUPLE>
        class Penalty;

        //----------------------------------------------------------------------
        /** Convenience function for the LHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD,
                 typename PARAMETER>
        void penaltyLHS( const SURFACEQUADRATURE& surfaceQuadrature,
                         SOLVER&                  solver, 
                         const BOUNDFIELD&        boundField,
                         const PARAMETER&         parameter,
                         const double             multiplier )
        {

            // object to compute the LHS penalty term
            typedef base::nitsche::Penalty<typename SURFACETUPLEBINDER::Tuple>
                Penalty;
            Penalty penalty( multiplier );

            // integrator and assembler object
            typedef base::asmb::StiffnessMatrix<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::Tuple> SysMat;
            
            typename SysMat::Kernel kernel =
                boost::bind( &Penalty::tangentStiffness, &penalty, _1, _2, _3, _4 );

            SysMat systemMat( kernel, surfaceQuadrature, solver, false );

            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {

                // get local penalty factor
                const double factor = multiplier * parameter.penaltyWeight( iter );

                penalty.setFactor( factor );
                
                systemMat( SURFACETUPLEBINDER::makeTuple( *iter ) );
            }
            
            return;
        }

        //----------------------------------------------------------------------
        namespace detail_{
            template<typename SURFACETUPLEBINDER,
                     typename EVALUATIONPOLICY,
                     typename SURFACEQUADRATURE,
                     typename SOLVER,
                     typename BOUNDFIELD,    
                     typename PARAMETER>
            void computePenaltyRHS( const SURFACEQUADRATURE& surfaceQuadrature,
                                    SOLVER& solver, BOUNDFIELD& boundField,
                                    const typename EVALUATIONPOLICY::Fun& bcFun,
                                    const PARAMETER& parameter,
                                    const double multiplier )
            {
                // object to compute the LHS penalty term
                typedef base::nitsche::Penalty<typename SURFACETUPLEBINDER::Tuple> Penalty;
                Penalty penalty( multiplier );

                // integrator and assembler object
                typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                                    typename SURFACETUPLEBINDER::Tuple>
                    SurfaceForceInt;
                
                typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                    boost::bind( &Penalty::template residualBoundary<EVALUATIONPOLICY>,
                                 &penalty, _1, _2, _3,
                                 boost::ref( bcFun ), _4 );
                SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                                  surfaceQuadrature, solver );
                
                // apply
                typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
                typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
                for ( ; iter != end; ++iter ) {

                    const double factor = multiplier * parameter.penaltyWeight( iter );
                    
                    penalty.setFactor( factor );
                    
                    surfaceForceInt( SURFACETUPLEBINDER::makeTuple( *iter ) );
                }
            
                return;

            }
        }
        //----------------------------------------------------------------------

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        //! Version with f(x)
        template<typename SURFACETUPLEBINDER,
                 typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD,        typename BCFUN,
                 typename PARAMETER>
        void penaltyRHS( const SURFACEQUADRATURE& surfaceQuadrature,
                         SOLVER&                  solver, 
                         const BOUNDFIELD&        boundField,
                         const BCFUN&             bcFun,
                         const PARAMETER&         parameter,
                         const double             multiplier )
        {
            typedef typename SURFACETUPLEBINDER::Tuple::GeomElement::DomainElement
                DomainElement;
            typedef base::auxi::EvaluateDirectly<DomainElement,BCFUN>
                EvaluationPolicy;

            detail_::computePenaltyRHS<SURFACETUPLEBINDER,
                                       EvaluationPolicy>( surfaceQuadrature,
                                                          solver,
                                                          boundField, bcFun,
                                                          parameter, multiplier );
            return;
        }

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        //! Version with f( Ep*, xi )
        template<typename SURFACETUPLEBINDER,
                 typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD,        typename BCFUN,
                 typename PARAMETER>
        void penaltyRHSViaElement( const SURFACEQUADRATURE& surfaceQuadrature,
                                   SOLVER&                  solver, 
                                   const BOUNDFIELD&        boundField,
                                   const BCFUN&             bcFun,
                                   const PARAMETER&         parameter,
                                   const double             multiplier )
        {
            typedef typename SURFACETUPLEBINDER::Tuple::GeomElement::DomainElement
                DomainElement;
            typedef base::auxi::EvaluateViaElement<DomainElement,BCFUN>
                EvaluationPolicy;

            detail_::computePenaltyRHS<SURFACETUPLEBINDER,
                                       EvaluationPolicy>( surfaceQuadrature,
                                                          solver,
                                                          boundField, bcFun,
                                                          parameter, multiplier );
            return;
        }

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD,        
                 typename PARAMETER>
        void penaltyRHSInterface( const SURFACEQUADRATURE& surfaceQuadrature,
                                  SOLVER&                  solver, 
                                  const BOUNDFIELD&        boundField,
                                  const PARAMETER&         parameter,
                                  const double             multiplier )
        {
            // object to compute the LHS penalty term
            typedef base::nitsche::Penalty<typename SURFACETUPLEBINDER::Tuple> Penalty;
            Penalty penalty( multiplier );

            // integrator and assembler object
            typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::Tuple>
                SurfaceForceInt;
                
            typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                boost::bind( &Penalty::residualInterface, &penalty, _1, _2, _3, _4 );
            SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                             surfaceQuadrature, solver );
                
            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {

                const double factor = multiplier * parameter.penaltyWeight( iter );
                    
                penalty.setFactor( factor );
                    
                surfaceForceInt( SURFACETUPLEBINDER::makeTuple( *iter ) );
            }
            
            return;
        }

        
        
        
    } // namespace asmb
} // namespace base


//------------------------------------------------------------------------------
/** Penalty method for the weak incorporation of Dirichlet boundary conditions.
 *  Consider a boundary value problem with solution \f$ u \f$ which has the
 *  Dirichlet boundary condition
 *  \f[
 *          u = \bar{u}  \quad x \in \Gamma_D
 *  \f]
 *  where \f$ \Gamma_D \f$ is the Dirichlet boundary. A weak incorporation of
 *  this condition via a penalty method leads to the surface integral
 *  \f[
 *         p(w,v) = \frac{\gamma}{h} \int_{\Gamma_D} w v d s
 *  \f]
 *  with a scalar multiplier \f$ \gamma \f$, a size measure of the domain mesh
 *  \f$ h \f$ and functions \f$ w \f$ and \f$ v \f$.
 *  With the penalty method, the weak statement of the boundary value problem
 *  is augmented by
 *  \f[
 *         p(u - \bar{u},v)
 *  \f]
 *  This class implements these two terms, one contributing to the system
 *  matrix and the other to the system force vector.
 *
 *  Note that simply augmenting the weak form by above terms, leads effectively
 *  to the Robin boundary condition
 *  \f[
 *       \frac{h}{\gamma} t(u) + u = \bar{u}
 *  \f]
 *  with the \a traction \f$ t \f$. This condition, as \f$ h \to 0 \f$ and/or
 *   \f$ \gamma \to \infty \f$ recovers the original boundary condition. But
 *  since both limits are never actually attained, the spurious term remains
 *  and optimal convergence is lost.
 *  
 *  \tparam SURFFIELDTUPLE  Tuple of surface element bound to field elements
 */
template<typename SURFFIELDTUPLE>
class base::nitsche::Penalty
{
public:
    //! Template parameter
    typedef SURFFIELDTUPLE SurfFieldTuple;

    //! @name Extract element types
    //@{
    typedef typename SurfFieldTuple::GeomElement  SurfaceElement;
    typedef typename SurfFieldTuple::TestElement  TestElement;
    typedef typename SurfFieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim GlobalVecDim;

    //! Type of domain element
    typedef typename SurfaceElement::DomainElement    DomainElement;

    //! Size of a DoF
    static const unsigned doFSize = TestElement::DegreeOfFreedom::size;

    //! Type of DoF value vector
    typedef typename base::Vector<doFSize,base::number>::Type VecDoF;

    //--------------------------------------------------------------------------
    Penalty( const double factor ) : factor_( factor ) { }

    void setFactor( const double factor ) { factor_ = factor; }
        

    //--------------------------------------------------------------------------
    /** Evaluation of the kernel function for the system matrix term
     *  \f$ p(u,v) \f$.
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[out]  result Result container (pre-sized and zero-initialised)
     */
    void tangentStiffness( const SurfFieldTuple& surfFieldTuple,
                           const LocalVecDim&    eta,
                           const double          weight,
                           base::MatrixD&        result ) const
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();
        const TestElement*    testEp  = surfFieldTuple.testElementPtr();
        const TrialElement*   trialEp = surfFieldTuple.trialElementPtr();

        // Get pointer to domain element
        const DomainElement* domainEp = surfEp -> getDomainElementPointer();

        // compute mesh size
        const double h = base::mesh::Size<DomainElement>::apply( domainEp );

        // Get surface metric
        GlobalVecDim dummy;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, dummy );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );
                
        // Evaluate the shape functions
        typename TestElement::FEFun::FunArray testFunValues;
        (testEp -> fEFun()).evaluate( domainEp, xi, testFunValues );

        typename TrialElement::FEFun::FunArray trialFunValues;
        (trialEp -> fEFun()).evaluate( domainEp, xi, trialFunValues );

        // deduce the size of every contribution
        const unsigned numRowBlocks = static_cast<unsigned>( testFunValues.size()  );
        const unsigned numColBlocks = static_cast<unsigned>( trialFunValues.size() );

        // scalar multiplier of the whole entry
        const double aux = weight * detG * (factor_ / h);

        // Loop over shape functions
        for ( unsigned i = 0; i < numRowBlocks; i++ ) {

            for ( unsigned j = 0; j < numColBlocks; j++ ) {

                const double entry =
                    testFunValues[i] * trialFunValues[j] * aux;
                
                for ( unsigned d = 0; d < doFSize; d++ ) {
                    result( i*doFSize + d, j*doFSize + d ) += entry;
                }
                
            }
        }

        return;
    }

    //--------------------------------------------------------------------------
    /** Compute the residual forces due to a Penalty method for BC application.
     *  The residual force term due to the the penalty method has the form
     *  \f[
     *     R_{pen} = -\frac{\gamma}{h} \int_{\Gamma_D} (u-\bar{u}) \cdot v ds
     *  \f]
     *  and this function provides the integral kernel of this term.
     *  \tparam BCFUN Type of BC evaluation policy
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[in]   bcFun  Function describing the boundary condition
     *  \param[out]  result Result container (pre-sized and zero-initialised)
     */
    template<typename EVALPOL>
    void residualBoundary( const SurfFieldTuple& surfFieldTuple,
                           const LocalVecDim& eta,
                           const double       weight,
                           const typename EVALPOL::Fun& bcFun,
                           base::VectorD&        result ) const
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();
        const TrialElement*   trialEp = surfFieldTuple.trialElementPtr();
        
        // Get pointer to domain element
        const DomainElement* domainEp = surfEp -> getDomainElementPointer();

        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );
                
        // Evaluate boundary coondition
        const VecDoF bc = EVALPOL::apply( domainEp, xi, bcFun );
        const VecDoF  u = base::post::evaluateField( domainEp, trialEp, xi );
        const VecDoF  funResidual = u - bc;

        this -> evaluateResidual_( surfFieldTuple, eta, weight, funResidual, result );

        return;
    }

    //--------------------------------------------------------------------------
    /** Compute the residual forces due to a Penalty method for interfaces.
     *  The residual force term due to the the penalty method has the form
     *  \f[
     *     R_{pen} = -\frac{\gamma}{h} \int_{\Gamma_D} u \cdot v ds
     *  \f]
     *  and this function provides the integral kernel of this term.
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[in]   bcFun  Function describing the boundary condition
     *  \param[out]  result Result container (pre-sized and zero-initialised)
     */
    void residualInterface( const SurfFieldTuple& surfFieldTuple,
                            const LocalVecDim&    eta,
                            const double          weight,
                            base::VectorD&        result ) const
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();
        const TrialElement*   trialEp = surfFieldTuple.trialElementPtr();
        
        // Get pointer to domain element
        const DomainElement* domainEp = surfEp -> getDomainElementPointer();

        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );
                
        // Evaluate boundary coondition
        const VecDoF  u = base::post::evaluateField( domainEp, trialEp, xi );

        this -> evaluateResidual_( surfFieldTuple, eta, weight, u, result );

        return;
    }

private:

    // Helper to reduce implementation
    void evaluateResidual_( const SurfFieldTuple& surfFieldTuple,
                            const LocalVecDim& eta,
                            const double       weight,
                            const VecDoF&      funResidual,
                            base::VectorD&     result ) const
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();
        const TestElement*    testEp  = surfFieldTuple.testElementPtr();
        
        // Get pointer to domain element
        const DomainElement* domainEp = surfEp -> getDomainElementPointer();

        // compute mesh size
        const double h = base::mesh::Size<DomainElement>::apply( domainEp );

        // Get surface metric
        GlobalVecDim dummy;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, dummy );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );
                
        // Evaluate the shape function
        typename TestElement::FEFun::FunArray testFunValues;
        (testEp -> fEFun()).evaluate( domainEp, xi, testFunValues );

        // deduce the size of every contribution
        const unsigned numRowBlocks = static_cast<unsigned>( testFunValues.size() );

        // scalar multiplier of the whole entry
        const double aux = weight * detG * (factor_ / h);
        
        // Loop over shape functions
        for ( unsigned i = 0; i < numRowBlocks; i++ ) {

            for ( unsigned d = 0; d < doFSize; d++ ) {

                result[ i*doFSize + d ] -= testFunValues[i] * funResidual[d] * aux;
            }
        }
    }

private:
    double factor_;
};


#endif

