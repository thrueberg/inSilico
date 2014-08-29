//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Energy.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_nitsche_energy_hpp
#define base_nitsche_energy_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/asmb includes
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
// base/auxi includes
#include <base/auxi/FunEvaluationPolicy.hpp>
// base/post includes
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace nitsche{

        template<typename KERNEL, typename SURFFIELDTUPLE>
        class Energy;

        //----------------------------------------------------------------------
        /** Convenience functions for the LHS terms of the penalty method */

        // PRIMAL
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL, typename SURFACEQUADRATURE,
                 typename SOLVER, typename BOUNDFIELD,
                 typename PARAMETER>
        void primalEnergyLHS( const KERNEL&            kernel, 
                              const SURFACEQUADRATURE& surfaceQuadrature,
                              SOLVER&                  solver, 
                              const BOUNDFIELD&        boundField,
                              const PARAMETER&         parameter, 
                              const bool               inOut = true,
                              const bool               plusMinus = true )
        {
            // object to compute the LHS penalty term
            typedef base::nitsche::Energy<KERNEL,
                                          typename SURFACETUPLEBINDER::Tuple> Energy;
            Energy energy( kernel );
            
            // primal
            typedef base::asmb::StiffnessMatrix<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::Tuple> SysMat;

            typename SysMat::Kernel kernelFunP =
                boost::bind( &Energy::primal, &energy, _1, _2, _3, _4 );
            
            SysMat primal( kernelFunP, surfaceQuadrature, solver, false );
            
            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {

                const double sign  = ( plusMinus ? 1.0 : -1.0 );
                const double kappa = sign * parameter.energyWeight( iter, inOut );
                energy.setKappa( kappa );

                primal( SURFACETUPLEBINDER::makeTuple(           *iter ) );
            }

            return;
        }

        //----------------------------------------------------------------------
        // DUAL
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL, typename SURFACEQUADRATURE,
                 typename SOLVER, typename BOUNDFIELD,
                 typename PARAMETER>
        void dualEnergyLHS( const KERNEL&            kernel, 
                            const SURFACEQUADRATURE& surfaceQuadrature,
                            SOLVER&                  solver, 
                            const BOUNDFIELD&        boundField,
                            const PARAMETER&         parameter, 
                            const bool               inOut = true,
                            const bool               plusMinus = true )
        {
            // object to compute the LHS penalty term
            typedef base::nitsche::Energy<KERNEL,
                                          typename SURFACETUPLEBINDER::Tuple> Energy;
            Energy energy( kernel );
            
            // dual
            typedef base::asmb::StiffnessMatrix<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::TransposedTuple>
                SysMatT;

            typename SysMatT::Kernel kernelFunD =
                boost::bind( &Energy::dual, &energy, _1, _2, _3, _4 );
            
            SysMatT dual( kernelFunD, surfaceQuadrature, solver, false );

            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {

                const double sign  = ( plusMinus ? 1.0 : -1.0 );
                const double kappa = sign * parameter.energyWeight( iter, inOut );
                energy.setKappa( kappa );

                dual(   SURFACETUPLEBINDER::makeTransposedTuple( *iter ) );
            }

            return;
        }

        //----------------------------------------------------------------------
        // BOTH
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL, typename SURFACEQUADRATURE,
                 typename SOLVER, typename BOUNDFIELD,
                 typename PARAMETER>
        void energyLHS( const KERNEL&            kernel, 
                        const SURFACEQUADRATURE& surfaceQuadrature,
                        SOLVER&                  solver, 
                        const BOUNDFIELD&        boundField,
                        const PARAMETER&         parameter, 
                        const bool               inOut = true,
                        const bool               plusMinus = true )
        {
            primalEnergyLHS<SURFACETUPLEBINDER>(
                kernel, surfaceQuadrature, solver, boundField,
                parameter, inOut, plusMinus );
            dualEnergyLHS<SURFACETUPLEBINDER>(
                kernel, surfaceQuadrature, solver, boundField,
                parameter, inOut, plusMinus );
        }

        //----------------------------------------------------------------------
        namespace detail_{
            template<typename SURFACETUPLEBINDER, 
                     typename EVALUATIONPOLICY,
                     typename KERNEL, typename SURFACEQUADRATURE,
                     typename SOLVER, typename BOUNDFIELD,    
                     typename PARAMETER>
            void computeEnergyRHS( const SURFACEQUADRATURE& surfaceQuadrature,
                                   const KERNEL& kernel, 
                                   SOLVER& solver, BOUNDFIELD& boundField,
                                   const typename EVALUATIONPOLICY::Fun& bcFun,
                                   const PARAMETER& parameter,
                                   const bool inOut,
                                   const bool plusMinus )
            {
                // object to compute the LHS penalty term
                typedef base::nitsche::Energy<KERNEL,
                                              typename SURFACETUPLEBINDER::Tuple> Energy;
                Energy energy( kernel );

                // integrator and assembler object
                typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                                    typename SURFACETUPLEBINDER::TransposedTuple>
                    SurfaceForceInt;
                typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                    boost::bind( &Energy::template rhs<EVALUATIONPOLICY>, &energy,
                                 _1, _2, _3, boost::ref( bcFun ), _4 );
            
                SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                                 surfaceQuadrature, solver );

                // apply
                typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
                typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
                for ( ; iter != end; ++iter ) {
                    
                    const double sign  = (plusMinus ? 1.0 : -1.0 );
                    const double kappa = sign * parameter.energyWeight( iter, inOut );
                    energy.setKappa( kappa );

                    surfaceForceInt( SURFACETUPLEBINDER::makeTransposedTuple( *iter ) );
                }

                return;

            }
        }
        //----------------------------------------------------------------------

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL,     typename SURFACEQUADRATURE,
                 typename SOLVER,     typename BOUNDFIELD,
                 typename BCFUN,      typename PARAMETER>
        void energyRHS( const KERNEL&            kernel,
                        const SURFACEQUADRATURE& surfaceQuadrature,
                        SOLVER&                  solver, 
                        const BOUNDFIELD&        boundField,
                        const BCFUN&             bcFun,
                        const PARAMETER&         parameter,
                        const bool               inOut = true,
                        const bool               plusMinus = true )
        {
            typedef typename SURFACETUPLEBINDER::Tuple::GeomElement::DomainElement
                DomainElement;
            typedef base::auxi::EvaluateDirectly<DomainElement,BCFUN>
                EvaluationPolicy;

            detail_::computeEnergyRHS<SURFACETUPLEBINDER,
                                      EvaluationPolicy>( surfaceQuadrature, kernel,
                                                         solver, boundField, bcFun,
                                                         parameter,
                                                         inOut, plusMinus );

            return;
        }

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL,     typename SURFACEQUADRATURE,
                 typename SOLVER,     typename BOUNDFIELD,
                 typename BCFUN,      typename PARAMETER>
        void energyRHS2( const KERNEL&            kernel,
                         const SURFACEQUADRATURE& surfaceQuadrature,
                         SOLVER&                  solver, 
                         const BOUNDFIELD&        boundField,
                         const BCFUN&             bcFun,
                         const PARAMETER&         parameter,
                         const bool               inOut = true,
                         const bool               plusMinus = true )
        {
            typedef typename SURFACETUPLEBINDER::Tuple::GeomElement::DomainElement
                DomainElement;
            typedef base::auxi::EvaluateViaElement<DomainElement,BCFUN>
                EvaluationPolicy;

            detail_::computeEnergyRHS<SURFACETUPLEBINDER,
                                      EvaluationPolicy>( surfaceQuadrature, kernel,
                                                         solver, boundField, bcFun,
                                                         parameter,
                                                         inOut, plusMinus );

            return;
        }

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL,     typename SURFACEQUADRATURE,
                 typename SOLVER,     typename BOUNDFIELD,
                 typename PARAMETER>
        void energyResidual( const KERNEL&            kernel,
                             const SURFACEQUADRATURE& surfaceQuadrature,
                             SOLVER&                  solver, 
                             const BOUNDFIELD&        boundField,
                             const PARAMETER&         parameter,
                             const bool               inOut = true,
                             const bool               plusMinus = true )
        {
            // object to compute the LHS penalty term
            typedef Energy<KERNEL,
                           typename SURFACETUPLEBINDER::Tuple> Energy;
            Energy energy( kernel );

            // integrator and assembler object
            typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::Tuple>
                SurfaceForceInt;
            typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                boost::bind( &Energy::residual, &energy, _1, _2, _3, _4 );
            
            SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                             surfaceQuadrature, solver );
            
            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {

                const double sign  = (plusMinus ? 1.0 : -1.0 );
                const double kappa = sign * parameter.energyWeight( iter, inOut );
                energy.setKappa( kappa );

                surfaceForceInt( SURFACETUPLEBINDER::makeTuple( *iter ) );
            }

            return;
        }

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL,     typename SURFACEQUADRATURE,
                 typename SOLVER,     typename BOUNDFIELD,
                 typename PARAMETER>
        void energyResidualTransposed( const KERNEL&            kernel,
                                       const SURFACEQUADRATURE& surfaceQuadrature,
                                       SOLVER&                  solver, 
                                       const BOUNDFIELD&        boundField,
                                       const PARAMETER&         parameter,
                                       const bool               inOut = true,
                                       const bool               plusMinus = true )
        {
            // object to compute the LHS penalty term
            typedef Energy<KERNEL,
                           typename SURFACETUPLEBINDER::Tuple> Energy;
            Energy energy( kernel );

            // integrator and assembler object
            typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::TransposedTuple>
                SurfaceForceInt;
            typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                boost::bind( &Energy::transposedLinearisedResidual,
                             &energy, _1, _2, _3, _4 );
            
            SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                             surfaceQuadrature, solver );
            
            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {

                const double sign  = (plusMinus ? 1.0 : -1.0 );
                const double kappa = sign * parameter.energyWeight( iter, inOut );
                energy.setKappa( kappa );

                surfaceForceInt( SURFACETUPLEBINDER::makeTransposedTuple( *iter ) );
            }

            return;
        }

    }
}

//------------------------------------------------------------------------------
/** Computation of the variationally consistent terms of the Nitsche Method.
 *  Consider, as an example, the problem
 *  \f[
 *       -\nabla^2 u = 0, x \in \Omega
 *  \f]
 *  with Dirichlet boundary condition
 *  \f[
 *        u = g, x \in \Gamma
 *  \f]
 *  Taking the residual of the PDE, multiplied by a suitable test function
 *  \f$ v \f$ and integrated by parts gives
 *  \f[
 *      0 = \int_\Omega -\nabla^2 u v d x
 *        = \int_\Omega \nabla u \nabla v d x - \int_\Gamma v(\nabla u) n d s
 *  \f]
 *  Now, zero is added in the form of \f$ 0 = u - g \f$ and multiplied by
 *  the normal-derivative of \f$ v \f$,
 *  \f[
 *      \int_\Omega \nabla u \nabla v d x
 *     -\int_\Gamma v (\nabla u) n + u (\nabla v ) n =
 *     -\int_\Gamma g (\nabla v) n
 *  \f]
 *  This gives rise to a (yet unstable) method to compute the Dirichlet
 *  problem without essential boundary conditions. A penalty term is used
 *  to stablised this formulation.
 *  This object provides the additional terms on the left- and right-hand
 *  sides of the above expression.
 *  \tparam KERNEL  Type of kernel giving the bilinear form and normal
 *                  derivatives
 *  \tparam SURFFIELDTUPLE Tuple with surface element
 * 
 */
template<typename KERNEL, typename SURFFIELDTUPLE>
class base::nitsche::Energy
{
public:
    //! @name Template parameter
    //@{
    typedef KERNEL         Kernel;
    typedef SURFFIELDTUPLE SurfFieldTuple;
    //@}

    //! @name Extract element types
    //@{
    typedef typename SurfFieldTuple::GeomElement  SurfaceElement;
    typedef typename SurfFieldTuple::TestElement  TestElement;
    typedef typename SurfFieldTuple::TrialElement TrialElement;
    //@}

    //! Type of field element tuple generated from the surface tuple
    typedef typename
    base::asmb::DomainFieldElementPointerTuple<SurfFieldTuple>::Type DomainFieldTuple;


    typedef typename SurfFieldTuple::TransposedTuple TransposedSurfFieldTuple;

    typedef typename
    base::asmb::DomainFieldElementPointerTuple<TransposedSurfFieldTuple>::Type
    TransposedDomainFieldTuple;

    //! Local coordinate
    typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim GlobalVecDim;

    //! Type of domain element
    typedef typename SurfaceElement::DomainElement    DomainElement;

    //! Size of a DoF
    static const unsigned doFSize = TestElement::DegreeOfFreedom::size;

    //! Type of DoF value vector
    typedef typename base::Vector<doFSize,base::number>::Type VecDof;

    //! Type of BC function
    typedef boost::function<VecDof( const GlobalVecDim& )> BCFun;

    //! Constructor with kernel object and multiplier
    Energy( const Kernel& kernel, const double kappa = 1.0 )
        : kernel_( kernel ), kappa_( kappa ) { }

    void setKappa( const double kappa ) { kappa_ = kappa; }

    //--------------------------------------------------------------------------
    /** Compute the system matrix contribution from the boundary energy.
     *  This primal term is the very one arising from integration by parts of
     *  the weighted residual. Let \f$ B_N u \f$ denote the boundary operator
     *  (possibly linearised), this function computes the discrete counterpart
     *  of
     *  \f[
     *        - \int_\Gamma (B_N u) v ds
     *  \f]
     *  The implementation is deferred to the function lhsHelper_() in order to
     *  share the functions with the dual() which computes the transposed of
     *  the above turn, i.e. \f$ u \f$ and \f$ v \f$ interchanged.
     *
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[out]  result Result container (pre-sized and zero-initialised)
     */
    void primal( const SurfFieldTuple& surfFieldTuple,
                 const LocalVecDim&    eta,
                 const double          weight,
                 base::MatrixD&        result )
    {
        // extract surface geometry element
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();

        // construct a field element pointer tuple with the domain elements
        const DomainFieldTuple domainFieldTuple =
            base::asmb::DomainFieldElementPointerTuple<SurfFieldTuple>::
            convert( surfFieldTuple );

        // call implementation
        const bool isDual = false;
        this -> lhsHelper_( domainFieldTuple, surfEp, eta, weight,
                            isDual, result );
    }
    

    //--------------------------------------------------------------------------
    /** Transposed energy terms due to Nitsche's method.
     *  In Nitsche's method, the Dirichlet boundary condition is incorporated
     *  by weak imposition weighted with the \em dual co-normal derivative
     *  \f[
     *       - \int_\Gamma (B_N v)( u - \bar{u} ) ds
     *  \f]
     *  Here, only the system matrix contribution is computed in its discrete
     *  form. In order to reduce the implementation, this function is called
     *  with the transposed tuple (i.e. test and trial fields interchanged)
     *  and in this function this transposition is reversed and the function
     *  lhsHelper_() called accordingly.
     *
     *  \param[in]  surfFieldTupleT Transposed tuple of field element pointers 
     *  \param[in]  eta    Local evaluation coordinate
     *  \param[in]  weight Corresponding quadrature weight
     *  \param[out] result Result container (pre-sized and zero-initialised)
     */
    void dual( const TransposedSurfFieldTuple& surfFieldTupleT,
               const LocalVecDim&              eta,
               const double                    weight,
               base::MatrixD&                  result )
    {
        // extract surface geometry element of the tuple
        const SurfaceElement* surfEp  = surfFieldTupleT.geomElementPtr();

        // extract domain field tuple (also transposed)
        const TransposedDomainFieldTuple domainFieldTupleT =
            base::asmb::DomainFieldElementPointerTuple<TransposedSurfFieldTuple>::
            convert( surfFieldTupleT );

        // transpose the domain field tuple in order to recover the normal one
        const DomainFieldTuple domainFieldTuple =
            domainFieldTupleT.transpose();

        // call implementation
        const bool isDual = true;
        this -> lhsHelper_( domainFieldTuple, surfEp, eta, weight,
                            isDual, result );
    }

private:
    //--------------------------------------------------------------------------
    /** Helper implementation.
     *  Evaluates the test field functions, evaluates the co-normal
     *  derivative of the kernel object and assembles into the result container.
     *  In the dual case, the test-field is actually the trial field of the
     *  transposed constellation and the result is assmebled in a transposed
     *  manner.
     *  \param[in]  domainFieldTuple  Tuple of domain elements
     *  \param[in]  surfEp            Pointer to surface geometry element
     *  \param[in]  eta               Local surface coordinate
     *  \param[in]  weight            Quadrature weight
     *  \param[in]  isDual            Flag for the dual case
     *  \param[out] result            Output
     */
    void lhsHelper_( const DomainFieldTuple&  domainFieldTuple,
                     const SurfaceElement*    surfEp,
                     const LocalVecDim&       eta,
                     const double             weight,
                     const bool               isDual, 
                     base::MatrixD&           result )
    {
        // Get pointer to domain element
        const DomainElement* domainEp = surfEp -> getDomainElementPointer();

        // Get surface metric
        GlobalVecDim normal;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, normal );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );

        // compute co-normal from kernel using domain tuple
        base::MatrixD coNormal;
        kernel_.coNormalDerivative( domainFieldTuple,
                                    xi,
                                    normal, coNormal );

        // evaluate test and trial functions
        typename TestElement::FEFun::FunArray testFunValues;
        (domainFieldTuple.testElementPtr() -> fEFun()).evaluate( domainEp, xi, testFunValues );

        // compute matrix entries
        const unsigned numFieldBlocks = static_cast<unsigned>(testFunValues.size() );
        const unsigned otherSize      = static_cast<unsigned>(coNormal.cols() );

        const double scalar = -1. * detG * weight * kappa_;

        for ( unsigned i = 0; i < numFieldBlocks; i++ ) {
            for ( unsigned d = 0; d < doFSize; d++ ) {
                for ( unsigned c = 0; c < otherSize; c++ ) {

                    if ( isDual ) 
                        result( c, i*doFSize+d ) +=
                            testFunValues[i] * coNormal(d,c) * scalar;
                    else 
                        result( i*doFSize+d, c ) +=
                            testFunValues[i] * coNormal(d,c) * scalar;
                    
                }
            }
        }

        return;
    }

public:
    //--------------------------------------------------------------------------
    /**
     *  \tparam BCFUN Type of BC function as a function of \f$ x \f$
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[in]   bcFun  Function describing the boundary condition
     *  \param[out]  result Result container (pre-sized and zero-initialised)
     */
    template<typename EVALPOL>
    void rhs( const TransposedSurfFieldTuple& surfFieldTuple, 
              const LocalVecDim&              eta,
              const double                    weight,
              const typename EVALPOL::Fun&    bcFun, 
              base::VectorD&                  result )
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();

        // Get surface metric
        GlobalVecDim normal;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, normal );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );

        // get domain field element pointer tuple (still transposed)
        const TransposedDomainFieldTuple domainFieldTupleT =
            base::asmb::DomainFieldElementPointerTuple<TransposedSurfFieldTuple>::
            convert( surfFieldTuple );

        // un-transpose the tuple
        const DomainFieldTuple domainFieldTuple =
            domainFieldTupleT.transpose();

        // compute co-normal from kernel using domain tuple
        base::MatrixD coNormal;
        kernel_.coNormalDerivative( domainFieldTuple, xi,
                                    normal, coNormal );

        // Evaluate boundary coondition
        const VecDof bc = EVALPOL::apply( surfEp -> getDomainElementPointer(),
                                          xi, bcFun );

        // Evaluate solution
        const VecDof  u = base::post::evaluateField( surfEp -> getDomainElementPointer(),
                                                     surfFieldTuple.trialElementPtr(),
                                                     xi );

        // scalar multiplier
        const double scalar = -1. * detG * weight * kappa_;

        const unsigned otherSize = static_cast<unsigned>(coNormal.cols() );
        
        for ( unsigned i = 0; i < otherSize; i++ ) {
            number sum = 0.;
            for ( unsigned d = 0; d < doFSize; d++ ) {
                sum -= scalar * coNormal( d, i ) * (u[d]-bc[d]);
            }
            result[i] += sum;
        }

        return;
    }

public:
    //--------------------------------------------------------------------------
    /**
     *  \tparam BCFUN Type of BC function as a function of \f$ x \f$
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[out]  result Result container (pre-sized and zero-initialised)
     */
    void residual( const SurfFieldTuple& surfFieldTuple,
                   const LocalVecDim&    eta,
                   const double          weight,
                   base::VectorD&        result )
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();

        // Get surface metric
        GlobalVecDim normal;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, normal );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );

        // convert tupe to domain tuple
        const DomainFieldTuple domainFieldTuple =
            base::asmb::DomainFieldElementPointerTuple<SurfFieldTuple>::
            convert( surfFieldTuple );

        // compute co-normal from kernel using domain tuple
        base::VectorD coNormalResidual;
        kernel_.boundaryResidual( domainFieldTuple, xi,
                                  normal, coNormalResidual );

        //
        const double scalar = detG * weight * kappa_;

        result += scalar * coNormalResidual;

        return;
    }

    //--------------------------------------------------------------------------
    /**
     *  \tparam BCFUN Type of BC function as a function of \f$ x \f$
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[out]  result Result container (pre-sized and zero-initialised)
     */
    void transposedLinearisedResidual( const TransposedSurfFieldTuple& surfFieldTuple,
                                       const LocalVecDim&    eta,
                                       const double          weight,
                                       base::VectorD&        result )
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();

        // Get surface metric
        GlobalVecDim normal;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, normal );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );

        // get domain field element pointer tuple (still transposed)
        const TransposedDomainFieldTuple domainFieldTupleT =
            base::asmb::DomainFieldElementPointerTuple<TransposedSurfFieldTuple>::
            convert( surfFieldTuple );

        // un-transpose the tuple
        const DomainFieldTuple domainFieldTuple =
            domainFieldTupleT.transpose();

        // compute co-normal from kernel using domain tuple
        base::MatrixD coNormal;
        kernel_.coNormalDerivative( domainFieldTuple, xi,
                                    normal, coNormal );

        // Evaluate solution
        const VecDof  u = base::post::evaluateField( surfEp -> getDomainElementPointer(),
                                                     surfFieldTuple.trialElementPtr(),
                                                     xi );

        // scalar multiplier
        const double scalar = -1. * detG * weight * kappa_;

        const unsigned otherSize = static_cast<unsigned>(coNormal.cols() );
        
        for ( unsigned i = 0; i < otherSize; i++ ) {
            double sum = 0.;
            for ( unsigned d = 0; d < doFSize; d++ ) {
                sum += scalar * coNormal( d, i ) * u[d];
            }
            result[i] += sum;
        }

        return;
    }

private:
    const Kernel&  kernel_; //!< Delivers conormal and dual conormal derivative
    double         kappa_;  //!< Factor for interface binding
};

#endif
