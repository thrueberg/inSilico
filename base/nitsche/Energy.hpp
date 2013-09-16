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

//------------------------------------------------------------------------------
namespace base{
    namespace nitsche{

        template<typename KERNEL, typename SURFFIELDTUPLE>
        class Energy;

        //----------------------------------------------------------------------
        /** Convenience function for the LHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL, typename SURFACEQUADRATURE,
                 typename SOLVER, typename BOUNDFIELD>
        void energyLHS( const KERNEL&            kernel, 
                        const SURFACEQUADRATURE& surfaceQuadrature,
                        SOLVER&                  solver, 
                        const BOUNDFIELD&        boundField,
                        const double kappa = 1.0 )
        {
            // object to compute the LHS penalty term
            typedef Energy<KERNEL,
                           typename SURFACETUPLEBINDER::Tuple> Energy;
            Energy energy( kernel, kappa );

            // integrator and assembler objects
            typedef base::asmb::StiffnessMatrix<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::Tuple> SysMat;

            typedef base::asmb::StiffnessMatrix<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::TransposedTuple>
                SysMatT;

            // primal 
            typename SysMat::Kernel kernelFunP =
                boost::bind( &Energy::primal, &energy, _1, _2, _3, _4 );

            SysMat primal( kernelFunP, surfaceQuadrature, solver );

            // dual
            typename SysMatT::Kernel kernelFunD =
                boost::bind( &Energy::dual, &energy, _1, _2, _3, _4 );
            
            SysMatT dual( kernelFunD, surfaceQuadrature, solver );

            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {
                primal( SURFACETUPLEBINDER::makeTuple(           *iter ) );
                dual(   SURFACETUPLEBINDER::makeTransposedTuple( *iter ) );
            }

            
            return;
        }

        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename KERNEL, typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD, typename BCFUN>
        void energyRHS( const KERNEL&            kernel,
                        const SURFACEQUADRATURE& surfaceQuadrature,
                        SOLVER&                  solver, 
                        const BOUNDFIELD&        boundField,
                        const BCFUN&             bcFun )
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
                boost::bind( &Energy::template rhs<BCFUN>, &energy,
                             _1, _2, _3, boost::ref( bcFun ), _4 );
            
            SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                             surfaceQuadrature, solver );
            
            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {
                surfaceForceInt( SURFACETUPLEBINDER::makeTuple( *iter ) );
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
    /** 
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
        const TestElement*    testEp  = surfFieldTuple.testElementPtr();
        
        this -> lhsHelper_( surfFieldTuple, eta, weight,
                            testEp, false, result );
    }

    void dual( const SurfFieldTuple& surfFieldTuple,
               const LocalVecDim&    eta,
               const double          weight,
               base::MatrixD&        result )
    {
        const TrialElement*    trialEp  = surfFieldTuple.trialElementPtr();
        
        this -> lhsHelper_( surfFieldTuple, eta, weight,
                            trialEp, true, result );
    }

private:
    template<typename FIELDELEMENT>
    void lhsHelper_( const SurfFieldTuple& surfFieldTuple,
                     const LocalVecDim&    eta,
                     const double          weight,
                     const FIELDELEMENT*   fieldEp,
                     const bool            dual, 
                     base::MatrixD&        result )
    {
        // convert tupe to domain tuple
        const DomainFieldTuple domainFieldTuple =
            base::asmb::DomainFieldElementPointerTuple<SurfFieldTuple>::
            convert( surfFieldTuple );

        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();

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
        if ( dual )            
            kernel_.dualCoNormalDerivative( domainFieldTuple, xi,
                                            normal, coNormal );
        else
            kernel_.coNormalDerivative( domainFieldTuple, xi,
                                        normal, coNormal );

        // evaluate test and trial functions
        typename FIELDELEMENT::FEFun::FunArray fieldFunValues;
        (fieldEp -> fEFun()).evaluate( domainEp, xi, fieldFunValues );

        // compute matrix entries
        const unsigned numFieldBlocks = static_cast<unsigned>(fieldFunValues.size() );
        const unsigned otherSize      = static_cast<unsigned>(coNormal.cols() );

        const double scalar = -1. * detG * weight * kappa_;

        for ( unsigned i = 0; i < numFieldBlocks; i++ ) {
            for ( unsigned d = 0; d < doFSize; d++ ) {
                for ( unsigned c = 0; c < otherSize; c++ ) {

                    if ( dual ) 
                        result( c, i*doFSize+d ) +=
                            fieldFunValues[i] * coNormal(d,c) * scalar;
                    else 
                        result( i*doFSize+d, c ) +=
                            fieldFunValues[i] * coNormal(d,c) * scalar;
                    
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
    template<typename BCFUN>
    void rhs( const SurfFieldTuple& surfFieldTuple,
              const LocalVecDim&    eta,
              const double          weight,
              const BCFUN&          bcFun,
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
        base::MatrixD coNormal;
        kernel_.dualCoNormalDerivative( domainFieldTuple, xi,
                                        normal, coNormal );

        // Evaluate element geometry
        const GlobalVecDim x = base::Geometry<SurfaceElement>()( surfEp, eta );

        // Evaluate boundary coondition
        const VecDof bc = bcFun( x );

        //
        const double scalar = -1. * detG * weight * kappa_;

        const unsigned otherSize = static_cast<unsigned>(coNormal.cols() );
        
        for ( unsigned i = 0; i < otherSize; i++ ) {
            number sum = 0.;
            for ( unsigned d = 0; d < doFSize; d++ ) {
                sum += scalar * coNormal( d, i ) * bc[d];
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
