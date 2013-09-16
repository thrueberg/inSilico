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

//------------------------------------------------------------------------------
namespace base{
    namespace nitsche{

        template<typename SURFFIELDTUPLE>
        class Penalty;

        //----------------------------------------------------------------------
        /** Convenience function for the LHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD>
        void penaltyLHS( const SURFACEQUADRATURE& surfaceQuadrature,
                         SOLVER&                  solver, 
                         const BOUNDFIELD&        boundField,
                         const double factor )
        {

            // object to compute the LHS penalty term
            typedef base::nitsche::Penalty<typename SURFACETUPLEBINDER::Tuple>
                Penalty;
            Penalty penalty( factor );

            // integrator and assembler object
            typedef base::asmb::StiffnessMatrix<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::Tuple> SysMat;
            
            typename SysMat::Kernel kernel =
                boost::bind( &Penalty::lhs, &penalty, _1, _2, _3, _4 );

            SysMat systemMat( kernel, surfaceQuadrature, solver );

            // apply
            typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
            typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
            for ( ; iter != end; ++iter ) {
                systemMat( SURFACETUPLEBINDER::makeTuple( *iter ) );
            }
            
            return;
        }


        //----------------------------------------------------------------------
        /** Convenience function for the RHS term of the penalty method */
        template<typename SURFACETUPLEBINDER,
                 typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD, typename BCFUN>
        void penaltyRHS( const SURFACEQUADRATURE& surfaceQuadrature,
                         SOLVER&                  solver, 
                         const BOUNDFIELD&        boundField,
                         const BCFUN&             bcFun,
                         const double factor )
        {
            // object to compute the LHS penalty term
            typedef Penalty<typename SURFACETUPLEBINDER::Tuple> Penalty;
            Penalty penalty( factor );
            
            // integrator and assembler object
            typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                                typename SURFACETUPLEBINDER::Tuple>
                SurfaceForceInt;
            typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                boost::bind( &Penalty::template rhs<BCFUN>, &penalty, _1, _2, _3,
                             boost::ref( bcFun ), _4 );
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
 *         p(u,v) - p(\bar{u},v)
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
    typedef typename base::Vector<doFSize,base::number>::Type VecDof;

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
    void lhs( const SurfFieldTuple& surfFieldTuple,
              const LocalVecDim&    eta,
              const double          weight,
              base::MatrixD&        result )
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();
        const TestElement*    testEp  = surfFieldTuple.testElementPtr();
        const TrialElement*   trialEp = surfFieldTuple.trialElementPtr();

        // check for bubnov case via pointer identity
        const bool isBubnov = (testEp == trialEp);
        
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
        if ( isBubnov ) trialFunValues = testFunValues;
        else (trialEp -> fEFun()).evaluate( domainEp, xi, trialFunValues );

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
    /** Evaluation of the kernel function for the system force term
     *  \f$ p( \bar{u}, v ) \f$.
     *
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

        // Evaluate element geometry
        const GlobalVecDim x = base::Geometry<SurfaceElement>()( surfEp, eta );

        // Evaluate boundary coondition
        const VecDof bc = bcFun( x );

        // scalar multiplier of the whole entry
        const double aux = weight * detG * (factor_ / h);

        // Loop over shape functions
        for ( unsigned i = 0; i < numRowBlocks; i++ ) {

            for ( unsigned d = 0; d < doFSize; d++ ) {

                result[ i*doFSize + d ] +=
                    testFunValues[i] * bc[d] * aux;
            }
        }

        return;
    }

private:
    double factor_;
};


#endif

