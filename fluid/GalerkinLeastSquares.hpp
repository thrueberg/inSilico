//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   GalerkinLeastSquares.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef galerkinleastsquares_hpp
#define galerkinleastsquares_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
#include <base/mesh/Size.hpp>
// base/aux includes
#include <base/auxi/EqualPointers.hpp>
// base/post includes
#include <base/post/evaluateField.hpp>
#include <base/asmb/BodyForce.hpp>
//------------------------------------------------------------------------------
namespace fluid{
    namespace gls{
    
        template<typename T> class GalerkinLeastSquares;
        template<typename FIELDTUPLE> class StressDivergence;
        template<typename FIELDTUPLE> class VectorLaplace;
        template<typename FIELDTUPLE> class PressureGradient;
        template<typename FIELDTUPLE> class VelocityDivergence;
        template<typename FIELDTUPLE> class PressureGradient2;
        template<typename FIELDTUPLE> class VelocityDivergence2;
        template<typename FIELDTUPLE> class PressureLaplace;

        template<typename FIELDTUPLE, typename EVALUATEPOLICY> class BodyForce1;
        template<typename FIELDTUPLE, typename EVALUATEPOLICY> class BodyForce2;


        //----------------------------------------------------------------------
        // Computation of body force terms with given force function f(x)
        template<typename FIELDTUPLEBINDER1,
                 typename FIELDTUPLEBINDER2,
                 typename QUADRATURE,
                 typename SOLVER, typename FIELDBINDER, typename FUN>
        void bodyForceComputation(
            const QUADRATURE& quadrature,
            SOLVER& solver,
            const FIELDBINDER& fieldBinder,
            const FUN& forceFun,
            const double stabil, const double viscosity, const bool dw )
        {
            typedef base::auxi::EvaluateDirectly<
                typename FIELDTUPLEBINDER1::Tuple::GeomElement,FUN> Evaluate;
            
            // body force wrapper
            BodyForce1<typename FIELDTUPLEBINDER1::Tuple,Evaluate>
                bodyForce1( forceFun, stabil, viscosity, dw );
            BodyForce2<typename FIELDTUPLEBINDER2::Tuple,Evaluate>
                bodyForce2( forceFun, stabil );

            // force integrator 1
            typedef base::asmb::ForceIntegrator<QUADRATURE,SOLVER,
                                                typename FIELDTUPLEBINDER1::Tuple> ForceInt1;
            typename ForceInt1::ForceKernel forceKernel1 =
                boost::bind( bodyForce1, _1, _2, _3, _4 );
            ForceInt1 forceInt1( forceKernel1, quadrature, solver );

            // force integrator 2
            typedef base::asmb::ForceIntegrator<QUADRATURE,SOLVER,
                                                typename FIELDTUPLEBINDER2::Tuple> ForceInt2;
            typename ForceInt2::ForceKernel forceKernel2 =
                boost::bind( bodyForce2, _1, _2, _3, _4 );
            ForceInt2 forceInt2( forceKernel2, quadrature, solver );
            

            // Apply to all elements
            typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
            typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
            for ( ; iter != end; ++iter ) {
                forceInt1( FIELDTUPLEBINDER1::makeTuple( *iter ) );
                forceInt2( FIELDTUPLEBINDER2::makeTuple( *iter ) );
            }
            return;

        }

    }
}

//------------------------------------------------------------------------------
/** Stabilisation in order to circumvent the inf-sup condition.
 *  Stokes' system reads
 *  \f[
 *       -\nabla \cdot \sigma(u,p) = f, \quad \nabla \cdot u = 0
 *  \f]
 *  with the fluid Cauchy stress \f$\sigma = 2 \mu \varepsilon(u) - pI\f$
 *  base on the velocity field \f$ u \f$, the pressure \f$ p \f$ and
 *  the symmetric gradient \f$ \varepsilon(u) = \nabla u + (\nabla u)^T \f$.
 *  Alternatively, the momentum balance equation can be expressed as
 *  \f[
 *       -\mu \nabla^2 u = f
 *  \f]
 *  which is based on \f$ \nabla \cdot u = 0 \f$ but leads to different
 *  natural boundary conditions (!!).
 *  For notational simplicity, let \f$ L(u,p) = f \f$ represent the above system
 *  of PDEs and \f$ L^* \f$ be the adjoint to \f$ L \f$. More precisely,
 *  \f[
 *       L^*(u,p) = -2 \mu \nabla \cdot \varepsilon(u) + \nabla p
 *  \f]
 *  By means of the weighted residual method, one gets the bilinear form
 *  \f[
 *      b(u,p;v,q) = (2\mu \varepsilon(u), \varepsilon(v)) -
 *                   (p, \nabla \cdot v) - (q, \nabla \cdot u )
 *  \f]
 *  and for the RHS \f$ l(v) = (f,v) \f$.
 *
 *  A FE discretisation of the above system typically leads to a block-matrix
 *  system and the solvability of this is given by the inf-sup condition.
 *  In order to circumvent this condition and allow for a wider range of
 *  admissable shape functions, a stabilisation is added to the bilinear forms
 *  \f[
 *        b(u,p;v,q) - l(v) - \sum_K \alpha_K (L(u,p) - f, S(v,q) )_K = 0
 *  \f]
 *  This extra term is evaluated element-wise, \f$ K \f$ refers to the
 *  elements. Now, there are two choices for the operator \f$ S \f$
 *  \f[
 *        S(v,q) = L(v,q),  \quad S(v,q) = -L^*(v,q)
 *  \f]
 *  These choices correspond to
 *  -  Franca & Hughes (Galerkin Least Squares, GLS+), upper bound for
 *     \f$ \alpha \f$
 *  -  Douglas & Wang (sometimes GLS-), no upper bound for \f$ \alpha \f$
 *
 *  The implementation of the individual terms follows
 *  - \f$(2\mu \nabla \cdot \varepsilon(u),
 *        2\mu \nabla \cdot \varepsilon(v))_K\f$ in fluid::gls::StressDivergence
 *  - \f$ (\mu \nabla^2 u, \mu \nabla^2 v)_K \f$ in fluid::gls::VectorLaplace
 *  - \f$ (\mu \nabla^2 v, \nabla p )_K \f$ in fluid::gls::PressureGradient
 *  - \f$ (\mu \nabla \cdot \varepsilon(v),
 *         \nabla p )_K \f$ in fluid::gls::PressureGradient2 
 *  - \f$ (\mu \nabla^2 u, \nabla q )_K \f$ in fluid::gls::VelocityDivergence
 *  - \f$ (\mu \nabla \cdot \varepsilon(u),
 *         \nabla q )_K \f$ in fluid::gls::VelocityDivergence2
 *  - \f$ -(\nabla p, \nabla q)_K \f$ in fluid::gls::PressureLaplace
 *  - \f$ (2\mu \nabla \cdot \varepsilon(v), f)_K \f$ in fluid::gls::BodyForce1
 *  - \f$ (\nabla q, f)_K \f$ in fluid::gls::BodyForce2
 *
 *  All these terms receive a global stabilisation mulitplier \f$ \alpha \f$, 
 *  some a boolean for the choice between GLS+ and GLS- options.
 *
 */
template<typename T>
class fluid::gls::GalerkinLeastSquares
{
    STATIC_ASSERT_MSG( (sizeof(T) == 0), 
                       "This class exists for documentation purpose only" );

};

//------------------------------------------------------------------------------
/** Stress-divergence correspondent of GLS stabilisation.
 */
template<typename FIELDTUPLE>
class fluid::gls::StressDivergence
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    typedef typename base::Matrix<globalDim,globalDim>::Type        Hessian;

    //! Constructor with form and test functions 
    StressDivergence( const double stabil,
                      const double viscosity,
                      const bool douglasWang = false )
        : stabil_(      stabil ),
          viscosity_(   viscosity ),
          douglasWang_( douglasWang )
    { }

    //--------------------------------------------------------------------------
    void tangentStiffness( const FieldTuple&  fieldTuple,
                           const LocalVecDim& xi,
                           const double       weight,
                           base::MatrixD&     matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // if pointers are identical, Galerkin-Bubnov scheme
        const bool isBubnov =
            base::auxi::EqualPointers<TestElement,TrialElement>::apply( testEp,
                                                                       trialEp );

        // Evaluate gradient of test and trial functions
        std::vector<Hessian> testHessX, trialHessX;
        const double detJ =
            (testEp -> fEFun()).evaluateHessian( geomEp, xi, testHessX );
        
        if ( isBubnov ) trialHessX = testHessX;
        else
            (trialEp -> fEFun()).evaluateHessian( geomEp, xi, trialHessX );
        
        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testHessX.size()  );
        const unsigned numColBlocks = static_cast<unsigned>( trialHessX.size() );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * globalDim );
        assert( static_cast<unsigned>( matrix.cols() ) == numColBlocks * globalDim );

        // mesh size
        const double h = base::mesh::elementSize( geomEp );

        // scalar multiplier
        const double scalar =
            (douglasWang_ ? 1.0 : -1.0) *
            stabil_ * viscosity_ * viscosity_ * h * h * detJ * weight;
        
        //
        for ( unsigned M = 0; M < numRowBlocks; M ++ ) {
            for ( unsigned N = 0; N < numColBlocks; N ++ ) {

                // dot-product of the laplacian of the shape functions
                double laplaceTest  = 0.;
                double laplaceTrial = 0.;
                Hessian gradGrad = base::constantMatrix<globalDim,globalDim>( 0. );
                
                for ( unsigned k = 0; k < globalDim; k ++ ) {
                    laplaceTest  +=  testHessX[M](k,k);
                    laplaceTrial += trialHessX[N](k,k);

                    for ( unsigned i = 0; i < globalDim; i++ ) {
                        for ( unsigned j = 0; j < globalDim; j++ ) {
                            gradGrad(i,j) += testHessX[M](k,i) * trialHessX[N](k,j);
                        }
                    }
                }

                // add to matrix block's diagonal entries
                for ( unsigned i = 0; i < globalDim; i++ ) {
                    for ( unsigned j = 0; j < globalDim; j++ ) {

                        matrix( M*globalDim + i, N*globalDim + j ) +=
                            scalar * ( (i==j? laplaceTest * laplaceTrial : 0.) +
                                       laplaceTest * trialHessX[N](i,j) +
                                       testHessX[M](i,j) * laplaceTrial +
                                       gradGrad(i,j) );
                    }
                }
            }
        }

        return;
    }
    
private:
    const double stabil_;
    const double viscosity_;
    const bool   douglasWang_;
};


//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class fluid::gls::VectorLaplace
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    typedef typename base::Matrix<globalDim,globalDim>::Type        Hessian;

    //! Constructor with form and test functions 
    VectorLaplace( const double stabil,
                   const double viscosity,
                   const bool douglasWang = false )
        : stabil_(      stabil ),
          viscosity_(   viscosity ),
          douglasWang_( douglasWang )
    { }

    //--------------------------------------------------------------------------
    void tangentStiffness( const FieldTuple&  fieldTuple,
                           const LocalVecDim& xi,
                           const double       weight,
                           base::MatrixD&     matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // if pointers are identical, Galerkin-Bubnov scheme
        const bool isBubnov =
            base::auxi::EqualPointers<TestElement,TrialElement>::apply( testEp,
                                                                       trialEp );

        // Evaluate gradient of test and trial functions
        std::vector<Hessian> testHessX, trialHessX;
        const double detJ =
            (testEp -> fEFun()).evaluateHessian( geomEp, xi, testHessX );
        
        if ( isBubnov ) trialHessX = testHessX;
        else
            (trialEp -> fEFun()).evaluateHessian( geomEp, xi, trialHessX );
        
        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testHessX.size()  );
        const unsigned numColBlocks = static_cast<unsigned>( trialHessX.size() );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * globalDim );
        assert( static_cast<unsigned>( matrix.cols() ) == numColBlocks * globalDim );

        // mesh size
        const double h = base::mesh::elementSize( geomEp );

        // scalar multiplier
        const double scalar =
            (douglasWang_ ? 1.0 : -1.0) *
            stabil_ * viscosity_ * viscosity_ * h * h * detJ * weight;
        
        //
        for ( unsigned M = 0; M < numRowBlocks; M ++ ) {
            for ( unsigned N = 0; N < numColBlocks; N ++ ) {

                // dot-product of the laplacian of the shape functions
                double laplaceTest  = 0.;
                double laplaceTrial = 0.;
                
                for ( unsigned k = 0; k < globalDim; k ++ ) {
                    laplaceTrial += trialHessX[N](k,k);
                    laplaceTest  +=  testHessX[M](k,k);
                }

                // add to matrix block's diagonal entries
                for ( unsigned d = 0; d < globalDim; d++ )
                    matrix( M*globalDim + d, N*globalDim + d ) +=
                        scalar * laplaceTest * laplaceTrial;
            }
        }

        return;
    }
    
private:
    const double stabil_;
    const double viscosity_;
    const bool   douglasWang_;
};

//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class fluid::gls::PressureGradient
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    typedef typename base::Matrix<globalDim,globalDim>::Type        Hessian;

    //! Constructor with form and test functions
    PressureGradient( const double stabil,
                      const double viscosity,
                      const bool douglasWang = false )
        : stabil_(      stabil ),
          viscosity_(   viscosity ),
          douglasWang_( douglasWang )
    { }

    //--------------------------------------------------------------------------
    void tangentStiffness( const FieldTuple&  fieldTuple,
                           const LocalVecDim& xi,
                           const double       weight,
                           base::MatrixD&     matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of test and trial functions
        std::vector<Hessian> testHessX;
        (testEp -> fEFun()).evaluateHessian( geomEp, xi, testHessX );
        
        std::vector<GlobalVecDim> trialGradX;
        const double detJ =
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );
        
        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testHessX.size() );
        const unsigned numCols      = static_cast<unsigned>( trialGradX.size()  );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * globalDim );
        assert( static_cast<unsigned>( matrix.cols() ) == numCols );

        // mesh size
        const double h = base::mesh::elementSize( geomEp );

        // scalar multiplier
        const double scalar =
            (douglasWang_? -1.0 : 1.0 ) * stabil_ * h * h * viscosity_ * detJ * weight;

        // compute entries
        for ( unsigned M = 0; M < numRowBlocks; M++ ) {
            double laplaceTest = 0.;
            for ( unsigned k = 0; k < globalDim; k++ ) laplaceTest += testHessX[M](k,k);

            for ( unsigned N = 0; N < numCols; N++ ) {
                for ( unsigned i = 0; i < globalDim; i++ )
                    matrix( M*globalDim + i, N ) +=
                        scalar * laplaceTest * trialGradX[N][i];
            }
        }

        return;
    }
    
private:
    const double stabil_;
    const double viscosity_;
    const bool   douglasWang_;
};

//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class fluid::gls::VelocityDivergence
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! Sanity check
    STATIC_ASSERT_MSG( FieldTuple::numFields >= 2,
                       "Minimum number of fields violated" );

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement      GeomElement;
    typedef typename FieldTuple::TestElement      TestElement;
    typedef typename FieldTuple::TrialElement     TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    
    VelocityDivergence( const double stabil, const double viscosity )
        : stabil_( stabil ), viscosity_( viscosity )
    { }

    //! Call PressureGradient with transposed setting
    void tangentStiffness( const FieldTuple&     fieldTuple,
                           const LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        base::MatrixD aux( matrix.cols(), matrix.rows() );
        aux = base::MatrixD::Zero( matrix.cols(), matrix.rows() );

        // call pressure gradient term with transposed setting
        fluid::gls::PressureGradient<typename
                                     FieldTuple::TransposedTuple>(
                                         stabil_, viscosity_ ).tangentStiffness(
                                             fieldTuple.transpose(), xi,
                                             weight, aux );
        
        matrix += aux.transpose();
    }

private:
    const double stabil_;
    const double viscosity_;
};

//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class fluid::gls::PressureGradient2
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    typedef typename base::Matrix<globalDim,globalDim>::Type        Hessian;

    //! Constructor with form and test functions
    PressureGradient2( const double stabil,
                       const double viscosity,
                       const bool douglasWang = false )
        : stabil_(      stabil ),
          viscosity_(   viscosity ),
          douglasWang_( douglasWang )
    { }

    //--------------------------------------------------------------------------
    void tangentStiffness( const FieldTuple&  fieldTuple,
                           const LocalVecDim& xi,
                           const double       weight,
                           base::MatrixD&     matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of test and trial functions
        std::vector<Hessian> testHessX;
        (testEp -> fEFun()).evaluateHessian( geomEp, xi, testHessX );
        
        std::vector<GlobalVecDim> trialGradX;
        const double detJ =
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );
        
        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testHessX.size() );
        const unsigned numCols      = static_cast<unsigned>( trialGradX.size()  );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * globalDim );
        assert( static_cast<unsigned>( matrix.cols() ) == numCols );

        // mesh size
        const double h = base::mesh::elementSize( geomEp );

        // scalar multiplier
        const double scalar =
            (douglasWang_? -1.0 : 1.0 ) * stabil_ * h * h * viscosity_ * detJ * weight;

        // compute entries
        for ( unsigned M = 0; M < numRowBlocks; M++ ) {
            double laplaceTest = 0.;
            for ( unsigned k = 0; k < globalDim; k++ ) laplaceTest += testHessX[M](k,k);

            for ( unsigned N = 0; N < numCols; N++ ) {

                GlobalVecDim aux = base::constantVector<globalDim>( 0. );
                for ( unsigned i = 0; i < globalDim; i++ ) {
                    for ( unsigned j = 0; j < globalDim; j++ ) {
                        aux[i] += testHessX[M](i,j) * trialGradX[N][j];
                    }
                }
                
                for ( unsigned i = 0; i < globalDim; i++ )
                    matrix( M*globalDim + i, N ) +=
                        scalar * (laplaceTest * trialGradX[N][i] + aux[i]);
            }
        }

        return;
    }
    
private:
    const double stabil_;
    const double viscosity_;
    const bool   douglasWang_;
};

//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class fluid::gls::VelocityDivergence2
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! Sanity check
    STATIC_ASSERT_MSG( FieldTuple::numFields >= 2,
                       "Minimum number of fields violated" );

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement      GeomElement;
    typedef typename FieldTuple::TestElement      TestElement;
    typedef typename FieldTuple::TrialElement     TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    
    VelocityDivergence2( const double stabil, const double viscosity )
        : stabil_( stabil ), viscosity_( viscosity )
    { }

    //! Call PressureGradient with transposed setting
    void tangentStiffness( const FieldTuple&     fieldTuple,
                           const LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        base::MatrixD aux( matrix.cols(), matrix.rows() );
        aux = base::MatrixD::Zero( matrix.cols(), matrix.rows() );

        // call pressure gradient term with transposed setting
        fluid::gls::PressureGradient2<typename
                                      FieldTuple::TransposedTuple>(
                                          stabil_, viscosity_ ).tangentStiffness(
                                              fieldTuple.transpose(), xi,
                                              weight, aux );
        
        matrix += aux.transpose();
    }

private:
    const double stabil_;
    const double viscosity_;
};

//------------------------------------------------------------------------------
/** Implementation of body force times divergence of symmetric gradient of
 *  velocity test function.
 */
template<typename FIELDTUPLE>
class fluid::gls::PressureLaplace
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //! Constructor with form and test functions 
    PressureLaplace( const double stabil )
        : stabil_( stabil )
    { }

    //--------------------------------------------------------------------------
    void tangentStiffness( const FieldTuple&  fieldTuple,
                           const LocalVecDim& xi,
                           const double       weight,
                           base::MatrixD&     matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // if pointers are identical, Galerkin-Bubnov scheme
        const bool isBubnov =
            base::auxi::EqualPointers<TestElement,TrialElement>::apply( testEp,
                                                                       trialEp );

        // Evaluate gradient of test and trial functions
        std::vector<GlobalVecDim> testGradX, trialGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );
        
        if ( isBubnov ) trialGradX = testGradX;
        else
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );

        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testGradX.size()  );
        const unsigned numColBlocks = static_cast<unsigned>( trialGradX.size() );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks  );
        assert( static_cast<unsigned>( matrix.cols() ) == numColBlocks  );

        // mesh size
        const double h = base::mesh::elementSize( geomEp );
        // scalar multiplier
        const double scalar = -stabil_ * detJ * weight * h * h;
        
        //matrix += scalar * (testGradX.transpose() * formGradX);
        for ( unsigned M = 0; M < numRowBlocks; M ++ ) {
            for ( unsigned N = 0; N < numColBlocks; N ++ ) {

                double entry = 0.;
                // dot-product of function gradients
                for ( unsigned k = 0; k < globalDim; k ++ )
                    entry += testGradX[M][k] * trialGradX[N][k];
                // multiply with weight and material
                entry *= scalar;

                // add to matrix block's diagonal entries
                matrix( M, N ) += entry;
            }
        }

        return;
    }

private:
    const double stabil_;
};

//------------------------------------------------------------------------------
template<typename FIELDTUPLE, typename EVALUATEPOLICY>
class fluid::gls::BodyForce1
    : public boost::function<void( const FIELDTUPLE&,
                                   const typename
                                   base::GeomTraits<typename FIELDTUPLE::GeomElement>::
                                   LocalVecDim&,
                                   const double, base::VectorD& ) >
{
    public:
    //! @name Template parameter
    //@{
    typedef FIELDTUPLE     FieldTuple;
    typedef EVALUATEPOLICY EvaluatePolicy;
    //@}

    //! @name Extract element types
    //@{
    typedef typename FieldTuple::GeomElement GeomElement;
    typedef typename FieldTuple::TestElement TestElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Size of a DoF
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //! Type of result vector
    typedef typename base::Vector<globalDim,base::number>::Type VecDof;

    typedef typename base::Matrix<globalDim,globalDim>::Type        Hessian;

    //! How to evaluate the force function
    typedef typename EvaluatePolicy::Fun ForceFun;
    
    //--------------------------------------------------------------------------
    //! Constructor setting all references
    BodyForce1( const ForceFun& forceFun,
               const double stabil, const double viscosity, const bool douglasWang  )
        : forceFun_(    forceFun ),
          stabil_(      stabil ),
          viscosity_(   viscosity ),
          douglasWang_( douglasWang )
    { }

    //--------------------------------------------------------------------------
    /** Evaluation of the kernel function due to a body force term.
     *  \param[in]   fieldTuple Tuple of field elements
     *  \param[in]   xi         Local evaluation coordinate
     *  \param[in]   weight     Corresponding quadrature weight
     *  \param[out]  result     Result container (pre-sized and zero-initialised)
     */
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     base::VectorD&     result ) const
    {
        // extract test and trial elements from tuple
        const GeomElement* geomEp = fieldTuple.geomElementPtr();
        const TestElement* testEp = fieldTuple.testElementPtr();

        // Evaluate force function
        const typename EvaluatePolicy::result_type f
            = EvaluatePolicy::apply( geomEp, xi, forceFun_ );

        // Evaluate the shape function hessian
        std::vector<Hessian> testHessX;
        const double detJ =
            (testEp -> fEFun()).evaluateHessian( geomEp, xi, testHessX );

        // deduce the size of every contribution
        const unsigned numFun = static_cast<unsigned>( testHessX.size() );

        // mesh size
        const double h = base::mesh::elementSize( geomEp );

        const double scalar =
            (douglasWang_ ? -1.0 : 1.0) * viscosity_ * stabil_ * weight * detJ * h * h;

        // Loop over shape functions
        for ( unsigned M = 0; M < numFun; M++ ) {

            double laplaceTest = 0.;
            for ( unsigned k = 0; k < globalDim; k ++ )
                laplaceTest  +=  testHessX[M](k,k);
                    
            for ( unsigned i = 0; i < globalDim; i++ ) {
                for ( unsigned j = 0; j < globalDim; j++ ) {

                    result[ M*globalDim + i ]
                        += scalar * ( (i==j ? laplaceTest : 0. ) +
                                      testHessX[M](j,i) ) * f[j];
                }
            }
            
        }
                
        return;
        
    }
    
    //--------------------------------------------------------------------------
private:
    const ForceFun&      forceFun_;  //!< Function producing body force
    const double         stabil_;
    const double         viscosity_;
    const bool           douglasWang_;

};

//------------------------------------------------------------------------------
/** Implementation of the body force times gradient of pressure test function.
 */
template<typename FIELDTUPLE, typename EVALUATEPOLICY>
class fluid::gls::BodyForce2
    : public boost::function<void( const FIELDTUPLE&,
                                   const typename
                                   base::GeomTraits<typename FIELDTUPLE::GeomElement>::
                                   LocalVecDim&,
                                   const double, base::VectorD& ) >
{
    public:
    //! @name Template parameter
    //@{
    typedef FIELDTUPLE     FieldTuple;
    typedef EVALUATEPOLICY EvaluatePolicy;
    //@}

    //! @name Extract element types
    //@{
    typedef typename FieldTuple::GeomElement GeomElement;
    typedef typename FieldTuple::TestElement TestElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim   LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim  GlobalVecDim;

    //! Size of a DoF
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;


    //! Type of result vector
    typedef typename base::Vector<globalDim,base::number>::Type VecDof;

    typedef typename base::Matrix<globalDim,globalDim>::Type        Hessian;

    //! How to evaluate the force function
    typedef typename EvaluatePolicy::Fun ForceFun;
    
    //--------------------------------------------------------------------------
    //! Constructor setting all references
    BodyForce2( const ForceFun& forceFun, const double stabil )
        : forceFun_(    forceFun ),
          stabil_(      stabil )
    { }

    //--------------------------------------------------------------------------
    /** Evaluation of the kernel function due to a body force term.
     *  \param[in]   fieldTuple Tuple of field elements
     *  \param[in]   xi         Local evaluation coordinate
     *  \param[in]   weight     Corresponding quadrature weight
     *  \param[out]  result     Result container (pre-sized and zero-initialised)
     */
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     base::VectorD&     result ) const
    {
        // extract test and trial elements from tuple
        const GeomElement* geomEp = fieldTuple.geomElementPtr();
        const TestElement* testEp = fieldTuple.testElementPtr();

        // Evaluate force function
        const typename EvaluatePolicy::result_type f
            = EvaluatePolicy::apply( geomEp, xi, forceFun_ );

        std::vector<GlobalVecDim> testGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        // deduce the size of every contribution
        const unsigned numFun = static_cast<unsigned>( testGradX.size() );

        // mesh size
        const double h = base::mesh::elementSize( geomEp );

        const double scalar = -stabil_ * weight * detJ * h * h;

        // Loop over shape functions
        for ( unsigned M = 0; M < numFun; M++ ) {

            for ( unsigned i = 0; i < globalDim; i++ ) {
                result[ M ] += scalar * testGradX[M][i] * f[i];
            }
        }
                
        return;
    }
    
    //--------------------------------------------------------------------------
private:
    const ForceFun&      forceFun_;  //!< Function producing body force
    const double         stabil_;
};


#endif
