//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Stokes.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef fluid_stokes_hpp
#define fluid_stokes_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
#include <base/kernel/Laplace.hpp>
// fluid includes
#include <fluid/evaluations.hpp>

//------------------------------------------------------------------------------
namespace fluid{

    template<typename FIELDTUPLE> class VectorLaplace;
    template<typename FIELDTUPLE> class StressDivergence;
    template<typename FIELDTUPLE> class PressureGradient;
    template<typename FIELDTUPLE> class VelocityDivergence;
}

//------------------------------------------------------------------------------
/** Computation of the vector Laplacian of the velocity field.
 *  Stokes' system can be converted, based on the fact that the velocity vector
 *  field is solenoidal, i.e.,  \f$ \nabla \cdot u = 0 \f$, 
 *  \f[
 *       \mu \nabla \cdot (\nabla u + \nabla^T u) = \mu \nabla \cdot \nabla u
 *  \f]
 *  The corresponding bilinear form of this result is implemented in this
 *  functor.
 *  \tparam FIELDTUPLE  Type of tuple of elements for evaluation
 */
template<typename FIELDTUPLE>
class fluid::VectorLaplace
    : public base::kernel::Laplace<FIELDTUPLE>
{
public:
    typedef FIELDTUPLE FieldTuple;
    
    typedef typename base::kernel::Laplace<FieldTuple> Base;

    //! Constructor with fluid viscosity
    VectorLaplace( const double viscosity )
        : base::kernel::Laplace<FIELDTUPLE>( viscosity ),
          viscosity_( viscosity )
    { }

    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&  fieldTuple,
                        const typename Base::LocalVecDim& xi,
                        const double       weight,
                        base::VectorD&     vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple, xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Compute the residual forces due to a given velocity field.
     *  
     *  \f[
     *      F[M*d+i] = \int_\Omega \mu \phi^M_{,k} (u^{n-s}_{,k}) dx
     *  \f]
     */
    template<unsigned HIST>
    void residualForceHistory( const FieldTuple&   fieldTuple,
                               const typename Base::LocalVecDim&  xi,
                               const double        weight,
                               base::VectorD&      vector ) const
    {
        // Extract element pointer from tuple
        const typename Base::GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const typename Base::TestElement*  testEp  = fieldTuple.testElementPtr();
        const typename Base::TrialElement* trialEp = fieldTuple.trialElementPtr();
        
        // Evaluate gradient of test and trial functions
        std::vector<typename Base::GlobalVecDim> testGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        // get velocity gradient
        const typename base::Matrix<Base::globalDim,Base::nDoFs>::Type  gradU
            = fluid::velocityGradientHistory<HIST>( geomEp, trialEp, xi );
        
        for ( unsigned M = 0; M < testGradX.size(); M++ ) {
            for ( unsigned i = 0; i < Base::nDoFs; i++ ) {
            
                double dotProd = 0.;
                for ( unsigned k = 0; k < Base::globalDim; k++ )
                    dotProd += gradU( k, i ) * testGradX[M][k];

                vector[M*Base::globalDim + i] += viscosity_ * dotProd * detJ * weight;
            }
        }
    }

private:
    const double viscosity_;
};

//------------------------------------------------------------------------------
/** 
 *  \tparam FIELDTUPLE  Type of tuple of elements for evaluation
 */
template<typename FIELDTUPLE>
class fluid::StressDivergence
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
    
    //! Number of DoFs per vectorial entry
    static const unsigned nDoFs = TestElement::DegreeOfFreedom::size;

    //! Spatial dimension
    static const unsigned globalDim = GeomElement::Node::dim;
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    StressDivergence( const double viscosity )
        : viscosity_( viscosity ) { }
    
    //--------------------------------------------------------------------------
    /** 
     */
    void tangentStiffness( const FieldTuple&     fieldTuple,
                           const LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // if pointers are identical, Galerkin-Bubnov scheme
        const bool isBubnov =
            base::aux::EqualPointers<TestElement,TrialElement>::apply( testEp,
                                                                       trialEp );

        // Evaluate gradient of test functions
        std::vector<GlobalVecDim> testGradX, trialGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        if ( isBubnov ) trialGradX = testGradX;
        else
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );


        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testGradX.size() );
        const unsigned numColBlocks = static_cast<unsigned>( trialGradX.size() );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * nDoFs );
        assert( static_cast<unsigned>( matrix.cols() ) == numColBlocks * nDoFs );

        // compute entries
        for ( unsigned M = 0; M < numRowBlocks; M++ ) {
            for ( unsigned N = 0; N < numColBlocks; N++ ) {

                for ( unsigned i = 0; i < nDoFs; i++ ) {

                    // psi^M,k * phi^N_k
                    const double diag = testGradX[M].dot( trialGradX[N] );
                    
                    for ( unsigned j = 0; j < nDoFs; j++ ) {

                        matrix( M * nDoFs + i, N * nDoFs + j ) +=
                            ( testGradX[M][j] * trialGradX[N][i] +
                              (i==j ? diag : 0.) ) *
                            weight * detJ * viscosity_;
                    }
                }

            }
        }

        return;
    }

    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&  fieldTuple,
                        const LocalVecDim& xi,
                        const double       weight,
                        base::VectorD&     vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple, xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Compute the residual forces due to a given velocity field
     */
    template<unsigned HIST>
    void residualForceHistory( const FieldTuple&   fieldTuple,
                               const LocalVecDim&  xi,
                               const double        weight,
                               base::VectorD&      vector ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of test and trial functions
        std::vector<GlobalVecDim> testGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        // get velocity gradient
        const typename base::Matrix<nDoFs,nDoFs>::Type
            gradU = base::post::evaluateFieldGradientHistory<HIST>( geomEp,
                                                                    trialEp, xi );


        // assemble
        for ( unsigned M = 0; M < testGradX.size(); M++ ) {
            for ( unsigned k = 0; k < nDoFs; k++ ) {

                double entry = 0.;
                for ( unsigned i = 0; i < nDoFs; i++ )
                    entry += testGradX[M][i] * (gradU(k,i) + gradU(i,k) );

                vector[M * nDoFs + k] +=
                    entry * detJ * viscosity_ * weight;
                
            }
        }
    }

    //--------------------------------------------------------------------------
    /*
     */
    void coNormalDerivative( const FieldTuple&  fieldTuple,
                             const LocalVecDim& xi,
                             const GlobalVecDim& normal,
                             base::MatrixD& result ) const
    {
        // Extract element pointer from tuple
        const GeomElement*   geomEp  = fieldTuple.geomElementPtr();
        const TrialElement*  trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of trial functions
        std::vector<GlobalVecDim> trialGradX;
        (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );

        // number of trial functions
        const unsigned numColBlocks = static_cast<unsigned>( trialGradX.size() );

        // initiate the result with zeros
        result = base::MatrixD::Zero( +nDoFs, numColBlocks * nDoFs );

        for ( unsigned M = 0; M < trialGradX.size(); M++ ) {

            const double trialGradN = trialGradX[M].dot( normal );

            for ( unsigned i = 0; i < nDoFs; i++ ) {
                
                for ( unsigned k = 0; k < nDoFs; k++ ) {

                    result( i, M * nDoFs + k ) =
                        viscosity_ * ( trialGradX[M][i] * normal[k] +
                                       (i==k ? trialGradN : 0.) );

                }
            }
        }

        
        return;
    }

private:
    const double viscosity_;
};

//------------------------------------------------------------------------------
/** Computation of the pressure gradient term of Stokes' system.
 *  Integration by parts yields
 *  \f[
 *     \int_{\Omega} v \cdot \nabla p dx = \int_{\Gamma} v \cdot (p n) ds
 *               - \int_{\Omega} p \nabla \cdot v d x
 *  \f]
 *  The boundary term is commonly included in the natural boundary condition,
 *  i.e., \f$ t = -p n + \mu u_{,n} \f$ (or similar in the divergence form).
 *  The domain integral term is implemented in this functor.
 *  \tparam FIELDTUPLE  Type of tuple of elements for evaluation
 */
template<typename FIELDTUPLE>
class fluid::PressureGradient
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
    
    //! Number of DoFs per vectorial entry
    static const unsigned nDoFs = TestElement::DegreeOfFreedom::size;

    //! Spatial dimension
    static const unsigned globalDim = GeomElement::Node::dim;
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //--------------------------------------------------------------------------
    /** Implementation of the mixed term in Stokes' system.
     *  \f[
     *      B[ Md +i, N ] = - \int_{\Omega} \phi^M_{,i} \psi^N ds
     *  \f]
     */
    void tangentStiffness( const FieldTuple&     fieldTuple,
                           const LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of test functions
        std::vector<GlobalVecDim> testGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        // Evaluate trial functions
        typename TrialElement::FEFun::FunArray trialFun;
        (trialEp ->  fEFun()).evaluate( geomEp, xi, trialFun );

        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testGradX.size() );
        const unsigned numCols      = static_cast<unsigned>( trialFun.size()  );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * nDoFs );
        assert( static_cast<unsigned>( matrix.cols() ) == numCols );

        // compute entries
        for ( unsigned M = 0; M < numRowBlocks; M++ ) {
            for ( unsigned N = 0; N < numCols; N++ ) {

                matrix.block( M*nDoFs, N, nDoFs, 1 ) +=
                    - detJ * weight * testGradX[M] * trialFun[N];
            }
        }

        return;
    }

    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&  fieldTuple,
                        const LocalVecDim& xi,
                        const double       weight,
                        base::VectorD&     vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple, xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Compute the residual forces due to a given pressure field.
     *  
     *  \f[
     *      F[M*d+i] = - \int_\Omega \phi^M_{,i} p^{n-s} dx
     *  \f]
     */
    template<unsigned HIST>
    void residualForceHistory( const FieldTuple&   fieldTuple,
                               const LocalVecDim&  xi,
                               const double        weight,
                               base::VectorD&      vector ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of test and trial functions
        std::vector<GlobalVecDim> testGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        // get pressure
        const double p = fluid::pressureHistory<HIST>( geomEp, trialEp, xi );
        
        for ( unsigned M = 0; M < testGradX.size(); M++ ) {
            for ( unsigned i = 0; i < nDoFs; i++ ) {
            
                vector[M*globalDim + i] += -testGradX[M][i] * p * detJ * weight;
            }
        }
    }

    //--------------------------------------------------------------------------
    /** Boundary term associated with the pressure gradient term.
     *  Integration by parts gives
     *  \f[
     *      \int_\Omega \nabla p \cdot v d x =
     *      \int_\Gamma p n \cdot v d s - \int_\Omega p \nabla \cdot v d x
     *  \f]
     *  and the boundary term is represented by this function. 
     */
    void coNormalDerivative( const FieldTuple&  fieldTuple,
                             const LocalVecDim& xi,
                             const GlobalVecDim& normal,
                             base::MatrixD& result ) const
    {
        // Extract element pointer from tuple
        const GeomElement*   geomEp  = fieldTuple.geomElementPtr();
        const TrialElement*  trialEp = fieldTuple.trialElementPtr();

        typename TrialElement::FEFun::FunArray trialFun;
        (trialEp ->  fEFun()).evaluate( geomEp, xi, trialFun );

        const unsigned numCols      = static_cast<unsigned>( trialFun.size() );

        result = base::MatrixD::Zero( globalDim, numCols );

        for ( unsigned i = 0; i < numCols; i++ ) {
        
            for ( unsigned d = 0; d < globalDim; d++)
                result( d, i) = -trialFun[i] * normal[d];

        }
        
        return;
    }

};

//------------------------------------------------------------------------------
/** Computation of the velocity divergence term of Stokes' system.
 *  \f[
 *      - \int_{\Omega} q \nabla \cdot u d x
 *  \f]
 *  This term is the transposed (and possibly negative) of
 *  fluid::PressureGradient. Therefore, this functor simply delegates the call
 *  with a transposed setting.
 *  \tparam FIELDTUPLE  Type of tuple of elements for evaluation
 */
template<typename FIELDTUPLE>
class fluid::VelocityDivergence
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

    
    VelocityDivergence( const bool changeSign = false )
        : changeSign_( changeSign )
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
        fluid::PressureGradient<typename FieldTuple::TransposedTuple>().tangentStiffness(
            fieldTuple.transpose(), xi, weight, aux );
        
        if ( changeSign_ ) aux *= -1.0;
        
        matrix += aux.transpose();
    }

    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&  fieldTuple,
                        const LocalVecDim& xi,
                        const double       weight,
                        base::VectorD&     vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple, xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Compute the residual forces due to a given velocity field.
     *  
     *  \f[
     *      F[M*d+i] = - \int_\Omega \psi^M (\nabla \cdot u^{n-s}) dx
     *  \f]
     */
    template<unsigned HIST>
    void residualForceHistory( const FieldTuple&   fieldTuple,
                               const LocalVecDim&  xi,
                               const double        weight,
                               base::VectorD&      vector ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate test functions
        typename TestElement::FEFun::FunArray testFun;
        (testEp ->  fEFun()).evaluate( geomEp, xi, testFun );

        const double divU = fluid::velocityDivergenceHistory<HIST>( geomEp,
                                                                    trialEp, xi );

        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );        
        
        for ( unsigned M = 0; M < testFun.size(); M++ ) {
            vector[M] += (changeSign_ ? -1.0 : +1.0) * 
                testFun[M] * divU * detJ * weight;
        }
    }


private:
    const bool changeSign_; //!< Flag for desired sign change
};

#endif

