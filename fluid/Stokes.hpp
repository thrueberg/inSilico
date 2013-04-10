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

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //! Flag for equal test and form functions --> Bubnov-Galerkin
    static const bool bubnov = boost::is_same<TrialElement,
                                              TestElement>::value;
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Constructor with fluid viscosity
    VectorLaplace( const double viscosity )
    : viscosity_( viscosity )
    { }

    //--------------------------------------------------------------------------
    /** Implementation of the vector Laplacian.
     *  \f[
     *      K[Md +i, Nd +j] = \int_\Omega \mu \phi^M_{,i} \phi^N_{,j}
     *                                   \delta_{ij} dx
     *  \f]
     *  
     */
    void tangentStiffness( const FieldTuple&     fieldTuple,
                           const LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        // call generic laplace kernel
        base::kernel::Laplace<FieldTuple> laplace( viscosity_ );
        laplace( fieldTuple, xi, weight, matrix );
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
     *      F[M*d+i] = \int_\Omega \mu \phi^M_{,k} (u^{n-s}_{,k}) dx
     *  \f]
     */
    template<typename HIST>
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
        const typename base::MatrixType<globalDim, nDoFs>::Type  gradU
            = fluid::velocityGradientHistory<HIST>( geomEp, trialEp, xi );
        
        for ( unsigned M = 0; M < testGradX.size(); M++ ) {
            for ( unsigned i = 0; i < nDoFs; i++ ) {
            
                double dotProd = 0.;
                for ( unsigned k = 0; k < globalDim; k++ )
                    dotProd += gradU( k, i ) * testGradX[M][k];

                vector[M*globalDim + i] += viscosity_ * dotProd * detJ * weight;
            }
        }
    }

private:
    const double viscosity_; //!< Fluid viscosity
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
        const unsigned numRowBlocks = testGradX.size();
        const unsigned numCols      = trialFun.size();
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
    template<typename HIST>
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
    template<typename HIST>
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
            vector[M] += (changeSign_ ? +1.0: -1.0) * 
                testFun[M] * divU * detJ * weight;
        }
    }


private:
    const bool changeSign_; //!< Flag for desired sign change
};

#endif

