//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   StressDivergence.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef fluid_stressdivergence_hpp
#define fluid_stressdivergence_hpp

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

    template<typename FIELDTUPLE> class StressDivergence;
}

//------------------------------------------------------------------------------
/** Computation of the fluid stiffness matrices refering to stress divergence.
 *  The velocity part of the stress divergence in weak form becomes
 *  \f[
 *       \int_\Omega 2 \mu \varepsilon(u) : \varepsilon(v) d x -
 *       \int_\Gamma 2 \mu \varepsilon(u) \cdot n \cdot v  d s
 *  \f]
 *   
 *
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
    /** Tangent stiffness matrix in symmetric form.
     *  From
     *  \f[
     *      \mu (\nabla u + \nabla u^T) : \nabla v =
     *      \mu (\nabla u + \nabla u^T) : \frac{1}{2}(\nabla v + \nabla v^T)
     *  \f]
     *  one gets for the entries of the stiffness matrix
     *  \f[
     *    K[Md+i, Nd+j] = \int_\Omega \mu
     *      (\psi^M_{,j} \phi^N_{,i} + \delta_{ij} \psi^M_{,k} \phi^N_{,k}) d s
     *  \f]
     *  
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
            base::auxi::EqualPointers<TestElement,TrialElement>::apply( testEp,
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

                // psi^M,k * phi^N_k
                double dotProd = 0.;
                for ( unsigned k = 0; k < nDoFs; k++ )
                    dotProd += testGradX[M][k] * trialGradX[N][k];

                for ( unsigned i = 0; i < nDoFs; i++ ) {
                    for ( unsigned j = 0; j < nDoFs; j++ ) {

                        matrix( M * nDoFs + i, N * nDoFs + j ) +=
                            ( testGradX[M][j] * trialGradX[N][i] +
                              (i==j ? dotProd : 0.) ) *
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
    /*  The co-normal derivative operator reads
     *  \f[
     *       B_n(u) = \mu (\nabla u + \nabla u^T) \cdot n
     *  \f]
     *  and its discrete counterpart is given by the matrix coefficients
     *  \f[
     *      B[i,Md+j] = \mu ( \phi^M_{,i} n_j + \delta_{ij} \phi^M_{,k} n_k )
     *  \f]
     *
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

            // normal derivative of trial function
            double dPhiDN = 0.;
            for ( unsigned k = 0; k < nDoFs; k++ )
                dPhiDN += trialGradX[M][k] * normal[k];

            for ( unsigned i = 0; i < nDoFs; i++ ) {
                for ( unsigned j = 0; j < nDoFs; j++ ) {

                    result( i, M * nDoFs + j ) =
                        viscosity_ * ( trialGradX[M][i] * normal[j] +
                                       (i==j ? dPhiDN : 0.) );

                }
            }
        }
        
        return;
    }

    //--------------------------------------------------------------------------
    /** The boundary term due to integration by parts.
     *  This term reads
     *  \f[
     *       \int_\Gamma \mu v \cdot [(\nabla u + \nabla u^T) \cdot n] d s
     *  \f]
     *  and is discretised as
     *  \f[
     *     F[Md+i] =  \int_\Gamma \mu
     *                     \psi^M ((\nabla u + \nabla u^T) \cdot n)[i] d s
     *  \f]
     *  
     *
     *  \param[in]  fieldTuple Tuple of field element pointers
     *  \param[in]  xi         Local evaluation coordinate
     *  \param[in]  normal     Surface normal \f$ n \f$
     *  \param[out] result     Result container
     */
    void boundaryResidual( const FieldTuple&   fieldTuple,
                           const LocalVecDim&  xi,
                           const GlobalVecDim& normal,
                           base::VectorD&      result ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate test functions
        typename TestElement::FEFun::FunArray testFun;
        (testEp -> fEFun()).evaluate( geomEp, xi, testFun );

        const std::size_t numRowBlocks = testFun.size();

        // initialise result container
        result = base::VectorD::Zero( numRowBlocks * nDoFs );

        const typename base::Matrix<nDoFs,nDoFs>::Type gradU =
            base::post::evaluateFieldGradient( geomEp, trialEp, xi );

        typename base::Vector<nDoFs>::Type dUDN;
        dUDN.noalias() = (gradU + gradU.transpose()) * normal;
        
        //
        for ( std::size_t M = 0; M < numRowBlocks; M++ ) {
            for ( unsigned i = 0; i < nDoFs; i++ ) {

                result[M*nDoFs+i] = viscosity_ * dUDN[i] * testFun[M];
            }
        }

        return;
    }

private:
    const double viscosity_;
};

#endif
