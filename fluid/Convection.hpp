//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   fluid/Convection.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef fluid_convection_hpp
#define fluid_convection_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// fluid includes
#include <fluid/evaluations.hpp>

//------------------------------------------------------------------------------
namespace fluid{

    template<typename FIELDTUPLE>
    class Convection;
}

//#define NEWTON
//------------------------------------------------------------------------------
/** System matrix contribution due to convection.
 *  The convection term in conservative form reads
 *  \f[
 *        \int_\Omega \rho v \cdot ( w \cdot \nabla u)
 *               + \rho \frac{1}{2} (\nabla \cdot w) v \cdot u dx
 *  \f]
 *  Here, \f$ w \f$ is the advection velocity which in this implementation is
 *  set to the current (latest computational result) velocity field.
 *  \tparam FIELDTUPLE  Type of tuple of elements for evaluation
 */
template<typename FIELDTUPLE>
class fluid::Convection
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! Sanity check
    STATIC_ASSERT_MSG( FieldTuple::numFields >= 3,
                       "Minimum number of fields violated" );


    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement      GeomElement;
    typedef typename FieldTuple::TestElement      TestElement;
    typedef typename FieldTuple::TrialElement     TrialElement;
    typedef typename FieldTuple::AuxField1Element VelocElement;
    //@}


    //! Number of DoFs per vectorial entry
    static const unsigned nDoFs = TestElement::DegreeOfFreedom::size;

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Constructor with fluid density
    Convection( const double density )
    : density_( density )
    { }

    //--------------------------------------------------------------------------
    /** Implementation of the convection term in conservative form.
     *  \f[
     *      C[Md+i,Nd+j] = \int_{\Omega} \rho \phi^M \delta_{ij}
     *          (\sum_k u_k \phi^N_{,k} + \frac{1}{2} \nabla \cdot u \phi^N) dx
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
        const VelocElement* velocEp = fieldTuple.auxField1ElementPtr();

        // Evaluate the velocity field
        const typename base::Vector<nDoFs,double>::Type uAdv =
            fluid::velocity( geomEp, velocEp, xi );

        // Evaluate its divergence
        const double divUAdv = fluid::velocityDivergence( geomEp, trialEp, xi );

#ifdef NEWTON
        // Evaluate its gradient
        const typename base::Matrix<nDoFs,nDoFs>::Type gradU =
            base::post::evaluateFieldGradient( geomEp, trialEp, xi );
#endif
        
        // Evaluate gradient of trial functions 
        std::vector<GlobalVecDim> trialGradX;
        const double detJ =
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );

        // Evaluate the test and trial functions
        typename TestElement::FEFun::FunArray testFun;
        (testEp -> fEFun()).evaluate( geomEp, xi, testFun );
        typename TrialElement::FEFun::FunArray trialFun;
        (trialEp -> fEFun()).evaluate( geomEp, xi, trialFun );

        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testFun.size()    );
        const unsigned numColBlocks = static_cast<unsigned>( trialGradX.size() );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * nDoFs );
        assert( static_cast<unsigned>( matrix.cols() ) == numColBlocks * nDoFs );

        // loop over the test and trial functions
        for ( unsigned M = 0; M < numRowBlocks; M++ ) { // test functions
            for ( unsigned N = 0; N < numColBlocks; N++ ) { // trial functions

                // sum_k  w_k phi^N_{,k}
                double advTrial = 0.;
                for ( unsigned k = 0; k < nDoFs; k++ )
                    advTrial += uAdv[k] * trialGradX[N][k];

                const double entry =
                    testFun[M] * (advTrial
                                  + 0.5 * divUAdv * trialFun[N]
                        ) * density_ * detJ * weight;

                for ( unsigned k = 0; k < nDoFs; k++ )
                    matrix( M*nDoFs+k, N*nDoFs+k ) += entry;

#ifdef NEWTON
                for ( unsigned i = 0; i < nDoFs; i++ ) {
                    for ( unsigned j = 0; j < nDoFs; j++ ) {
                
                        const double entry2 =
                            testFun[M] *
                            ( trialFun[N] * gradU(j,i)
                              + 0.5 * trialGradX[N][j] * uAdv[i] ) *
                            density_ * detJ * weight;
                        
                        matrix( M*nDoFs+i, N*nDoFs+j ) += entry2;
                    }
                }
#endif 
                
            } // N
        } // M
        return;
    }

    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&     fieldTuple,
                        const LocalVecDim&    xi,
                        const double          weight,
                        base::VectorD&        vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple, xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Contribution of the convection term to residual forces.
     *  \f[
     *      F[Md+i] = \int_\Omega \rho \phi^M (u \cdot \nabla u_i^{n-s}) dx
     *  \f]
     */
    template<unsigned HIST>
    void residualForceHistory( const FieldTuple&      fieldTuple,
                               const LocalVecDim&     xi,
                               const double           weight,
                               base::VectorD&         vector ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();
        const VelocElement* velocEp = fieldTuple.auxField1ElementPtr();

        // Evaluate the test and trial functions
        typename TestElement::FEFun::FunArray testFun;
        (testEp -> fEFun()).evaluate( geomEp, xi, testFun );

        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );
        
        // Evalute velocity gradient
        const typename base::Matrix<GeomElement::Node::dim,
                                    TrialElement::DegreeOfFreedom::size,
                                    double>::Type
            gradU = fluid::velocityGradientHistory<HIST>( geomEp, velocEp, xi );

        // Evalute velocity field
        const typename base::Vector<TrialElement::DegreeOfFreedom::size,
                                    double>::Type
            U = fluid::velocityHistory<HIST>( geomEp, trialEp, xi );

        for ( unsigned i = 0; i < nDoFs; i++ ) {
            // (u \nabla) u
            double convectiveDeriv = 0.;
            for ( unsigned k = 0; k < nDoFs; k++ )
                convectiveDeriv += U[k] * gradU(k,i);
            
            for ( unsigned M = 0; M < testFun.size(); M++ ) {
                
                vector[M*nDoFs + i] += density_ * convectiveDeriv * testFun[M]
                    * detJ * weight;
            }
        }

        return;
    }

    
private:
    const double density_; //!< Fluid density
};

#endif
