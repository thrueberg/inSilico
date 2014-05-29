#ifndef convection_hpp
#define convection_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// fluid includes
#include <fluid/evaluations.hpp>

//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class Convection
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement      GeomElement;
    typedef typename FieldTuple::TestElement      TestElement;
    typedef typename FieldTuple::TrialElement     TrialElement;
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

        // Evaluate the velocity field
        const typename base::Vector<nDoFs,double>::Type uAdv =
            fluid::velocityHistory<1>( geomEp, trialEp, xi );

        // Evaluate its divergence
        const double divUAdv =
            fluid::velocityDivergenceHistory<1>( geomEp, trialEp, xi );

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

            } // N
        } // M
        return;
    }

    
private:
    const double density_; //!< Fluid density
};

#endif
