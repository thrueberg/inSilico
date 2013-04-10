//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Static.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef heat_static_hpp
#define heat_static_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// heat includes
#include <heat/Temperature.hpp>

//------------------------------------------------------------------------------
namespace heat{

    template<typename MATERIAL,
             typename FIELDTUPLE>
    class Static;

}

//------------------------------------------------------------------------------
/** Static heat distribution with a user-defined material behaviour.
 *  Let the conductivity \f$ \kappa(u) \f$ be a given tensor-valued function of
 *  the current temperature \f$ u \f$.  Then a step in a Newton method reads
 *  \f[
 *      \int_{\Omega} \nabla v \cdot \left(
 *         \kappa(u^n) \cdot \nabla (\Delta u) +
 *         \kappa^{\prime}(u^n) \cdot \nabla u^n \Delta u \right) dx =
 *    - \int_{\Omega} \nabla v \cdot \kappa(u^n) \cdot \nabla u^n dx + f^{ext}
 *  \f]
 *  for the unknown temperature increment \f$ \Delta u \f$ and all evaluations
 *  at the given temperature state \f$ u^n \f$.
 *  This class implements the left hand side terms in tangentStiffness() and the
 *  right hand side term in residualForce(). Note that the resulting element
 *  matrix is in general not symmetric, unless the derivative
 *  \f$ \kappa^{\prime} \f$ is zero.
 *  \tparam MATERIAL   Type of material providing the conductivity and gradient
 *  \tparam FIELDTUPLE Tuple of field element pointers
 */
template<typename MATERIAL, typename FIELDTUPLE>
class heat::Static
{
public:
    //! @name Template parameters
    //@{
    typedef MATERIAL     Material;
    typedef FIELDTUPLE   FieldTuple;
    //@}

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement      GeomElement;
    typedef typename FieldTuple::TestElement      TestElement;
    typedef typename FieldTuple::TrialElement     TrialElement;
    //@}


    //! Flag for equal test and form functions --> Bubnov-Galerkin
    static const bool bubnov = boost::is_same<TrialElement,
                                              TestElement>::value;
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //! Constructor with material object
    Static( const Material& material )
        : material_( material )
    { }

    //--------------------------------------------------------------------------
    /**  Contribution to the element stiffness matrix in a quadrature rule.
     *   The element stiffness matrix for the static heat transfer reads
     *   \f[
     *     K[i, j] = \int_\Omega \nabla \phi^i \dot \kappa \cdot \nabla \phi^j
     *               + \nabla \phi^i \kappa^\prime \nabla u \phi^j dx
     *   \f]
     *   \param[in]  fieldTuple Tuple of field element pointers
     *   \param[in]  xi       Local coordinate: quadrature point
     *   \param[in]  weight   Weight corresponding to the quadrature point
     *   \param[out] matrix   Result storage
     */
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
        std::vector<GlobalVecDim> testGradX, trialGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        if ( bubnov ) trialGradX = testGradX;
        else
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );

        // evaluate trial functions
        typename TrialElement::FEFun::FunArray funValues;
        (trialEp -> fEFun()).evaluate( geomEp, xi, funValues );

        // Sizes and sanity checks
        const unsigned numRows = testGradX.size();
        const unsigned numCols = trialGradX.size();
        assert( static_cast<unsigned>( matrix.rows() ) == numRows );
        assert( static_cast<unsigned>( matrix.cols() ) == numCols );

        // evaluate temperature and its gradient
        const double      u     = heat::temperature(         geomEp, trialEp, xi );
        const mat::Vector gradU = heat::temperatureGradient( geomEp, trialEp, xi );

        // get material behaviour
        mat::Tensor conductivity, conductivityGradient;
        material_.conductivity(          u, gradU, conductivity         );
        material_.conductivityGradient(  u, gradU, conductivityGradient );

        //
        for ( unsigned i = 0; i < numRows; i ++ ) {
            for ( unsigned j = 0; j < numCols; j ++ ) {

                double entry = 0.;
                for ( unsigned K = 0; K < globalDim; K ++ ) {
                    for ( unsigned L = 0; L < globalDim; L++ ) {

                        // contribution of conductivity
                        entry +=
                            testGradX[i][K] * conductivity(K,L) * trialGradX[j][L];

                        
                        // contribution of conductivity derivative
                        entry +=
                            testGradX[i][K] * conductivityGradient(K,L) *
                            gradU[L] * funValues[j];
                    }
                }

                // matrix entry
                matrix( i, j ) += detJ * weight * entry;
            }
        }
        return;
    }

    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&   fieldTuple,
                        const LocalVecDim&  xi,
                        const double        weight,
                        base::VectorD&      vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple, xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Compute the internal forces given the current state of temperature.
     *  \f[
     *     F^{int}[K] = -\int_{\Omega}
     *             \nabla \phi^K \cdot \kappa(u^n) \cdot \nabla u^n dx
     *  \f]
     *   \param[in]  fieldTuple Tuple of field element pointers
     *   \param[in]  xi       Local coordinate: quadrature point
     *   \param[in]  weight   Weight corresponding to the quadrature point
     *   \param[out] vector   Result storage
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

        const unsigned numRows = testGradX.size();
        assert( static_cast<unsigned>( vector.size() ) == numRows );

        // evaluate temperature and its gradient
        const double      u     = heat::temperatureHistory<HIST>( geomEp, trialEp, xi );
        const mat::Vector gradU = heat::temperatureGradientHistory<HIST>( geomEp,
                                                                          trialEp, xi );

        // get material behaviour
        mat::Tensor conductivity;
        material_.conductivity( u, gradU, conductivity );
        
        // loop over the test functions
        for ( unsigned i = 0; i < numRows; i++ ) { // test functions

            double sum = 0.;
            for ( unsigned K = 0; K < globalDim; K++ ) {
                for ( unsigned L = 0; L < globalDim; L++ ) {
                    sum += testGradX[i][K] * conductivity(K,L) * gradU[L];
                }
            }

            vector( i ) += sum * detJ * weight;
            
        }
        return;
    }


private:
    const Material& material_; //!< Material object
};

#endif
