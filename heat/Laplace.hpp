//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   heat/Laplace.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef heat_laplace_hpp
#define heat_laplace_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
#include <base/types.hpp>
#include <base/aux/functions.hpp>
#include <base/kernel/Laplace.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>
// head includes
#include <heat/Temperature.hpp>

//------------------------------------------------------------------------------
namespace heat{

    template<typename FIELDTUPLE>
    class Laplace;

}

//------------------------------------------------------------------------------
/** Representation of the bilinear form for the Laplace operator.
 *  This bilinear form reads
 *  \f[
 *      a(u,v) = \int_\Omega \kappa \nabla u \cdot \nabla v d x
 *  \f]
 *  \tparam FIELDTUPLE  Type of tuple of elements for evaluation
 */
template<typename FIELDTUPLE>
class heat::Laplace
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! Sanity check
    STATIC_ASSERT_MSG( FieldTuple::numFields >= 2,
                       "Minimum number of fields violated" );

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
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

    //! Type of material function
    typedef boost::function< double( const GeomElement *,
                                     const LocalVecDim & ) > ConductivityFun;

    //! Constructor with form and test functions 
    Laplace( const double materialFactor )
        : conductivityFun_( base::aux::ConstantFun<double>( materialFactor ) )
    { }

    //! Set other material functions
    void setConductivityFunction( ConductivityFun& cf )
    {
        conductivityFun_ = cf;
    }


    //--------------------------------------------------------------------------
    /**  Contribution to the element stiffness matrix in a quadrature rule
     *   The element stiffness matrix for the Laplace operator reads
     *   \f[
     *       K[m, n] = \int_\Omega \mu(x) \nabla phi^n \cdot \nabla phi^m dx
     *   \f]
     *   This object adds the weighted integrand evaluated at a local coordinate
     *   \f$\xi\f$ to a provided storage.
     *   \param[in]  fieldTuple Tuple of field element pointers
     *   \param[in]  xi         Local coordinate: quadrature point
     *   \param[in]  weight     Weight corresponding to the quadrature point
     *   \param[out] matrix     Result storage
     */
    void tangentStiffness( const FieldTuple&  fieldTuple, 
                           const LocalVecDim& xi,
                           const double       weight,
                           base::MatrixD&     matrix ) const
    {
        // Evaluate conductivity at local coordinate of element
        const double conductivity = conductivityFun_( fieldTuple.geomElementPtr(),
                                                      xi );

        // call generic laplace kernel
        base::kernel::Laplace<FieldTuple> laplace( conductivity );
        laplace.tangentStiffness( fieldTuple, xi, weight, matrix );

    }
    
    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&  fieldTuple,
                        const LocalVecDim& xi,
                        const double       weight,
                        base::VectorD&     vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple,
                                         xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Compute the residual forces due to a given temperature field.
     *  
     *  \f[
     *      F[i] = \int_\Omega \nabla \phi^i \kappa \nabla u^{n-s} d x
     *  \f]
     *  \tparam   HIST Number of history term to use, \f$ s \f$
     *  \param[in] fieldTuple   Tuple of element pointers
     *  \param[in] xi           Evaluation coordinate
     *  \param[in] weight       Quadrature weight
     *  \param[in] vector       Result container
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
        
        // Evaluate conductivity at local coordinate of element
        const double conductivity = conductivityFun_( geomEp, xi );

        // Evaluate gradient of test and trial functions
        std::vector<GlobalVecDim> testGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        // Evalute temperature gradient 
        const mat::Vector gradU
            = heat::temperatureGradientHistory<HIST>( geomEp, trialEp, xi );

        // times conductivity
        const GlobalVecDim flux = conductivity * gradU.head( globalDim );

        for ( unsigned i = 0; i < testGradX.size(); i++ ) {

            const double dotProd = base::dotProduct( flux, testGradX[i] );

            vector[i] += dotProd * detJ * weight;
        }
        
    }
    
private:
    //! Functor representing the material conductivity
    ConductivityFun conductivityFun_;
};

#endif
