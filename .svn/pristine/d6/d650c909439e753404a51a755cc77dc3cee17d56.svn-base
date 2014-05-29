//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   heat/Convection.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef heat_convection_hpp
#define heat_convection_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// head includes
#include <heat/Velocity.hpp>

//------------------------------------------------------------------------------
namespace heat{

    template<typename FIELDTUPLE>
    class Convection;

}

//------------------------------------------------------------------------------
/** Convection term of heat transport equation.
 *  Implements the bilinear form (linear in \f$ u \f$ and \f$ v \f$)
 *  \f[
 *       c( V; u, v) = \int_\Omega \rho V \cdot (\nabla u ) v d x
 *  \f]
 *  with a given vector field \f$ V \f$ which is provided via a third field
 *  in the FIELDTUPLE.
 *  \tparam FIELDTUPLE  Type of tuple of elements for evaluation
 */
template<typename FIELDTUPLE>
class heat::Convection
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

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //! Constructor with density
    Convection( const double density )
        : density_( density )
    { }

    //--------------------------------------------------------------------------
    /**  Contribution to the element stiffness matrix in a quadrature rule.
     *   The element stiffness matrix for the static heat transfer reads
     *   \f[
     *          C[i,j] = \int_\Omega (V \cdot \nabla \phi^j) \phi^i dx
     *   \f]
     *   \param[in]  fieldTuple Tuple of field element pointers
     *   \param[in]  xi       Local coordinate: quadrature point
     *   \param[in]  weight   Weight corresponding to the quadrature point
     *   \param[out] matrix   Result storage
     */
    void tangentStiffness( const FieldTuple&      fieldTuple,
                           const LocalVecDim&     xi,
                           const double           weight,
                           base::MatrixD&         matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();
        const VelocElement* velocEp = fieldTuple.auxField1ElementPtr();

        // Evaluate gradient of trial functions
        std::vector<GlobalVecDim> funGradX;
        const double detJ =
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, funGradX );

        // evaluate test functions
        typename TestElement::FEFun::FunArray funValues;
        (testEp -> fEFun()).evaluate( geomEp, xi, funValues );

        // Sizes and sanity checks
        const unsigned numRows = static_cast<unsigned>( funValues.size() );
        const unsigned numCols = static_cast<unsigned>( funGradX.size()  );
        assert( static_cast<unsigned>( matrix.rows() ) == numRows );
        assert( static_cast<unsigned>( matrix.cols() ) == numCols );

        // evalute velocity field
        const mat::Vector V = heat::velocity( geomEp, velocEp, xi );

        //
        for ( unsigned i = 0; i < numRows; i ++ ) {
            for ( unsigned j = 0; j < numCols; j ++ ) {

                double entry = 0.;
                for ( unsigned K = 0; K < globalDim; K ++ ) {

                    entry += funGradX[j][K] * V[K];
                }

                // matrix entry
                matrix( i, j ) +=
                    density_ * detJ * weight * entry * funValues[i];
            }
        }
        return;
    }

    //--------------------------------------------------------------------------
    //! Residual force computation, delegate to residualForceHistory<0>()
    void residualForce( const FieldTuple&      fieldTuple,
                        const LocalVecDim&     xi,
                        const double           weight,
                        base::VectorD&         vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple, 
                                         xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    /** Contribution of the convection term to residual forces.
     *  \f[
     *      F[M] = \int_\Omega \phi^M ( v \cdot \nabla u^{n-s} ) dx
     *  \f]
     *  \tparam   HIST Number of history term to use, \f$ s \f$
     *  \param[in] fieldTuple   Tuple of element pointers
     *  \param[in] xi           Evaluation coordinate
     *  \param[in] weight       Quadrature weight
     *  \param[in] vector       Result container
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
        
        // Evalute temperature gradient 
        const mat::Vector gradU
            = heat::temperatureGradientHistory<HIST>( geomEp, trialEp, xi );

        // evalute velocity field
        const mat::Vector V = heat::velocityHistory<HIST>( geomEp, velocEp, xi );

        const double convectDeriv = (gradU.transpose() * V)[0];
        
        for ( unsigned M = 0; M < testFun.size(); M++ ) {

            vector[M] += convectDeriv * testFun[M] * detJ * weight;
        }

    }

private:
    const double density_; //!< Mass density
};

#endif
