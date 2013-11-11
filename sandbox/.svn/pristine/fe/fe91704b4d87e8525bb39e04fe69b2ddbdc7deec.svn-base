#ifndef convection_hpp
#define convection_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>
#include <base/time/derivative.hpp>

//------------------------------------------------------------------------------
template<typename FIELDTUPLE> class Convection;

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
class Convection
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
    typedef typename FieldTuple::AuxField1Element DispElement;
    typedef typename FieldTuple::AuxField2Element PressElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //! Constructor 
    Convection( const double phi, const double k, const double deltaT,
                const unsigned step = 0)
        : phi_( phi ), k_( k ), deltaT_( deltaT ), step_( step )
    { }

    void incrementTimeStep() { step_++; }

    GlobalVecDim advectionVelocity( const FieldTuple& fieldTuple,
                                    const LocalVecDim& xi ) const
    {
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const DispElement*  dispEp  = fieldTuple.auxField1ElementPtr();
        const PressElement* pressEp = fieldTuple.auxField2ElementPtr();

        // evaluate pressure gradient
        const typename base::Matrix<globalDim,1,double>::Type gradP =
            base::post::evaluateFieldGradient( geomEp, pressEp, xi );

        // evaluate time derivative of displacement
        const GlobalVecDim dotU =
            base::time::evaluateTimeDerivative( geomEp, dispEp, xi,
                                                deltaT_, step_ );

        // advection velocity
        GlobalVecDim V;
        for ( unsigned d = 0; d < globalDim; d ++ )
            V[d] = phi_ * dotU[d] - k_ * gradP( d, 0 );

        return V;
    }
    

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

        // advection velocity
        const GlobalVecDim V = this -> advectionVelocity( fieldTuple, xi );

        //
        for ( unsigned i = 0; i < numRows; i ++ ) {
            for ( unsigned j = 0; j < numCols; j ++ ) {

                double entry = 0.;
                for ( unsigned K = 0; K < globalDim; K ++ ) {

                    entry += funGradX[j][K] * V[K];
                }

                // matrix entry
                matrix( i, j ) +=
                    detJ * weight * entry * funValues[i];
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
    template<typename HIST>
    void residualForceHistory( const FieldTuple&      fieldTuple,
                               const LocalVecDim&     xi,
                               const double           weight,
                               base::VectorD&         vector ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate the test and trial functions
        typename TestElement::FEFun::FunArray testFun;
        (testEp -> fEFun()).evaluate( geomEp, xi, testFun );

        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );
        
        // Evalute concentration gradient
        const GlobalVecDim gradC
            = base::post::evaluateFieldGradientHistory<HIST>( geomEp, trialEp, xi );

        // evalute velocity field
        const GlobalVecDim V
            = this -> advectionVelocity( fieldTuple, xi ); // History !!
            
        const double convectDeriv = (gradC.transpose() * V)[0];
        
        for ( unsigned M = 0; M < testFun.size(); M++ ) {

            vector[M] += convectDeriv * testFun[M] * detJ * weight;
        }

    }

private:
    const double phi_;    //!< Porosity
    const double k_;      //!< Permeability
    const double deltaT_; //!< Time step size

    unsigned step_; //!< Number of time step
};

#endif
