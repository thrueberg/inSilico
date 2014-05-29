//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   VelocityDivergence.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef fluid_velocitydivergence_hpp
#define fluid_velocitydivergence_hpp

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

    template<typename FIELDTUPLE> class VelocityDivergence;
}


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

