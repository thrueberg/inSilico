//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   kernel/Mass.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_kernel_mass_hpp
#define base_kernel_mass_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/aux includes
#include <base/aux/functions.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace kernel{

        template<typename FIELDTUPLE>
        class Mass;
    }
}

//------------------------------------------------------------------------------
/** Representation of the bilinear form for the scalar product.
 *  This term, usually arising from time discretisations, has the form
 *  \f[
 *      m(u,v) = \int_\Omega \alpha u \cdot v dx
 *  \f]
 *  where \f$ \alpha \f$ is some scalar multiplier, e.g.
 *  \f$ \alpha = \frac{\rho}{\Delta t} \f$ with material density and time step
 *  size.
 *
 *  \tparam FIELDTUPLE  Tuple of field elements
 */
template<typename FIELDTUPLE>
class base::kernel::Mass
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    typedef void result_type;

    //! Flag for equal test and form functions --> Bubnov-Galerkin
    static const bool bubnov = boost::is_same<TrialElement,
                                              TestElement>::value;
    
    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Size of DOFs
    static const unsigned doFSize   = TrialElement::DegreeOfFreedom::size;

    //! Constructor with form and test functions 
    Mass( const double factor )
    : factor_( factor )
    { }

    //--------------------------------------------------------------------------
    /**  Contribution to the element stiffness matrix in a quadrature rule
     *   The element stiffness matrix for the \f$ L_2 \f$ product reads
     *   \f[
     *       M[M*d+i,N*d+i] = \int_\Omega \alpha \phi^M \phi^N dx
     *   \f]
     *   This object adds the weighted integrand evaluated at a local coordinate
     *   \f$\xi\f$ to a provided storage.
     *   \param[in]  fieldTuple Tuple of field element pointers
     *   \param[in]  xi       Local coordinate: quadrature point
     *   \param[in]  weight   Weight corresponding to the quadrature point
     *   \param[out] matrix   Result storage
     */
    void operator()( const FieldTuple&  fieldTuple,
                     const LocalVecDim& xi,
                     const double       weight,
                     base::MatrixD&     matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate trial functions
        typename TrialElement::FEFun::FunArray trialFun;
        (trialEp ->  fEFun()).evaluate( geomEp, xi, trialFun );

        // Evaluate test functions
        typename TestElement::FEFun::FunArray  testFun;
        if ( bubnov ) testFun = trialFun;
        else
            (testEp -> fEFun()).evaluate( geomEp, xi, testFun );
        
        // element jacobian
        const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );
        
        // Sizes and sanity checks
        const unsigned numRowBlocks =  testFun.size();
        const unsigned numColBlocks = trialFun.size();
        assert( static_cast<unsigned>( matrix.rows() ) * doFSize == numRowBlocks );
        assert( static_cast<unsigned>( matrix.cols() ) * doFSize == numColBlocks );

        // scalar multiplier
        const double scalar = factor_ * detJ * weight;
        
        //matrix += scalar * (testGradX.transpose() * formGradX);
        for ( unsigned M = 0; M < numRowBlocks; M++ ) {
            for ( unsigned N = 0; N < numColBlocks; N++ ) {

                const double entry = scalar * testFun[M] * trialFun[N];

                // add to matrix block's diagonal entries
                for ( unsigned d = 0; d < doFSize; d++ )
                    matrix( M*doFSize + d, N*doFSize + d ) += entry;
            }
        }
        
    }

private:
    const double factor_; //!< Scalar multiplier (e.g. density over step size)
};

#endif
