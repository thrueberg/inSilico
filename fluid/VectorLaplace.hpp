//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   VectorLaplace.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef fluid_vectorlaplace_hpp
#define fluid_vectorlaplace_hpp

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
    : public base::kernel::Laplace<FIELDTUPLE>
{
public:
    typedef FIELDTUPLE FieldTuple;
    
    typedef typename base::kernel::Laplace<FieldTuple> Base;
    typedef typename Base::LocalVecDim                 LocalVecDim;
    typedef typename Base::GlobalVecDim                GlobalVecDim;

    typedef typename Base::GeomElement                 GeomElement;
    typedef typename Base::TestElement                 TestElement;
    typedef typename Base::TrialElement                TrialElement;
    
    //! Constructor with fluid viscosity
    VectorLaplace( const double viscosity )
        : base::kernel::Laplace<FIELDTUPLE>( viscosity ),
          viscosity_( viscosity )
    { }

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
        const typename base::Matrix<Base::globalDim,Base::doFSize>::Type  gradU
            = fluid::velocityGradientHistory<HIST>( geomEp, trialEp, xi );
        
        for ( unsigned M = 0; M < testGradX.size(); M++ ) {
            for ( unsigned i = 0; i < Base::doFSize; i++ ) {
            
                double dotProd = 0.;
                for ( unsigned k = 0; k < Base::globalDim; k++ )
                    dotProd += gradU( k, i ) * testGradX[M][k];

                vector[M*Base::doFSize + i] += viscosity_ * dotProd * detJ * weight;
            }
        }
    }

    
private:
    const double viscosity_;
};

#endif
