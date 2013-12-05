//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   NeumannForce.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_asmb_neumannforce_hpp
#define base_asmb_neumannforce_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
#include <base/aux/algorithms.hpp>
// base/asmb includes
#include <base/asmb/ForceIntegrator.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        template<typename SURFFIELDTUPLE>
        class NeumannForce;


        //----------------------------------------------------------------------
        /** Convenicen function for the computation of the terms due to applied
         *  Neumann boundary conditions.
         */
        template<typename SURFACEQUADRATURE, typename SOLVER,
                 typename BOUNDFIELD>
        void neumannForceComputation(
            const SURFACEQUADRATURE& surfaceQuadrature,
            SOLVER&                  solver, 
            const BOUNDFIELD&        boundField,
            const typename
            NeumannForce<typename BOUNDFIELD::ElementPtrTuple>::ForceFun& ff )
        {

            // object to compute the neumann force
            typedef NeumannForce<typename BOUNDFIELD::ElementPtrTuple> NeumannForce;
            typename NeumannForce::ForceFun forceFun = ff;
            NeumannForce neumannForce( forceFun );

            // integrator object
            typedef ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                    typename BOUNDFIELD::ElementPtrTuple>
                SurfaceForceInt;
            typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                boost::bind( neumannForce, _1, _2, _3, _4 );
            SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                             surfaceQuadrature, solver );

            // apply
            std::for_each( boundField.elementsBegin(),
                           boundField.elementsEnd(), surfaceForceInt );
            return;
        }
        
    } // namespace asmb
} // namespace base

//------------------------------------------------------------------------------
/** Wrapper around a surface force function.
 *  Given any forcing function \f$ f(x,n)\f$, the corresponding linear form reads
 *  \f[
 *      L(\phi) = \int_\Gamma f(x,n) \phi(x) d s_x
 *  \f]
 *  Here, the integral kernel function \f$ f(x,n) \phi(x) \f$ is represented. The
 *  computation of the integral and the assembly is done outside in
 *  base::ForceIntegrator.
 *  \tparam SURFELEMENT  Type of Surface element
 *  \tparam TESTELEMENT  Type of element of the FE test space
 */
template<typename SURFFIELDTUPLE>
class base::asmb::NeumannForce
    : public boost::function<void( const SURFFIELDTUPLE&,
                                   const typename
                                   base::GeomTraits<typename SURFFIELDTUPLE::GeomElement>::
                                   LocalVecDim&, const double, base::VectorD& )>
{
public:
    //! Template parameter
    typedef SURFFIELDTUPLE SurfFieldTuple;

    //! @name Extract element types
    //@{
    typedef typename SurfFieldTuple::GeomElement  SurfaceElement;
    typedef typename SurfFieldTuple::TestElement  TestElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim GlobalVecDim;

    //! Type of domain element
    typedef typename SurfaceElement::DomainElement    DomainElement;

    //! Type of result vector
    typedef typename
    base::VectorType<TestElement::DegreeOfFreedom::size,number>::Type VecDof;

    //! Type of force function
    typedef boost::function<VecDof( const GlobalVecDim&,
                                    const GlobalVecDim&)> ForceFun;


    //! @name Sanity checks: force function is a function of x and n
    //{
    static base::TypeEquality<typename ForceFun::first_argument_type,
                              GlobalVecDim> sanity1;
    static base::TypeEquality<typename ForceFun::first_argument_type,
                              GlobalVecDim> sanity2;
    //@}
    
    //--------------------------------------------------------------------------
    //! Constructor setting all references
    NeumannForce( const ForceFun& forceFun  )
        : forceFun_( forceFun )
    { }

    //--------------------------------------------------------------------------
    /** Evaluation of the kernel function due to a surface force term.
     *  \param[in]   surfFieldTuple  Tuple of field element pointers
     *  \param[in]   eta    Local evaluation coordinate
     *  \param[in]   weight Corresponding quadrature weight
     *  \param[out]  vector Result container (pre-sized and zero-initialised)
     */
    void operator()( const SurfFieldTuple& surfFieldTuple,
                     const LocalVecDim&    eta,
                     const double          weight,
                     base::VectorD&        vector ) const
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp = surfFieldTuple.geomElementPtr();
        const TestElement*    testEp = surfFieldTuple.testElementPtr();
        
        // Get pointer to domain element
        const DomainElement* domainEp = surfEp -> getDomainElementPointer();
        
        // Evaluate element geometry
        const GlobalVecDim x = base::Geometry<SurfaceElement>()( surfEp, eta );

        // Get surface normal and metric
        GlobalVecDim normal;
        const double detG = base::SurfaceNormal<SurfaceElement>()( surfEp, eta,
                                                                   normal );

        // Evaluate force function
        const typename ForceFun::result_type f = forceFun_( x, normal );

        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );
                
        // Evaluate the shape function
        typename TestElement::FEFun::FunArray funValues;
        (testEp -> fEFun()).evaluate( domainEp, xi, funValues );

        // deduce the size of every contribution
        const unsigned numFun = funValues.size();
        const unsigned doFSize = vector.size() / numFun;

        // Loop over shape functions
        for ( unsigned s = 0; s < numFun; s++ ) {

            vector.segment( s * doFSize,
                            doFSize )  += f * funValues[s] * weight * detG;
        }
                
        return;
        
    }
    
    //--------------------------------------------------------------------------
private:
    const ForceFun&  forceFun_;  //!< Function representing surface force
};


#endif

