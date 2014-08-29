//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   NeumannForce2.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_asmb_neumannforce2_hpp
#define base_asmb_neumannforce2_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/asmb includes
#include <base/asmb/ForceIntegrator.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        template<typename SURFFIELDTUPLE>
        class NeumannForce2;


        //----------------------------------------------------------------------
        /** Convenicen function for the computation of the terms due to applied
         *  Neumann boundary conditions.
         */
        template<typename FIELDTUPLEBINDER,
                 typename SURFACEQUADRATURE, typename SOLVER, typename FIELDBINDER>
        void neumannForceComputation2(
            const SURFACEQUADRATURE& surfaceQuadrature,
            SOLVER&                  solver, 
            const FIELDBINDER&       fieldBinder,
            const typename
            NeumannForce2<typename FIELDTUPLEBINDER::Tuple>::ForceFun& ff )
        {

            // object to compute the neumann force
            typedef NeumannForce2<typename FIELDTUPLEBINDER::Tuple> NeumannForce;
            typename NeumannForce::ForceFun forceFun = ff;
            NeumannForce neumannForce( forceFun );

            // integrator object
            typedef ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                    typename FIELDTUPLEBINDER::Tuple>
                SurfaceForceInt;
            typename SurfaceForceInt::ForceKernel surfaceForceKernel =
                boost::bind( neumannForce, _1, _2, _3, _4 );
            SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                             surfaceQuadrature, solver );

            // Apply to all elements
            typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
            typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
            for ( ; iter != end; ++iter ) {
                surfaceForceInt( FIELDTUPLEBINDER::makeTuple( *iter ) );
            }
            
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
 *  \tparam SURFFIELDTUPLE Tuple of surface elements
 */
template<typename SURFFIELDTUPLE>
class base::asmb::NeumannForce2
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

    //! Size of a DoF
    static const unsigned doFSize = TestElement::DegreeOfFreedom::size;

    //! Type of result vector
    typedef typename base::Vector<doFSize,number>::Type VecDof;

    //! Type of force function
    typedef boost::function<VecDof( const SurfaceElement*, 
                                    const LocalVecDim& )> ForceFun;


    //--------------------------------------------------------------------------
    //! Constructor setting all references
    NeumannForce2( const ForceFun& forceFun  )
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
        
        GlobalVecDim normal;
        const double detG = base::SurfaceNormal<SurfaceElement>()( surfEp, eta,
                                                                   normal );

        // Evaluate force function
        const typename ForceFun::result_type f = forceFun_( surfEp, eta );

        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );
                
        // Evaluate the shape function
        typename TestElement::FEFun::FunArray funValues;
        (testEp -> fEFun()).evaluate( domainEp, xi, funValues );

        // deduce the size of every contribution
        const unsigned numFun = static_cast<unsigned>( funValues.size() );

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

