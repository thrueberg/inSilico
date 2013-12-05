//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ErrorNorm.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_post_errornorm_hpp
#define base_post_errornorm_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/linearAlgebra.hpp>
#include <base/geometry.hpp>
// base/asmb includes
#include <base/asmb/FieldElementPointerTuple.hpp>
// base/post includes
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace post{

        template<typename GEOMELEMENT, typename FIELDELEMENT,unsigned ORDER>
        class ErrorNorm;

        //----------------------------------------------------------------------
        /** Convenience function for the computation of the FE error in a norm.
         *  \tparam ORDER      Order of Sobolev norm
         *  \tparam QUADRATURE Type of quadrature to apply
         *  \tparam MESH       Type of geometry mesh
         *  \tparam FIELD      Type of field with the FE solution
         *  \param[in]   quadrature  Quadrature object
         *  \param[in]   mesh        Mesh object
         *  \param[in]   field       Field object with FE solution
         *  \param[in]   refSol      Function object providing reference solution
         *  \return                  Error in desired Sobolev norm
         */
        template<unsigned ORDER, typename QUADRATURE, typename MESH,
                 typename FIELD>
        double errorComputation( const QUADRATURE& quadrature,
                                 const MESH&       mesh,
                                 const FIELD&      field,
                                 const typename
                                 base::post::ErrorNorm<
                                     typename MESH::Element,
                                     typename FIELD::Element,
                                     ORDER>::Reference& refSol )
        {
            typedef ErrorNorm<typename MESH::Element,typename FIELD::Element,
                              ORDER> Error;
            Error error( refSol );

            // compute for every element the error
            typename MESH::ElementPtrConstIter  elemIter  = mesh.elementsBegin();
            typename MESH::ElementPtrConstIter  elemLast  = mesh.elementsEnd();

            double errorSquared = 0.;
    
            for ( ; elemIter != elemLast; ++elemIter ) {

                // get corresponding field element
                typename FIELD::Element* fieldElem =
                    field.elementPtr( (*elemIter) -> getID() );
                    

                // Construct a field element pointer tuple 
                base::asmb::FieldElementPointerTuple<
                    typename MESH::Element*,
                    typename FIELD::Element*> fept( *elemIter, fieldElem );

                // apply quadrature
                quadrature.apply( error, fept, errorSquared );
            }
            
            // return the square root of the sum of the squares
            return  std::sqrt( errorSquared );
        }
        
        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            //! Return type of approximate & reference solution evaluation calls
            template<unsigned DIM, unsigned SIZE, unsigned ORDER>
            struct ValueBinder;

            template<unsigned DIM, unsigned SIZE>
            struct ValueBinder<DIM,SIZE,0>
            {
                typedef typename base::Vector<SIZE,number>::Type Type;
            };

            template<unsigned DIM, unsigned SIZE>
            struct ValueBinder<DIM,SIZE,1>
            {
                typedef typename base::Matrix<DIM,SIZE,number>::Type Type;
            };

            //------------------------------------------------------------------
            //! Redirect calls to evaluation of field or the gradient thereof
            template<unsigned DIM, unsigned SIZE,unsigned ORDER>
            struct FieldEvaluationBinder;

            template<unsigned DIM, unsigned SIZE>
            struct FieldEvaluationBinder<DIM,SIZE,0>
            {
                template<typename GEOMELEMENT, typename FIELDELEMENT>
                static typename ValueBinder<DIM,SIZE,0>::Type
                apply( const GEOMELEMENT* gep, const FIELDELEMENT* fep,
                       const typename FIELDELEMENT::FEFun::VecDim& xi )
                {
                    return base::post::evaluateField( gep, fep, xi );
                }
            };
        
            template<unsigned DIM, unsigned SIZE>
            struct FieldEvaluationBinder<DIM,SIZE,1>
            {
                template<typename GEOMELEMENT, typename FIELDELEMENT>
                static typename ValueBinder<DIM,SIZE,1>::Type
                apply( const GEOMELEMENT* gep, const FIELDELEMENT* fep,
                       const typename FIELDELEMENT::FEFun::VecDim& xi )
                {
                    return base::post::evaluateFieldGradient( gep, fep, xi );
                }
            };

            //------------------------------------------------------------------
            //! Case FIELDDIM != GEOMDIM, convert from surface to domain
            template<typename GEOMELEMENT, unsigned FIELDDIM, unsigned GEOMDIM>
            struct SurfacePolicy
            {
                typedef const typename GEOMELEMENT::DomainElement DomainElement;

                static DomainElement* domainElementPtr( const GEOMELEMENT* ge)
                {
                    return ge -> getDomainElementPointer();
                }

                typedef typename base::Vector<FIELDDIM>::Type VecDim;
                typedef typename base::Vector<GEOMDIM>::Type  VecLDim;

                static VecDim makeParametric( const GEOMELEMENT* ge,
                                              const VecLDim& eta )
                {
                    return ge -> localDomainCoordinate( eta );
                }

            };

            //! Case FIELDDIM == GEOMDIM, no conversion needed
            template<typename GEOMELEMENT, unsigned FIELDDIM>
            struct SurfacePolicy<GEOMELEMENT,FIELDDIM,FIELDDIM>
            {
                typedef const GEOMELEMENT DomainElement;

                static DomainElement* domainElementPtr( const GEOMELEMENT* ge)
                {
                    return ge;
                }

                typedef typename base::Vector<FIELDDIM>::Type VecDim;
                
                static VecDim makeParametric( const GEOMELEMENT* ge,
                                              const VecDim& eta )
                {
                    return eta;
                }

            };

        
        } // namespace detail_
    } // namespace post
} // namespace base

//------------------------------------------------------------------------------
/** Computation of the norm of the error on an element.
 *  Given the approximate (i.e. FE) solution and a reference solution, the
 *  norm of the error is defined as
 *  \f[
 *      | u - u^h |_{\Omega,s}
 *               = [ \int_\Omega (D^{(s)} u - D^{(s)} u^h)^2 dx ]^{1/2}
 *  \f]
 *  where \f$ s \f$ is the order of the derivative to control and \f$ D^{(s)}\f$
 *  the corresponding differential operator. In order to sum up the individual
 *  element contributions, we make use of
 *  \f[
 *      | u - u^h |^2_{\Omega,s} = \sum_{\tau \in \Omega} | u - u^h |^2_{\tau,s}
 *  \f]
 *  Therefore, this object returns the square of the norm of the error for a
 *  specific element. Note that in the case of \f$ s = 0 \f$ we can speak of 
 *  the \f$ L_2 \f$ norm and in case of \f$ s = 1 \f$ it will be the
 *  \f$ H^1 \f$ semi-norm.
 *
 *  \tparam GEOMELEMENT  Mesh element
 *  \tparam FIELDELEMENT Field element
 *  \tparam ORDER   Sobolev-Norm order \f$ s \f$
 */
template<typename GEOMELEMENT, typename FIELDELEMENT, unsigned ORDER>
class base::post::ErrorNorm
    : public boost::function<void( const
                                   base::asmb::FieldElementPointerTuple<GEOMELEMENT*,
                                   FIELDELEMENT*>&,
                                   const typename FIELDELEMENT::ShapeFun::VecDim&,
                                   const double,
                                   double&)>

{
 public:
    //! @name Template parameter
    //@{
    typedef GEOMELEMENT GeomElement;
    typedef FIELDELEMENT FieldElement;
    static const unsigned order   = ORDER;
    //@}
        
    //! @name Deduced attributes
    //@{    
    static const unsigned dim     = GeomElement::Node::dim;
    static const unsigned doFSize = FieldElement::DegreeOfFreedom::size;
    //@}

    //! Result type of the approximation and reference solution
    typedef typename detail_::ValueBinder<dim,doFSize,order>::Type ValueType;

    //! Arguemnt type of the reference solution
    typedef typename base::Vector<dim,double>::Type VecDim;

    //! Type of reference function
    typedef boost::function< ValueType( const VecDim& ) >     Reference;

    //! Type of approximation call
    typedef detail_::FieldEvaluationBinder<dim,doFSize,order> Approximation;
    
    //! Construction with approximate and reference solutions and quadrature
    ErrorNorm( const Reference& reference )
    : reference_( reference )
    { }

    //! Tuple of element pointers
    typedef
    base::asmb::FieldElementPointerTuple<GeomElement*, FieldElement*> FEPT;

    //! Conversion from surface to domain policy
    typedef detail_::SurfacePolicy<GeomElement,FieldElement::dim,GeomElement::dim>
    SurfacePolicy;

    //--------------------------------------------------------------------------
    /** Main function, evaluates solutions and sums upt the differences for each
     *  quadrature point. 
     *  \param[in] fept             Pair of mesh and field element pointers
     *  \param[in] xi               Local evaluation coordinate
     *  \param[in] weight           Quadrature weight
     *  \param[in,out] errorSquared Sum of squares of errors
     */
    void operator()( const FEPT& fept, 
                     const typename GeomElement::GeomFun::VecDim& xi,
                     const double weight,
                     double& errorSquared ) const
    {
        // Access to the mesh and field element pointer tuples
        const GeomElement*  geomEPtr  = fept.geomElementPtr();
        const FieldElement* fieldEPtr = fept.testElementPtr();
        
        // get global coordinate
        const typename GEOMELEMENT::Node::VecDim x =
            base::Geometry<GEOMELEMENT>()( geomEPtr, xi );

        // get jacobian
        const double detJ = base::Jacobian<GEOMELEMENT>()( geomEPtr, xi );

        // get reference
        const ValueType ref = reference_( x );
                    
        // get approximation
        const ValueType approx =
            Approximation::apply( SurfacePolicy::domainElementPtr( geomEPtr ),
                                  fieldEPtr,
                                  SurfacePolicy::makeParametric( geomEPtr, xi ) );
        // point-wise error
        const ValueType error = ref - approx;

        // sum-up
        errorSquared += base::dotProduct( error, error ) * detJ * weight;

        return;
    }
    
private:
    const Reference& reference_;     //!< Reference solution
};


#endif
