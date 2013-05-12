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
// base/post includes
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace post{

        template<typename QUAD, unsigned DIM, unsigned DOFSIZE, unsigned ORDER>
        class ErrorNorm;

        //----------------------------------------------------------------------
        /**
         */
        template<unsigned ORDER, typename QUADRATURE, typename MESH,
                 typename FIELD>
        double errorComputation( const QUADRATURE& quadrature,
                                 const MESH&       mesh,
                                 const FIELD&      field,
                                 const typename
                                 base::post::ErrorNorm<QUADRATURE,
                                                       MESH::Node::dim,
                                                       FIELD::DegreeOfFreedom::size,
                                                       ORDER>::Reference& refSol )
        {
            typedef ErrorNorm<QUADRATURE,
                              MESH::Node::dim,
                              FIELD::DegreeOfFreedom::size,ORDER> Error;
            Error error( refSol, quadrature );

            // compute for every element the error
            typename MESH::ElementPtrConstIter  elemIter  = mesh.elementsBegin();
            typename MESH::ElementPtrConstIter  elemLast  = mesh.elementsEnd();
            typename FIELD::ElementPtrConstIter fieldElem = field.elementsBegin();

            double errorSquared = 0.;
    
            for ( ; elemIter != elemLast; ++elemIter, ++fieldElem ) {
                errorSquared += error.compute( *elemIter, *fieldElem );
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
 *  \tparam REF     Functor providing the reference solution
 *  \tparam QUAD    Quadrature
 *  \tparam ORDER   Sobolev-Norm order \f$ s \f$
 */
template<typename QUAD, unsigned DIM, unsigned DOFSIZE, unsigned ORDER>
class base::post::ErrorNorm
{
 public:
    //! @name Template parameter
    //@{
    typedef QUAD    Quadrature;
    static const unsigned dim     = DIM;
    static const unsigned doFSize = DOFSIZE;
    static const unsigned order   = ORDER;
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
    ErrorNorm( const Reference&  reference,
               const Quadrature& quadrature )
        : reference_( reference ),
          quadrature_( quadrature )
    { }
    
    //--------------------------------------------------------------------------
    /** Main function, evaluates solutions and sums upt the differences for each
     *  quadrature point. 
     *  \param[in] geomElemPtr  Pointer to geometry element
     *  \param[in] fieldElemPtr Pointer to field element
     *  \return          Value of the error squared on this element
     */
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    double compute( const GEOMELEMENT*  geomElemPtr,
                    const FIELDELEMENT* fieldElemPtr ) const
    {
        // sanity checks
        {
            STATIC_ASSERT_MSG( (dim==GEOMELEMENT::Node::dim),
                               "Dimensions do not fit!");
            STATIC_ASSERT_MSG( (doFSize==FIELDELEMENT::DegreeOfFreedom::size),
                               "DoF sizes do not fit!");
        }

        // value of the error in this element
        double error2 = 0.;

        // do the quadrature loop
        typename Quadrature::Iter qIter = quadrature_.begin();
        typename Quadrature::Iter qEnd  = quadrature_.end();
        for ( ; qIter != qEnd; ++qIter ) {

            // local evaluation point
            const typename Quadrature::VecDim xi = qIter -> second;
            
            // get global coordinate
            const typename GEOMELEMENT::Node::VecDim x =
                base::Geometry<GEOMELEMENT>()( geomElemPtr, xi );

            // get jacobian
            const double detJ = base::Jacobian<GEOMELEMENT>()( geomElemPtr, xi );

            // get reference
            const ValueType ref = reference_( x );

            // get approximation
            const ValueType approx = Approximation::apply( geomElemPtr,
                                                           fieldElemPtr,
                                                           xi );
            // point-wise error
            const ValueType error = ref - approx;

            // sum-up
            error2 += base::dotProduct( error, error ) * detJ * (qIter -> first);
            
        }
            
        return error2;
    }
    
    
private:
    const Reference&     reference_;     //!< Reference   solution
    const Quadrature&    quadrature_;    //!< Quadrature object
};


#endif
