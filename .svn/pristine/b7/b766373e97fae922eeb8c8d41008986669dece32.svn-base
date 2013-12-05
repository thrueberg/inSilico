//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   sfun/TensorProduct.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_sfun_tensorproduct_hpp
#define base_sfun_tensorproduct_hpp

//------------------------------------------------------------------------------
// boost  includes
#include <boost/array.hpp>
#include <boost/type_traits.hpp>
// base   includes
#include <base/verify.hpp>
#include <base/meta.hpp>
// base/sfun includes
#include <base/sfun/ShapeFunTraits.hpp>
// base/mesh includes
#include <base/mesh/HierarchicOrder.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace sfun{

        template<typename SFUN, unsigned DIM, base::sfun::Ordering ORDERING>
        class TensorProduct;

        namespace detail_{

            //! Helper for the recursive evaluation
            template<typename SFUN, unsigned DIM>
            struct TensorProductEvaluation;

            //! Helper for lexicographic to hiearchic reordering
            template<typename SFUN, unsigned DIM>
            struct MakeHierarchic
            {
                static const base::Shape shape =
                    base::HyperCubeShape<DIM>::value;
                
                static const unsigned numFun =
                    base::MToTheN<SFUN::numFun,DIM>::value;
                
                typedef base::mesh::HierarchicOrder<shape,SFUN::degree> Hierarchic;

                template<typename ARRAY>
                static void apply( const ARRAY& in, ARRAY& out)
                {
                    for ( unsigned n = 0; n < numFun; n++ ) {
                        const unsigned hier = Hierarchic::apply( n );
                        out[ hier ] = in[ n];
                    }
                }
            };

            //! Just copy in- to output
            struct PlainCopy
            {
                template<typename ARRAY>
                static void apply( const ARRAY& in, ARRAY& out)
                {
                    out = in;
                }
            };
                
        }
        
    }
}

//------------------------------------------------------------------------------
/** Representation of tensor-product shape functions
 *
 *  \tparam SFUN Underlying one-dimensional function type
 *  \tparam DIM  Spatial dimension of the tensor-product
 */
template<typename SFUN, unsigned DIM, base::sfun::Ordering ORDERING>
class base::sfun::TensorProduct
{
    STATIC_ASSERT_MSG( (DIM>0), "Nonsense" );
    
public:
    //! @name Template parameter
    //@{
    typedef SFUN ShapeFun1D;
    static const unsigned dim                  = DIM;
    static const base::sfun::Ordering ordering = ORDERING;
    //@}
    
    //! Polynomial degree
    static const unsigned degree = ShapeFun1D::degree;

    //! Number of functions in one dimension
    static const unsigned numFun1D = ShapeFun1D::numFun;

    //! Number of shape functions of the resulting tensor product
    static const unsigned numFun   = base::MToTheN<numFun1D,dim>::value;

    //! Underlying shape of reference domain: a hypercube
    static const base::Shape shape = base::HyperCubeShape<dim>::value;

    //! Traits object for type defintions (hard-wired to scalar shape functions)
    typedef base::sfun::ShapeFunResultArrays<dim,numFun,
                                             base::sfun::ScalarShapeFunResult> SFRA;

    //! Hierarchic re-ordering for Lagrange shape functions
    typedef typename base::IfElse< ordering == base::sfun::HIERARCHIC, 
                                   detail_::MakeHierarchic<ShapeFun1D,dim>,
                                   detail_::PlainCopy >::Type   Reordering;

    //--------------------------------------------------------------------------
    //! @name Use types from traits object
    //!{
    typedef typename SFRA::VecDim                VecDim;
    typedef typename SFRA::FunArray              FunArray;    
    typedef typename SFRA::GradArray             GradArray;   
    typedef typename SFRA::HessianArray          HessianArray;
    //@}

    //--------------------------------------------------------------------------
    //! Constructor
    TensorProduct( const ShapeFun1D & shapeFun1D = ShapeFun1D() )
        : shapeFun1D_( shapeFun1D ) { }
    
    //--------------------------------------------------------------------------
    //! Evaluation functions
    //@{
    typedef detail_::TensorProductEvaluation<ShapeFun1D,DIM> TPEvaluator;
    
    //! Plain function evaluation
    void fun( const VecDim& xi, FunArray& values ) const
    {
        FunArray lexi;
        TPEvaluator::fun( shapeFun1D_, xi, lexi );
        Reordering::template apply<FunArray>( lexi, values );
    }
    
    //! Evaluation of the functions' gradients
    void gradient( const VecDim& xi, GradArray& values ) const
    {
        GradArray lexi;
        TPEvaluator::gradient( shapeFun1D_, xi, lexi );
        Reordering::template apply<GradArray>( lexi, values );
    }
    
    //! Evaluation of the functions' Hessians
    void hessian(  const VecDim& xi, HessianArray& values ) const
    {
        HessianArray lexi;
        TPEvaluator::hessian( shapeFun1D_, xi, lexi );
        Reordering::template apply<HessianArray>::template apply( lexi, values );
    }
    //@}

    //--------------------------------------------------------------------------
    //! Provide the support points of the functions
    static void supportPoints( boost::array<VecDim, numFun> & supportPoints )
    {
        boost::array<VecDim,numFun> lexi;

        TPEvaluator::supportPoints( lexi );

        // This line does not compile with g++ 4.7.2-2ubuntu1 ??
        //Reordering::template apply< boost::array<VecDim,numFun> >::apply( lexi, supportPoints );
        //Using ADL instead:
        Reordering::apply( lexi, supportPoints );
    }

private:
    //! Local object for one-dimensional shape functions
    const ShapeFun1D & shapeFun1D_; 
};

//------------------------------------------------------------------------------
namespace base{
    namespace sfun{
        namespace detail_{

            //------------------------------------------------------------------
            //! Evaluation of 1D-Tensor-Product (delegates calls to SFUN)
            template<typename SFUN>
            struct TensorProductEvaluation<SFUN,1>
            {
                static const unsigned numFun = SFUN::numFun;

                typedef base::sfun::ShapeFunResultArrays<1,numFun,
                                                         base::sfun::ScalarShapeFunResult> SFRA;
                
                
                static void fun( const SFUN& shapeFun1D,
                                 const typename SFRA::VecDim& xi,
                                 typename SFRA::FunArray& values )
                {
                    shapeFun1D.fun( xi, values );
                }

                static void gradient( const SFUN& shapeFun1D,
                                      const typename SFRA::VecDim& xi,
                                      typename SFRA::GradArray& values )
                {
                    shapeFun1D.gradient( xi, values );
                }

                static void hessian( const SFUN& shapeFun1D,
                                     const typename SFRA::VecDim& xi,
                                     typename SFRA::HessianArray& values )
                {
                    shapeFun1D.hessian( xi, values );
                }

                static void supportPoints( boost::array<typename SFRA::VecDim,
                                                        numFun> & supportPoints )
                {
                    SFUN::supportPoints( supportPoints );
                } 

            };

            //------------------------------------------------------------------
            //! Evaluate of D-Tensor-Product (recursive calls)
            template<typename SFUN, unsigned DIM>
            struct TensorProductEvaluation
            {
                //! Number of total shape functions
                static const unsigned numFun = base::MToTheN<SFUN::numFun,DIM>::value;

                //! Traits object for type definitions
                typedef base::sfun::ShapeFunResultArrays<DIM,numFun,
                                                         base::sfun::ScalarShapeFunResult> SFRA;
  
                //! Object for the evaluation of DIM-1 Tensor product
                typedef TensorProductEvaluation<SFUN,DIM-1> LowerDimEvaluation;

                //--------------------------------------------------------------
                //! Evaluate shape functions as tensor-products
                static void fun( const SFUN& shapeFun1D,
                                      const typename SFRA::VecDim& xi,
                                      typename SFRA::FunArray& values )
                {
                    // get values from lower dimension
                    typename LowerDimEvaluation::SFRA::FunArray lowerDimValues;
                    const typename LowerDimEvaluation::SFRA::VecDim lowerXi =
                        base::head<DIM-1>( xi );
                    
                    LowerDimEvaluation::fun( shapeFun1D, lowerXi, lowerDimValues );

                    // get values from one-dimensional quadrature
                    typename SFUN::FunArray oneDimValues;
                    const typename SFUN::VecDim oneDXi = base::tail<1>( xi );

                    shapeFun1D.fun( oneDXi, oneDimValues );

                    // construct results
                    unsigned ctr = 0;
                    for ( unsigned nOuter = 0;
                          nOuter < SFUN::numFun; nOuter++ ) {

                        for ( unsigned nInner = 0;
                              nInner < LowerDimEvaluation::numFun; nInner ++ ) {

                            // Evaluation of the tensor-product
                            values[ ctr++ ] = lowerDimValues[nInner] * oneDimValues[nOuter];
                        }
                    }

                    return;
                }

                //--------------------------------------------------------------
                //! Evaluate shape function gradients as tensor-products
                static void gradient( const SFUN& shapeFun1D,
                                      const typename SFRA::VecDim& xi,
                                      typename SFRA::GradArray& values )
                {
                    // get values from lower dimension
                    typename LowerDimEvaluation::SFRA::FunArray     lowerDimValues;
                    typename LowerDimEvaluation::SFRA::GradArray    lowerDimGradients;
                    const typename LowerDimEvaluation::SFRA::VecDim lowerXi
                        = base::head<DIM-1>( xi );
                    
                    LowerDimEvaluation::fun( shapeFun1D, lowerXi, lowerDimValues );
                    LowerDimEvaluation::gradient( shapeFun1D, lowerXi, lowerDimGradients );

                    // get values from one-dimensional quadrature
                    typename SFUN::FunArray     oneDimValues;
                    typename SFUN::GradArray    oneDimGradients;
                    const typename SFUN::VecDim oneDXi = base::tail<1>( xi );

                    shapeFun1D.fun( oneDXi, oneDimValues );
                    shapeFun1D.gradient( oneDXi, oneDimGradients );

                    // construct results
                    unsigned ctr = 0;
                    for ( unsigned nOuter = 0;
                          nOuter < SFUN::numFun; nOuter++ ) {

                        for ( unsigned nInner = 0;
                              nInner < LowerDimEvaluation::numFun; nInner ++ ) {

                            // Evaluation of the tensor-product
                            for ( unsigned d = 0; d < DIM-1; d ++ ) {
                                values[ctr][d] =
                                    oneDimValues[ nOuter ] *
                                    lowerDimGradients[nInner][d];
                            }

                            values[ctr][DIM-1] =
                                oneDimGradients[nOuter][0] *
                                lowerDimValues[ nInner ];

                            ctr++;
                        }
                    }
                    
                    return;
                }

                //--------------------------------------------------------------
                //! Evaluate shape function Hessians as tensor-products
                static void hessian( const SFUN & shapeFun1D,
                                     const typename SFRA::VecDim & xi,
                                     typename SFRA::HessianArray & values ) 
                {
                    // get values from lower dimension
                    typename LowerDimEvaluation::SFRA::FunArray     lowerDimValues;
                    typename LowerDimEvaluation::SFRA::GradArray    lowerDimGradients;
                    typename LowerDimEvaluation::SFRA::HessianArray lowerDimHessians;
                    const typename LowerDimEvaluation::SFRA::VecDim lowerXi
                        = base::head<DIM-1>( xi );
                    
                    LowerDimEvaluation::fun( shapeFun1D, lowerXi, lowerDimValues );
                    LowerDimEvaluation::gradient( shapeFun1D, lowerXi, lowerDimGradients );
                    LowerDimEvaluation::hessian(  shapeFun1D, lowerXi, lowerDimHessians );

                    // get values from one-dimensional quadrature
                    typename SFUN::FunArray     oneDimValues;
                    typename SFUN::GradArray    oneDimGradients;
                    typename SFUN::HessianArray oneDimHessians;
                    const typename SFUN::VecDim oneDXi = base::tail<1>( xi );

                    shapeFun1D.fun( oneDXi, oneDimValues );
                    shapeFun1D.gradient( oneDXi, oneDimGradients );
                    shapeFun1D.hessian(  oneDXi, oneDimHessians );

                    // construct results
                    unsigned ctr = 0;
                    for ( unsigned nOuter = 0;
                          nOuter < SFUN::numFun; nOuter++ ) {

                        for ( unsigned nInner = 0;
                              nInner < LowerDimEvaluation::numFun; nInner ++ ) {

                            // Evaluation of the tensor-product
                            for ( unsigned d1 = 0; d1 < DIM-1; d1 ++ ) {
                                
                                for ( unsigned d2 = 0; d2 < DIM-1; d2 ++ ) {

                                    values[ctr]( d1, d2 ) =
                                        oneDimValues[ nOuter ] *
                                        lowerDimHessians[ nInner ]( d1, d2 );
                                }

                                values[ctr]( d1, DIM-1 ) =
                                    oneDimGradients( 0, nOuter ) *
                                    lowerDimGradients( d1, nInner );

                                values[ctr]( DIM-1, d1 ) =
                                    oneDimGradients( 0, nOuter ) *
                                    lowerDimGradients( d1, nInner );

                            }

                            values[ctr]( DIM-1, DIM-1 ) =
                                oneDimHessians[nOuter](0,0) *
                                lowerDimValues( nInner );
                            
                            ctr++;
                        }
                    }

                    return;
                }

                //--------------------------------------------------------------
                //! Construct shape functions' support points via tensor-product
                static void supportPoints( boost::array<typename SFRA::VecDim,
                                                        numFun> & supportPoints )
                {
                    // get values from lower dimension
                    boost::array< typename LowerDimEvaluation::SFRA::VecDim,
                                  LowerDimEvaluation::numFun> lowerDimSupportPoints;
                    LowerDimEvaluation::supportPoints( lowerDimSupportPoints );

                    // get values from one-dimensional quadrature
                    boost::array< typename SFUN::VecDim,
                                  SFUN::numFun> oneDimSupportPoints;
                    SFUN::supportPoints( oneDimSupportPoints );

                    // construct results
                    unsigned ctr = 0;
                    for ( unsigned nOuter = 0;
                          nOuter < SFUN::numFun; nOuter++ ) {

                        for ( unsigned nInner = 0;
                              nInner < LowerDimEvaluation::numFun; nInner ++ ) {

                            for ( unsigned d = 0; d < DIM-1; d++ )
                                supportPoints[ctr][d] =
                                    lowerDimSupportPoints[nInner][d];

                            supportPoints[ctr][DIM-1] =
                                oneDimSupportPoints[nOuter][0];

                            ctr++;
                        }
                    }

                    return;
                }
                

            };
            
        } // end namespace detail_
        
    }
}

#endif
