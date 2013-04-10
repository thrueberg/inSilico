//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   quad/TensorProduct.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_quad_tensorproduct_hpp
#define base_quad_tensorproduct_hpp

//------------------------------------------------------------------------------
// system includes
#include <utility>
// boost  includes
#include <boost/array.hpp>
// base   includes
#include <base/linearAlgebra.hpp>
#include <base/meta.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace quad{

        template<typename QUAD, unsigned DIM>
        class TensorProduct;

        namespace detail_{
            
            template<typename QUAD, unsigned DIM>
            struct TensorProductConstructor;
        }
    }
}

//------------------------------------------------------------------------------
/** \brief   %Quadrature for tensor-product (0,1)^d elements (hypercubes)
 *  \details  Wrapper for any one-dimensional (0,1)-quadrature rule which 
 *   generates the points and weights for the corresponding d-fold tensor-product
 *  \param QUAD Underlying one-dimensional (0,1)-quadrature rule
 *  \param DIM  Dimension d of the d-fold tensor-product
 */
template<typename QUAD, unsigned DIM>
class base::quad::TensorProduct
{
    // Sanity check
    STATIC_ASSERT_MSG( (QUAD::dim == 1),
                       "Only one-dimensional rules can be used here" );
    
public:
    //! @name Template parameter
    //@{
    typedef QUAD  Quadrature1D;      //!< Underlying 1D quadrature
    static const unsigned dim = DIM; //!< dimension of the tensor-product
    //@}

    //! Number of points of 1D quadrature
    static const unsigned numPoints1D = QUAD::numPoints;
    //! Number of points of tensor-product quadrature
    static const unsigned numPoints   = base::MToTheN<numPoints1D,dim>::value;

    //! Local coordinate vector
    typedef typename base::VectorType<dim>::Type     VecDim;

    //! Iterator for external access
    typedef typename boost::array< std::pair<double, VecDim>,
                                   numPoints >::const_iterator  Iter;

    //! Constructor which generates the array of (weight,point)-pairs
    TensorProduct()
    {
        // Use constructor helper object defined below
        base::quad::detail_::TensorProductConstructor<Quadrature1D,
                                                      dim> tpc;
        tpc( weightsAndPoints_ );
    }

    //--------------------------------------------------------------------------
    //! @name Accessor
    //@{
    //! Begin of array iterator
    Iter begin() const  { return weightsAndPoints_.begin(); }
    //! End of array iterator
    Iter end()   const  { return weightsAndPoints_.end();   }
    //@}

private:
    //! Array of weight and point pairs for the tensor product rule
    boost::array<std::pair<double,VecDim>, numPoints> weightsAndPoints_;
};

//------------------------------------------------------------------------------
namespace base{
    namespace quad{
        namespace detail_{

            //------------------------------------------------------------------
            //! Constructor for one-dimensional 'tensor'-product rule 
            template<typename QUAD>
            struct TensorProductConstructor<QUAD,1>
            {
                typedef typename base::VectorType<1>::Type Vec1;

                //! Construct 1D-Tensor product rule (dummy)
                void operator()( boost::array<std::pair<double, Vec1>, 
                                              QUAD::numPoints> & weightsAndPoints  )
                {
                    //! construct underlying one-dimensional quadrature rule
                    QUAD quad1D;
                    
                    //! do a simple copy of the 1D rule
                    std::copy( quad1D.begin(), quad1D.end(),
                               weightsAndPoints.begin() );
                    
                    return;
                }
            };
            
            //------------------------------------------------------------------
            //! Constructor for general dimensions
            template<typename QUAD, unsigned DIM>
            struct TensorProductConstructor
            {
                typedef typename base::VectorType<DIM>::Type VecDim;

                //! Overloaded function call operator generates the quadrature rule
                void operator()( boost::array<std::pair<double, VecDim>, 
                                              base::MToTheN<QUAD::numPoints,
                                                              DIM>::value> &
                                 weightsAndPoints )
                {
                    // lower-dimensional quadrature
                    TensorProduct<QUAD,DIM-1> lowerQuad;
                    const unsigned lowerNumPoints = TensorProduct<QUAD,DIM-1>::numPoints;
                    
                    // copy lower-dimensional quadrature into arrays 
                    for ( unsigned outer = 0; outer < QUAD::numPoints; outer ++ ) {
                        //! Counter
                        unsigned index = outer * lowerNumPoints;
                    
                        for ( typename TensorProduct<QUAD,DIM-1>::Iter lower = lowerQuad.begin();
                              lower != lowerQuad.end(); ++lower ) {
                            
                            // copy weights
                            weightsAndPoints[ index ].first = lower -> first;
                        
                            // copy points into [0,DIM-1) subrange of the DIM-points
                            weightsAndPoints[index].second.template head<DIM-1>()
                                = lower -> second;
                        
                            index ++;
                        }
                    }

                    // underlying one-dimensional quadrature rule
                    QUAD quad1D;

                    // make outer product for the new rule
                    unsigned outer = 0;
                    for ( typename QUAD::Iter iter = quad1D.begin();
                          iter !=  quad1D.end();
                          ++iter, outer++ ) {
                    
                        for ( unsigned inner = 0; inner < lowerNumPoints; inner ++ ) {
                        
                            const unsigned index = outer * lowerNumPoints + inner;
                            weightsAndPoints[index].first *= iter -> first;
                            weightsAndPoints[index].second.template tail<1>()
                                = iter -> second;
                        
                        }
                    }
                }

            };
            // End of struct
            //------------------------------------------------------------------
        }
    }
}

#endif
