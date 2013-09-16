//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   evaluateField.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_post_evaluatefield_hpp
#define base_post_evaluatefield_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
// base includes
#include <base/linearAlgebra.hpp>
// base/mesh includes
#include <base/mesh/sampleStructured.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace post{

        //----------------------------------------------------------------------
        // Evaluate the field at a given point
        template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
        typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                  base::number>::Type
        evaluateFieldHistory( const GEOMELEMENT*  geomElemPtr,
                              const FIELDELEMENT* fieldElemPtr,
                              const typename FIELDELEMENT::FEFun::VecDim& xi );

        // Evaluate the field at a given point
        template<typename GEOMELEMENT, typename FIELDELEMENT>
        typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,
                                  base::number>::Type
        evaluateField( const GEOMELEMENT*  geomElemPtr,
                       const FIELDELEMENT* fieldElemPtr,
                       const typename FIELDELEMENT::FEFun::VecDim& xi )
        {
            return evaluateFieldHistory<0>( geomElemPtr, fieldElemPtr, xi );
        }

        //----------------------------------------------------------------------
        // Evaluate the gradient of the field at a given point
        template<unsigned HIST, typename GEOMELEMENT,typename FIELDELEMENT>
        static typename base::Matrix<GEOMELEMENT::Node::dim,
                                         FIELDELEMENT::DegreeOfFreedom::size,
                                         base::number>::Type
        evaluateFieldGradientHistory( const GEOMELEMENT*  geomElemPtr,
                                      const FIELDELEMENT* fieldElemPtr,
                                      const typename FIELDELEMENT::FEFun::VecDim& xi );

        
        // Evaluate the gradient of the field at a given point
        template<typename GEOMELEMENT,typename FIELDELEMENT>
        static typename base::Matrix<GEOMELEMENT::Node::dim,
                                         FIELDELEMENT::DegreeOfFreedom::size,
                                         base::number>::Type
        evaluateFieldGradient( const GEOMELEMENT*  geomElemPtr,
                               const FIELDELEMENT* fieldElemPtr,
                               const typename FIELDELEMENT::FEFun::VecDim& xi )
        {
            return evaluateFieldGradientHistory<0>( geomElemPtr, fieldElemPtr, xi );
        }

        //----------------------------------------------------------------------
        namespace detail_{

            template<typename GRID, typename FIELD, typename OUTITER, unsigned ORDER>
            struct FieldSampling;

            template<typename GRID, typename FIELD, typename OUTITER>
            struct FieldSampling<GRID,FIELD,OUTITER,0>
            {
                typedef typename GRID::Element Element;
                
                static void apply( const std::size_t eIndex,
                                   const typename base::GeomTraits<Element>::LocalVecDim& xi,
                                   const GRID& grid, const FIELD& field,
                                   OUTITER outIter )
                {
                    const Element*                 const gep = grid.elementPtr( eIndex );
                    const typename FIELD::Element* const fep = field.elementPtr( eIndex );
                    *outIter++ = evaluateField( gep, fep, xi );
                }
            };

            template<typename GRID, typename FIELD, typename OUTITER>
            struct FieldSampling<GRID,FIELD,OUTITER,1>
            {
                typedef typename GRID::Element Element;

                static void apply( const std::size_t eIndex,
                                   const typename base::GeomTraits<Element>::LocalVecDim& xi,
                                   const GRID& grid, const FIELD& field,
                                   OUTITER outIter )
                {
                    const Element*                 const gep = grid.elementPtr( eIndex );
                    const typename FIELD::Element* const fep = field.elementPtr( eIndex );
                    *outIter++ = evaluateFieldGradient( gep, fep, xi );
                }
            };

        } // namespace detail_

        //----------------------------------------------------------------------
        //! Special case: evaluate a field on the grid
        template<typename GRID, typename FIELD, typename OUTITER>
        void sampleField( const GRID&    grid, const FIELD&   field, 
                          OUTITER        outputIter,
                          const unsigned resolution = 1,
                          const bool discont = false )
        {
            typename base::mesh::detail_::GridSampler<typename GRID::Element>::Type
                sampler =
                boost::bind( &detail_::FieldSampling<GRID,FIELD,OUTITER,0>::apply,
                             _1, _2,
                             boost::ref( grid ), boost::ref( field ), outputIter );

            base::mesh::sampleStructured<GRID::dim>( grid.gridSizes(), sampler,
                                                     resolution, discont );
        }
        
        //! Special case: evaluate a field gradient on the grid
        template<typename GRID, typename FIELD, typename OUTITER>
        void sampleFieldGradient( const GRID&    grid, const FIELD&   field, 
                                  OUTITER        outputIter,
                                  const unsigned resolution = 1,
                                  const bool discont = false )
        {
            typename base::mesh::detail_::GridSampler<typename GRID::Element>::Type
                sampler =
                boost::bind( &detail_::FieldSampling<GRID,FIELD,OUTITER,1>::apply,
                             _1, _2,
                             boost::ref( grid ), boost::ref( field ), outputIter );

            base::mesh::sampleStructured<GRID::dim>( grid.gridSizes(), sampler,
                                                     resolution, discont );
        }

        
    }
}

//------------------------------------------------------------------------------
/** Evaluate the solution field in a given element at a local coordinate.
 *  The FE solution field (or any other field represented by FEBASIS),
 *  has the interpolation representation
 *  \f[
 *       u^h (x) = \sum u^i \phi^i (x(\xi))
 *  \f]
 *  which is evaluated by this object. Other than providing a global coordinate
 *  \f$ x \f$, the pair of a geometry element and its local evaluation point
 *  \f$ \xi \f$ are provided. 
 *  \tparam     GEOMELEMENT  Geometry element
 *  \tparam     FIELDELEMENT Field element
 *  \param[in]  geomElemPtr  Pointer to geometry element
 *  \param[in]  fieldElemPtr Pointer to field element
 *  \param[in]  xi    Local evaluation coordinate
 *  \returns          Value of the approximate field
 */
template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,base::number>::Type
base::post::evaluateFieldHistory( const GEOMELEMENT*  geomElemPtr,
                                  const FIELDELEMENT* fieldElemPtr,
                                  const typename FIELDELEMENT::FEFun::VecDim& xi ) 
{
    // Type of degree of freedom
    typedef typename FIELDELEMENT::DegreeOfFreedom DegreeOfFreedom;
        
    // Return type of this object
    typedef typename base::Vector<DegreeOfFreedom::size,
                                      base::number>::Type ReturnType;

    // get dof-objects for element
    std::vector<DegreeOfFreedom*> elementDoFs;
    std::copy( fieldElemPtr -> doFsBegin(), fieldElemPtr -> doFsEnd(),
               std::back_inserter( elementDoFs ) );

    // get values of the dofs
    std::vector<ReturnType> doFValues;
    for ( unsigned s = 0; s < elementDoFs.size(); s ++ ) {

        ReturnType doFValue;
        for ( unsigned d = 0; d < DegreeOfFreedom::size; d ++ )
            doFValue[d] =
                elementDoFs[s] ->template getHistoryValue<HIST>( d );
        
        doFValues.push_back( doFValue );
    }

    // Evaluate the shape function
    typename FIELDELEMENT::FEFun::FunArray funValues;
    (fieldElemPtr -> fEFun()).evaluate( geomElemPtr, xi, funValues );
        
    // sanity check
    assert( doFValues.size() == funValues.size() );

    // Linear combination yields the result
    ReturnType result = base::constantVector<DegreeOfFreedom::size>( 0. );
    for ( unsigned f = 0; f < funValues.size(); f ++ )
        result += funValues[f] * doFValues[f];

    return result;
}

//------------------------------------------------------------------------------
/** Evaluate the gradient of solution field in a given element at a
 *  local coordinate.
 *  The gradient of the FE solution field at time point \f$ t_{n-s}\f$
 *  has the interpolation representation
 *  \f[
 *       \nabla_x u_{n-1}^h (x) = \sum u_{n-s}^i \nabla_x \phi^i (x(\xi))
 *  \f]
 *  which is evaluated by this object. Other than providing a global coordinate
 *  \f$ x \f$, the pair of a geometry element and its local evaluation point
 *  \f$ \xi \f$ are provided.
 *  \tparam     HIST         Number of history term
 *  \tparam     GEOMELEMENT  Geometry element
 *  \tparam     FIELDELEMENT Field element
 *  \param[in]  geomElemPtr  Pointer to geometry element
 *  \param[in]  fieldElemPtr Pointer to field element
 *  \param[in]  xi           Local evaluation coordinate
 *  \returns                 Value of the approximate field gradient
 */
template<unsigned HIST,typename GEOMELEMENT,typename FIELDELEMENT>
typename base::Matrix<GEOMELEMENT::Node::dim,
                          FIELDELEMENT::DegreeOfFreedom::size,base::number>::Type
base::post::evaluateFieldGradientHistory( const GEOMELEMENT*  geomElemPtr,
                                          const FIELDELEMENT* fieldElemPtr,
                                          const typename FIELDELEMENT::FEFun::VecDim& xi ) 
{
    // Type of degree of freedom
    typedef typename FIELDELEMENT::DegreeOfFreedom DegreeOfFreedom;
        
    //! Value and gradient types
    typedef typename base::Vector<DegreeOfFreedom::size,
                                      base::number>::Type ValueType;

    // get dof-objects for element
    std::vector<DegreeOfFreedom*> elementDoFs;
    std::copy( fieldElemPtr -> doFsBegin(), fieldElemPtr -> doFsEnd(),
               std::back_inserter( elementDoFs ) );

    // get values of the dofs
    std::vector<ValueType> doFValues;
    for ( unsigned s = 0; s < elementDoFs.size(); s ++ ) {

        ValueType doFValue;
        for ( unsigned d = 0; d < DegreeOfFreedom::size; d ++ )
            doFValue[d] = elementDoFs[s] ->template getHistoryValue<HIST>( d );

        doFValues.push_back( doFValue );
    }

    // Evaluate the shape function
    std::vector<typename GEOMELEMENT::Node::VecDim> funGradValues;
    (fieldElemPtr -> fEFun()).evaluateGradient( geomElemPtr, xi, funGradValues );
        
    // sanity check
    assert( doFValues.size() == funGradValues.size() );

    // Linear combination yields the result
    typename base::Matrix<GEOMELEMENT::Node::dim,
                              DegreeOfFreedom::size,base::number>::Type
        result = base::constantMatrix<GEOMELEMENT::Node::dim,
                                      DegreeOfFreedom::size>( 0. );
    for ( unsigned f = 0; f < funGradValues.size(); f ++ )
        result += funGradValues[f] * ( doFValues[f].transpose() );

    return result;
}


#endif
