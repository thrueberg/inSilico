//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   evaluateOnCutCells.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_cut_evaluateoncutcells_hpp
#define base_cut_evaluateoncutcells_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// boost includes
#include <boost/bind.hpp>
#include <boost/function.hpp>
// base  includes
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>


//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        //-----------------------------------------------------------------------
        template<typename VALUETYPE, typename MESH, typename FIELD, typename OP,
                 typename CELL>
        void evaluateAtCutCellNodes( const MESH& mesh, const FIELD& field,
                                     const OP& fieldEvaluator,
                                     const std::vector<CELL>& cells,
                                     std::vector<VALUETYPE>& nodalValues );

        template<typename FIELDTUPLEBINDER, 
                 typename VALUETYPE, typename FIELDBINDER, typename OP,
                 typename CELL>
        void evaluateInCutCells( const FIELDBINDER& fieldBinder,
                                 const OP& fieldEvaluator,
                                 const std::vector<CELL>& cells,
                                 std::vector<VALUETYPE>& cellValues,
                                 const bool inOut = true );



        //-----------------------------------------------------------------------
        //! Short cut for pure field evaluation at the cut cell nodes
        template<typename MESH, typename FIELD, typename CELL>
        void evaluateFieldAtCutCellNodes(
            const MESH& mesh,
            const FIELD& field,
            const std::vector<CELL>& cells,
            std::vector<typename base::Vector<FIELD::DegreeOfFreedom::size,number>::Type>&
            nodalValues )
        {
            // Functor for the evaluation of the field
            boost::function<
                typename base::Vector<FIELD::DegreeOfFreedom::size,number>::Type
                ( const typename MESH::Element*,
                  const typename FIELD::Element*,
                  const typename FIELD::Element::FEFun::VecDim& )>
                evaluateField =
                boost::bind( &base::post::evaluateField<typename MESH::Element,
                             typename FIELD::Element>, _1, _2, _3 );

            // call evaluation function
            typedef typename
                base::Vector<FIELD::DegreeOfFreedom::size,number>::Type ValueType;
            evaluateAtCutCellNodes<ValueType>( mesh, field, evaluateField,
                                               cells, nodalValues );

            return;
        }

        // helper for signed distance function evaluation
        namespace detail_{
            template<typename MESH>
            struct EvaluateDistanceAtCutCellNodes
            {
                static double apply(
                    const typename MESH::Element* geomEp,
                    const typename MESH::Element* dummy,
                    const typename MESH::Element::GeomFun::VecDim& xi,
                    const std::vector<typename base::cut::LevelSet<MESH::Node::dim> >&
                    levelSets )
                {
                    return base::cut::signedDistance( geomEp, xi, levelSets );
                }

            };
        }

        // evaluate the signed distance at cut cell nodes
        template<typename MESH, typename CELL>
        void evaluateCutCellNodeDistances( const MESH& mesh,
                                           const std::vector<
                                           typename base::cut::LevelSet<MESH::Node::dim> >&
                                           levelSets,
                                           const std::vector<CELL>& cells,
                                           std::vector<double>& result )
        {
            boost::function<double( const typename MESH::Element*,
                                    const typename MESH::Element*,
                                    const typename MESH::Element::GeomFun::VecDim )>
                evaluateDistance
                = boost::bind(
                    &detail_::EvaluateDistanceAtCutCellNodes<MESH>::apply, _1, _2, _3,
                    boost::ref( levelSets ) );
            
            base::cut::evaluateAtCutCellNodes( mesh, mesh,
                                               evaluateDistance, cells,
                                               result );
        }
    
    }
}

//------------------------------------------------------------------------------
/** Evaluate a functor at the nodes of the cut cells.
 *  Provided is a functor of type
 *  \code{.cpp}
 *  VALUETYPE( const MESH::Element*, const FIELD::Element*,
 *             const FIELD::Element::FEFun::VecDim& )
 *  \endcode
 *  and it will be applied to all nodes of all cells that are cut. This
 *  functionality allows to evaluate the field or some other related entity
 *  at those nodes and visualise the datum in cut cells.
 *
 *  \tparam VALUETYPE  Result type of the functor
 *  \tparam MESH       Type of bulk mesh
 *  \tparam FIELD      Type of bulk field
 *  \tparam OP         Type of evaluation functor
 *  \tparam CELL       Type of cell
 *  \param[in]  mesh           Bulk mesh
 *  \param[in]  field          Bulk field
 *  \param[in]  fieldEvaluator Evaluation functor
 *  \param[in]  cells          Element cell structures
 *  \param[out] nodalValues    Result of the evaluation
 */
template<typename VALUETYPE, typename MESH, typename FIELD, typename OP,
         typename CELL>
void base::cut::evaluateAtCutCellNodes( const MESH& mesh, const FIELD& field,
                                        const OP& fieldEvaluator,
                                        const std::vector<CELL>& cells,
                                        std::vector<VALUETYPE>& nodalValues )
{

    // go through all elements and cells
    const std::size_t numElements = cells.size();
    for ( std::size_t e = 0; e < numElements; e++ ) {

        // only act on cut cells
        if ( cells[e].isCut() ) {

            // geometry and field elements
            const typename MESH::Element*  geomEp  = mesh.elementPtr(  e );
            const typename FIELD::Element* fieldEp = field.elementPtr( e );

            // extract nodes from cut cell
            std::vector<typename CELL::VecDim> cutCellNodes;
            cells[e].getNodes( cutCellNodes );

            // go through all nodes of the cut cell
            for ( std::size_t n = 0; n < cutCellNodes.size(); n++ ) {

                nodalValues.push_back(
                    fieldEvaluator( geomEp, fieldEp, 
                                    cutCellNodes[n] ) );

            }

        } // if cell is cut
    } // end loop over elements

    return;
}

//------------------------------------------------------------------------------
/** Evaluate a functor in cut cells.
 *  Provided is a functor of type
 *  \code{.cpp}
 *  VALUETYPE( const MESH::Element*, const FIELD::Element*,
 *             const FIELD::Element::FEFun::VecDim& )
 *  \endcode
 *  and it will be applied to all cells that are cut. This
 *  functionality allows to evaluate the field or some other related entity
 *  at the cells and write the result to a file. For every cut cell, all
 *  sub-elements on the in- or outside part will be queried and the functor
 *  applied at the geometric centre of that sub-element.
 *
 *  \tparam VALUETYPE  Result type of the functor
 *  \tparam OP         Type of evaluation functor
 *  \tparam CELL       Type of cell
 *  \param[in]  mesh           Bulk mesh
 *  \param[in]  field          Bulk field
 *  \param[in]  fieldEvaluator Evaluation functor
 *  \param[in]  cells          Element cell structures
 *  \param[out] cellValues     Result of the evaluation
 *  \param[in]  inOut          True for inside evaluation, otherwise outside
 */
template<typename FIELDTUPLEBINDER, 
         typename VALUETYPE, typename FIELDBINDER, typename OP,
         typename CELL>
void base::cut::evaluateInCutCells( const FIELDBINDER& fieldBinder,
                                    const OP& fieldEvaluator,
                                    const std::vector<CELL>& cells,
                                    std::vector<VALUETYPE>& cellValues,
                                    const bool inOut )
{
    // Simplex center
    const typename FIELDBINDER::Mesh::Element::GeomFun::VecDim centre =
        base::ShapeCentroid<base::SimplexShape<FIELDBINDER::Mesh::Node::dim>::value>::apply();
    
    // go through all elements and cells
    typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
    typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
    for ( std::size_t e = 0; iter != end; ++iter, e++ ) {

        // only act on cut cells
        if ( cells[e].isCut() ) {

            // number of sub-elements per cell
            const std::size_t numSubCells =
                (inOut ?
                 cells[e].numVolumeInElements() :
                 cells[e].numVolumeOutElements() );

            for ( unsigned c = 0; c < numSubCells; c++ ) {
                
                // create local coordinate
                const typename CELL::VecDim eta =
                    cells[e].mapVolumeCoordinate( centre, c, inOut );

                // evaluate
                cellValues.push_back( fieldEvaluator( FIELDTUPLEBINDER::makeTuple(*iter),
                                                      eta ) );
                
            } // end loop over sub-elements
            
        } // if cell is cut
    } // end loop over elements

    return;
}

#endif
    
