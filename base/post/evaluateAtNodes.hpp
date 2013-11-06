//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   evaluateAtNodes.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_post_evaluateatnodes_hpp
#define base_post_evaluateatnodes_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <iterator>
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/verify.hpp>
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace post{

        namespace detail_{
            template<typename FIELD>
            struct ValueTypeBinder
            {
                static const unsigned size = FIELD::DegreeOfFreedom::size;
                typedef typename base::Vector<size,number>::Type Type;
            };

        }

        template<unsigned HIST,typename MESH, typename FIELD>
        void evaluateHistoryAtNodes( const MESH& mesh,
                                     const FIELD& field,
                                     std::vector<typename detail_::ValueTypeBinder<FIELD>::Type>&
                                     nodalValues );

        template<typename MESH, typename FIELD>
        void evaluateAtNodes( const MESH& mesh,
                              const FIELD& field,
                              std::vector<typename detail_::ValueTypeBinder<FIELD>::Type>&
                              nodalValues )
        {
            evaluateHistoryAtNodes<0>( mesh, field, nodalValues );
        }

    }
}

//------------------------------------------------------------------------------
/** Evaluate the field at all nodes of the given mesh.
 *  \tparam MESH  Type of mesh
 *  \tparam FIELD Type of field
 */
template<unsigned HIST,typename MESH, typename FIELD>
void base::post::evaluateHistoryAtNodes( const MESH& mesh,
                                         const FIELD& field,
                                         std::vector<typename base::post::detail_::
                                         ValueTypeBinder<FIELD>::Type>& nodalValues )
{
    // DoF size
    static const unsigned doFSize = FIELD::DegreeOfFreedom::size;

    // Value type of the evaluation
    typedef typename base::post::detail_::ValueTypeBinder<FIELD>::Type ValueType;

    // resize the result container
    nodalValues.resize( std::distance( mesh.nodesBegin(), mesh.nodesEnd() ) );

    typename MESH::ElementPtrConstIter  geomElemIter = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter  geomElemLast = mesh.elementsEnd();
    typename FIELD::ElementPtrConstIter fieldElemIter = field.elementsBegin();

    // go through element range
    for ( ; geomElemIter != geomElemLast; ++geomElemIter, ++fieldElemIter ) {

        // storage of nodes IDs
        std::vector<std::size_t> nodeIDs;

        // Collect IDs of element's nodes
        for ( typename MESH::Element::NodePtrConstIter nIter =
                  (*geomElemIter) -> nodesBegin();
              nIter != (*geomElemIter) -> nodesEnd(); ++nIter ) {
            nodeIDs.push_back( (*nIter) -> getID() );
        }

        // Support points of the geometry shape function
        typedef typename MESH::Element::GeomFun::VecDim VecDim;
        boost::array<VecDim,MESH::Element::numNodes> supportPoints;
        MESH::Element::GeomFun::supportPoints( supportPoints );

        // Look up element DoF values
        std::vector<ValueType> doFValues;
        for ( typename FIELD::Element::DoFPtrConstIter dIter =
                  (*fieldElemIter) -> doFsBegin();
              dIter != (*fieldElemIter) -> doFsEnd();
              ++dIter ) {

            ValueType doFValue;
            for ( unsigned d = 0; d < doFSize; d++ )
                doFValue[d] = (*dIter) ->template getHistoryValue<HIST>( d );

            doFValues.push_back( doFValue );
        }

        // Nodal DoFs always by means of field evaluation
        for ( unsigned n = 0; n < MESH::Element::numNodes; n ++ ) {

            // Evaluation point
            const VecDim xi = supportPoints[ n ];

            // Evaluate the shape function
            typename FIELD::Element::FEFun::FunArray funValues;
            ((*fieldElemIter) -> fEFun() ).evaluate( *geomElemIter, xi, funValues );
                                
            // sanity check
            assert( doFValues.size() == funValues.size() );
                
            // Linear combination yields the result
            ValueType result = base::constantVector<FIELD::DegreeOfFreedom::size>( 0. );
            for ( unsigned f = 0; f < funValues.size(); f ++ ) 
                result += funValues[f] * doFValues[f];

            // Store result
            nodalValues[ nodeIDs[n] ] = result;
                
        } // end loop over secondary nodes

    } // end loop over elements

    return;
} 
#endif
