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
            // Bind a vector of DoF-size to the type of value
            template<typename FIELD>
            struct ValueTypeBinder
            {
                static const unsigned size = FIELD::DegreeOfFreedom::size;
                typedef typename base::Vector<size,number>::Type Type;
            };

        }

        //-----------------------------------------------------------------------
        template<typename VALUETYPE, typename MESH, typename FIELD, typename OP>
        void evaluateAtNodes( const MESH& mesh, const FIELD& field,
                              const OP& doFOperator,
                              std::vector<VALUETYPE>& nodalValues );

        //-----------------------------------------------------------------------
        /** Simply evaluate the field itself at some point in the past.
         */
        template<unsigned HIST,typename MESH, typename FIELD>
        void evaluateFieldHistoryAtNodes(
            const MESH& mesh, const FIELD& field,
            std::vector<typename detail_::ValueTypeBinder<FIELD>::Type>& nodalValues )
        {
            // Functor for the evaluation of the field
            boost::function<base::number( const typename FIELD::DegreeOfFreedom*,
                                          const unsigned)> getDoFHistoryComponent =
                boost::bind( &FIELD::DegreeOfFreedom::template getHistoryValue<HIST>,
                             _1, _2 );
            
            typedef typename detail_::ValueTypeBinder<FIELD>::Type ValueType;
            evaluateAtNodes<ValueType>( mesh, field, getDoFHistoryComponent,
                                        nodalValues );
        }

        //-----------------------------------------------------------------------
        /** Evaluate the field itself at the current time point.
         */
        template<typename MESH, typename FIELD>
        void evaluateFieldAtNodes(
            const MESH& mesh, const FIELD& field,
            std::vector<typename detail_::ValueTypeBinder<FIELD>::Type>& nodalValues )
        {
            // delegate call for HIST=0
            evaluateFieldHistoryAtNodes<0>( mesh, field, nodalValues );
        }

    }
}

//------------------------------------------------------------------------------
/** Apply a functor to the DoFs and interpolate the result at the nodes.
 *  Since nodes and DoF locations do not necessarily coincide, the post-
 *  processing requires to evaluate the field at the locations of the nodes of
 *  the mesh. This function does that job by going over all elements of field
 *  and mesh, and interpolating the field at the local coordinates that are
 *  the nodes of the geometry representation. Moreover, a functor is provided
 *  by the caller of type
 *  \code{.cpp}
 *       ValueComponent doFOperator( const DoFPtr, const unsigned )
 *  \endcode
 *  by means of which the some information from the components of the DoFs
 *  is retrieved.
 *  Commonly, this functor is just the access to the values of the field but
 *  could be something different too.
 *
 *  \tparam VALUETYPE  Type of the result of the field evaluation
 *  \tparam MESH       Type of geometry representation
 *  \tparam FIELD      Type of field to evaluate
 *  \tparam OP         Type of evaluation functor
 *  \param[in]  mesh        Geometry representation
 *  \param[in]  field       The field to evaluate
 *  \param[in]  doFOperator Evaluation functor 
 *  \param[out] nodalValues Results at the nodes of the mesh
 */
template<typename VALUETYPE, typename MESH, typename FIELD, typename OP>
void base::post::evaluateAtNodes( const MESH& mesh, const FIELD& field,
                                  const OP& doFOperator,
                                  std::vector<VALUETYPE>& nodalValues )
{
    // The type of the result values
    typedef VALUETYPE ValueType;

    // resize the result container
    const std::size_t numNodes =
        static_cast<std::size_t>( std::distance( mesh.nodesBegin(), mesh.nodesEnd() ) );
    nodalValues.resize( numNodes  );

    // Mesh- and Field-Element iterators
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

            // apply functor to get desired value
            ValueType doFValue;
            for ( unsigned d = 0; d < FIELD::DegreeOfFreedom::size; d++ ) 
                doFValue[d] = doFOperator( *dIter, d );

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
            ValueType result = funValues[0] * doFValues[0];
                //base::constantVector<FIELD::DegreeOfFreedom::size>( 0. );
            for ( unsigned f = 1; f < funValues.size(); f ++ ) 
                result += funValues[f] * doFValues[f];

            // Store result
            nodalValues[ nodeIDs[n] ] = result;
                
        } // end loop over secondary nodes

    } // end loop over elements

    return;
} 
#endif
