//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   generateDoFIndicesFromFaces.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_generatedofindicesstructured_hpp
#define base_dof_generatedofindicesstructured_hpp

//------------------------------------------------------------------------------
// std   includes
#include <algorithm>
// boost includes
#include <boost/array.hpp>
// base includs
#include <base/types.hpp>
#include <base/MultiIndex.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename FEBASIS, unsigned DIM, typename DOFINDEXARRAY>
        std::size_t generateDoFIndicesStructured( DOFINDEXARRAY& doFIndexArray,
                                                  const typename
                                                  base::MultiIndex<DIM>::Type&
                                                  gridDimensions );
    }
}

//------------------------------------------------------------------------------
/** Generate DoF index map for a structured grid.
 *  Generate for every element an array of DoF indices corresponding to the
 *  chosen Finite Element basis. Let the grid have \f$ N \f$ elements per
 *  Cartesian direction. Note that \f$ N \f$ is a multi-index and has for
 *  example in \f$ d = 3 \f$ the entries
 *  \f[
 *         N = ( n_1, n_2, n_3 )
 *  \f]
 *  The totalitiy of degrees of freedom can be represented by the multi-index
 *  \f[
 *         D = (p-c) N + c
 *  \f]
 *  where \f$ p \f$ the polynomial degree of the basis functions, \f$ c \f$ the
 *  inter-element continuity and both are interpreted as \f$ d \f$-dimensional
 *  multi-indices with constant entries. The array of supported DoFs on the
 *  element with multi-index \f$ E \f$ is denoted by \f$ [ I_1, I_2 \f$. The
 *  lower and upper multi-indices of this array are given by
 *  \f[
 *      I_1 = (p-c) E  \qquad I_2 = I_1 + p
 *  \f]
 *  In this function, this formulae are implemented by converting the linear
 *  element and local DoF indices to multi-indices, calculating the global
 *  index values, and converting them back to linear indices.
 *
 *  \tparam FEBASIS       Finite Element basis
 *  \tparam DIM           Spatial dimension
 *  \tparam DOFINDEXARRAY Array of DoF indices
 *  \param[out] doFIndexArray  Array of DoF indices (filled by this function)
 *  \param[in]  gridDimensions Dimensions of the grid (i.e. \f$ N \f$ above)
 *  
 */
template<typename FEBASIS, unsigned DIM, typename DOFINDEXARRAY>
std::size_t
base::dof::generateDoFIndicesStructured( DOFINDEXARRAY& doFIndexArray,
                                         const typename
                                         base::MultiIndex<DIM>::Type& gridDimensions )
{
    // Multi-index types
    typedef typename base::MultiIndex<DIM> MultiIndex;
    typedef typename MultiIndex::Type      MultiIndexType;

    // Introspect the FE basis
    static const unsigned degree      = FEBASIS::degree;
    static const unsigned continuity  = FEBASIS::continuity;
    static const unsigned numElemDoFs = FEBASIS::FiniteElement::numTotalDoFs;

    // Total DoF array dimensions
    const MultiIndexType doFArrayDim =
        (degree-continuity) * gridDimensions + (continuity + 1);

    // number of dofs per element per direction
    const MultiIndexType elemDoFDim = MultiIndex::constant( degree+1 );

    // total number of elements
    const std::size_t numElements = MultiIndex::length( gridDimensions );

    // verify the first dimension of the index array
    VERIFY_MSG( doFIndexArray.size() == numElements,
                "Wrong dimensions of index array" );
    
    // go through all elements of the grid
    for ( std::size_t e = 0; e < numElements; e++ ) {

        // generate element multi-index
        const MultiIndexType eM = MultiIndex::wrap( e, gridDimensions );

        // multi-index of begin of the element DoF array
        const MultiIndexType startM = (degree - continuity) * eM;

        // go through the element's DoFs
        for ( unsigned d = 0; d < numElemDoFs; d++ ) {

            // Local array multi-index
            const MultiIndexType localM = MultiIndex::wrap( d, elemDoFDim );

            // global multi-index of the current DoF
            const MultiIndexType globalM = startM + localM;

            // global linear index of the current DoF
            const std::size_t global = MultiIndex::unwrap( globalM, doFArrayDim );

            // set index array entry
            doFIndexArray[e][d] = global;
        }
        
    }

    // length of multi-index gives total number of DoFs
    return MultiIndex::length( doFArrayDim );
}


#endif
