//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   dof/generate.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_generate_hpp
#define base_dof_generate_hpp

//------------------------------------------------------------------------------
// std includes
#include <iterator>
#include <vector>
// boost includes
#include <boost/type_traits.hpp>

// base/dof includes
#include <base/dof/IndexMap.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename FEBASIS, typename MESH, typename FIELD>
        void generate( const MESH& mesh,
                       FIELD& field );

        template<typename FEBASIS, typename MESH, typename FIELD>
        void generateAndPreserveIndexMap( const MESH& mesh, FIELD& field,
                                          base::dof::IndexMap<FEBASIS>& indexMap );

    }
}

//------------------------------------------------------------------------------
/** Generate degrees of freedom topology using a given mesh.
 *  \tparam FEBASIS Type of the finite element basis 
 *  \tparam MESH    Type of geometry mesh
 *  \tparam FIELD   Type of field discretised with DoFs
 *  \param[in]   mesh   Mesh defining the problem's geometry
 *  \param[out]  field  Field to be discretised
 */
template<typename FEBASIS, typename MESH, typename FIELD>
void base::dof::generate( const MESH& mesh, FIELD& field )
{
    // Temporary index map object
    base::dof::IndexMap<FEBASIS> indexMap;
    // call function with index map
    base::dof::generateAndPreserveIndexMap( mesh, field, indexMap );
}


//------------------------------------------------------------------------------
/** Generate degrees of freedom topology using a given mesh.
 *  \tparam FEBASIS Type of the finite element basis 
 *  \tparam MESH    Type of geometry mesh
 *  \tparam FIELD   Type of field discretised with DoFs 
 *  \param[in]     mesh     Mesh defining the problem's geometry
 *  \param[out]    field    Field to be discretised
 *  \param[out]    indexMap DoF connectivity
 */
template<typename FEBASIS, typename MESH, typename FIELD>
void base::dof::generateAndPreserveIndexMap( const MESH& mesh,
                                             FIELD& field, 
                                             base::dof::IndexMap<FEBASIS>& indexMap )

{
    // Create the index map from elements
    indexMap.generateDoFIndices( mesh );

    // Number of DoFs
    field.addDoFs( indexMap.numDoFs() );
    
    // A DoF element for every element
    const std::size_t numMeshElements =
        std::distance( mesh.elementsBegin(), mesh.elementsEnd() );
    field.addElements( numMeshElements );

    // go through all elements of the mesh
    typename MESH::ElementPtrConstIter elementIter = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter lastElement = mesh.elementsEnd();
    typename FIELD::ElementPtrIter doFElementIter =
        field.elementsBegin();
    for( ; elementIter != lastElement; ++elementIter, ++doFElementIter ) {
        
        // look up the dof IDs of this element
        std::vector<std::size_t> elementDoFIDs;
        indexMap.lookUpElementDoFIndices( std::back_inserter( elementDoFIDs ),
                                          (*elementIter) -> getID() );

        // go  through elements node ptrs
        typename FIELD::Element::DoFPtrIter doF = (*doFElementIter) -> doFsBegin();

        for ( unsigned d = 0; d < elementDoFIDs.size(); d++, ++doF ) {
            // pass pointer to dof object with same ID to the dof element
            (*doF) = field.doFPtr( elementDoFIDs[d] );
        }
        
        // pass ID of mesh to dof element
        (*doFElementIter) -> setID( (*elementIter) -> getID() );
    }

    return;
}


#endif
