//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   copyConnectivity.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_copyconnectivity_hpp
#define base_dof_copyconnectivity_hpp

//------------------------------------------------------------------------------
// std   includes
#include <algorithm>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{
        namespace detail_{

            //------------------------------------------------------------------
            /** In the case that DoF connectivity structure is identical with
             *  the mesh (isoparametric case), a DoF object corresponds to
             *  every node of the mesh.
             *
             *  \tparam MESH           Type of mesh
             *  \tparam DOFINDEXARRAY  Type of array (2D) of elemnet dof indices
             *  \param[in]      mesh          Mesh from wich DoFs are generated
             *  \param[in,out]  doFIncesArray Array with the DoF numbers
             *  \return         numDoFs       Total number of unique DoFs
             */
            template<typename MESH, typename DOFINDEXARRAY>
            std::size_t copyConnectivity( const MESH& mesh,
                                          DOFINDEXARRAY& doFIndexArray )
            {
                // return number of detected dofs
                std::size_t numDoFs = 0;
                
                // deduce element type
                typedef typename MESH::Element Element;

                typename MESH::ElementPtrConstIter  eIter = mesh.elementsBegin();
                typename MESH::ElementPtrConstIter  last  = mesh.elementsEnd();

                for ( ; eIter != last; ++eIter ) {

                    // element ID
                    const std::size_t elemID = (*eIter) -> getID();

                    // go through nodes of elements and copy their IDs
                    typename Element::NodePtrConstIter nIter = (*eIter) -> nodesBegin();
                    typename Element::NodePtrConstIter nLast = (*eIter) -> nodesEnd();
                    for ( unsigned n = 0; nIter != nLast; ++nIter, n++ ) {
                        const std::size_t nodeID = (*nIter) -> getID();
                        doFIndexArray[ elemID ][ n ] = nodeID;
                        numDoFs = std::max( numDoFs, (nodeID+1) );
                    }
                }

                return numDoFs;
            }


        } // namespace detail_
    }
}
#endif
