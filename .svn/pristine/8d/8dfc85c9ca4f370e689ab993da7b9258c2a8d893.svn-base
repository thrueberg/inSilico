//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   generateSurfaceMesh.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_generatesurfacemesh_hpp
#define base_cut_generatesurfacemesh_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// base  includes
#include <base/linearAlgebra.hpp>
// base/cut includes
#include <base/cut/Cell.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{


        template<typename DOMAINMESH>
        struct SurfaceMeshBinder
        {
            typedef typename DOMAINMESH::Element DomainElement;
            static const base::Shape surfSimplex =
                base::SimplexShape<DomainElement::dim-1>::value;
            typedef base::mesh::SurfaceElement<DomainElement,
                                               surfSimplex>   SurfaceElement;
            typedef base::mesh::Unstructured<SurfaceElement>  SurfaceMesh;
        };


        template<typename DOMAINMESH, typename CELL>
        void generateSurfaceMesh(
            const DOMAINMESH& domainMesh,
            const std::vector<CELL> & cutCells,
            typename SurfaceMeshBinder<DOMAINMESH>::SurfaceMesh& surfaceMesh );
    }
}

//------------------------------------------------------------------------------
template<typename DOMAINMESH, typename CELL>
void base::cut::generateSurfaceMesh(
    const DOMAINMESH& domainMesh,
    const std::vector<CELL> & cutCells,
    typename base::cut::SurfaceMeshBinder<DOMAINMESH>::SurfaceMesh& surfaceMesh )
{
    typedef typename DOMAINMESH::Element DomainElement;
    typedef typename base::cut::SurfaceMeshBinder<DOMAINMESH>::SurfaceMesh SurfaceMesh;
    typedef typename SurfaceMesh::Element SurfaceElement;
            
    // count number of surface elements and nodes
    std::size_t numSurfElements = 0;
    std::size_t numSurfNodes    = 0;
    const unsigned numNodesPerElement = CELL::dim;
    for ( std::size_t c = 0; c < cutCells.size(); c++ ) {
        const std::size_t tmp = cutCells[c].numSurfaceElements();
        numSurfElements += tmp;
        numSurfNodes    += tmp * numNodesPerElement;
    }

    // prepare mesh
    surfaceMesh.allocate( numSurfNodes, numSurfElements );

    std::size_t nodeId = 0;

    // Go through all boundary elements and new surface elements
    typename SurfaceMesh::ElementPtrIter surfElemIter =
        surfaceMesh.elementsBegin();
    for ( std::size_t c = 0; c < cutCells.size(); c++ ) {

        if ( cutCells[c].isCut() ) {

            // Get relevant data from iterator
            DomainElement* domainElement = domainMesh.elementPtr( c );

            //
            std::vector<typename CELL::VecDim> nodes;
            cutCells[c].getNodes( nodes );

            //
            std::vector<typename CELL::SurfIndexSimplex> surfSimplices;
            cutCells[c].getSurface( surfSimplices );

            for ( unsigned s = 0; s < surfSimplices.size();
                  s++, ++surfElemIter ) {

                // Pass pointer to domain element to surface element
                (*surfElemIter) -> setDomainElementPointer( domainElement );

                // Access to surface elements geometry nodes
                typename SurfaceElement::NodePtrIter nodePtrIter =
                    (*surfElemIter) -> nodesBegin();

                // Go through the parameter points of the surface element
                typename SurfaceElement::ParamIter paramIter =
                    (*surfElemIter) -> parametricBegin();
                typename SurfaceElement::ParamIter paramEnd  =
                    (*surfElemIter) -> parametricEnd();

                for ( unsigned p = 0;
                      paramIter != paramEnd;
                      ++paramIter, ++nodePtrIter, p++ ) {

                    // Get local coordinate 
                    const typename CELL::VecDim xi = nodes[ surfSimplices[s][p] ];

                    // Assign to parameter of surface element
                    *paramIter = xi;

                    const typename CELL::VecDim x =
                        base::Geometry<DomainElement>()( domainElement, xi );

                    // get pointer to surface mesh node
                    typename SurfaceMesh::Node* np = surfaceMesh.nodePtr( nodeId );
                    // set node data
                    np -> setID( nodeId );
                    np -> setX( &(x[0]) );
                    // pass pointer to surface element
                    *nodePtrIter = np;

                    // pass node pointer from old to new mesh
                    nodeId++;
                                                
                } // end loop over the simplex nodes
            } // end loop over cut-cells surface simplices
            
        } // if cell is cut
            
    } // end loop over surface elements


            
}
        
#endif
