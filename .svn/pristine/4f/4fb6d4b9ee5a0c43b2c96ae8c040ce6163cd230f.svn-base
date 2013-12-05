//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   extractMeshFromCutCells.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_extractmeshfromcutcells_hpp
#define base_cut_extractmeshfromcutcells_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// boost includes
#include <boost/array.hpp>
#include <boost/bind.hpp>
// base  includes
#include <base/shape.hpp>
#include <base/geometry.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
// base/cut includes
#include <base/cut/Cell.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        //! Type definitions for a linear simplex mesh
        template<typename MESH,unsigned DIM>
        struct SimplexMesh
        {
            static const base::Shape simplexShape =
                base::SimplexShape<DIM>::value;
            
            typedef typename MESH::Node                    Node;
            typedef base::LagrangeShapeFun<1,simplexShape> SFun;
            typedef base::mesh::Element<Node,SFun>         Element;
            typedef base::mesh::Unstructured<Element>      Type;
        };

        //! Mesh of linear surface simplices
        template<typename MESH>
        struct SurfaceSimplexMesh
            : public SimplexMesh<MESH,MESH::Element::dim-1>
        { };

        //! Mesh of linear volume simplices
        template<typename MESH>
        struct VolumeSimplexMesh
            : public SimplexMesh<MESH,MESH::Element::dim>
        { };

        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            //! Count specific item of the cut cells, i.e. nodes or elements
            template<typename CUTCELLS, typename OBJNUM>
            std::size_t countCutCellObjects( const CUTCELLS& cutCells,
                                             const OBJNUM& objNum )
            {
                std::size_t result = 0;
                typename CUTCELLS::const_iterator cIter = cutCells.begin();
                typename CUTCELLS::const_iterator cEnd  = cutCells.end();
                for ( ; cIter != cEnd; ++cIter ) {
                    if ( cIter -> isCut() ) {
                        result += objNum( *cIter );
                    }
                }
                
                return result;
            }

            //------------------------------------------------------------------
            //! Generate the nodes from a cut cell mesh
            template<typename MESH, typename CUTCELLS, typename SIMPMESH>
            void generateNodes( const MESH& mesh,
                                const CUTCELLS& cutCells,
                                SIMPMESH& simplexMesh )
            {
                // begin and end of cut cells
                typename CUTCELLS::const_iterator cIter = cutCells.begin();
                typename CUTCELLS::const_iterator cEnd  = cutCells.end();

                // Access to geometry of the domain mesh needed
                typename MESH::ElementPtrConstIter eIter = mesh.elementsBegin();

                // node iteration
                typename SIMPMESH::NodePtrIter  nIter = simplexMesh.nodesBegin();
                std::size_t nodeID = 0;
                
                for ( ; cIter != cEnd; ++cIter, ++eIter ) {

                    // get nodes
                    std::vector<typename MESH::Element::GeomFun::VecDim> nodes;
                    cIter -> getNodes( nodes);

                    // pass all nodes to mesh
                    for ( std::size_t n = 0; n < nodes.size(); n++ ) {

                        // compute physical coordinate of local node
                        const typename SIMPMESH::Node::VecDim x =
                            base::Geometry<typename MESH::Element>()( *eIter,
                                                                          nodes[n] );

                        // set coordiantes and ID
                        (*nIter) -> setX( &(x[0]) );
                        (*nIter) -> setID( nodeID );

                        // increment
                        nodeID++;
                        ++nIter;
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            template<typename MESH,unsigned DIM> struct ExtractSimplexMesh;
            
            
        } // namespace detail_

    }
}



//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        namespace detail_{

            template<typename CELL, bool ISSURF> struct GetElements;
            
            template<typename CELL>
            struct GetElements<CELL,true>
            {
                static void apply( const CELL& cell, const bool dummy,
                                   std::vector<typename CELL::SurfIndexSimplex>&
                                   elements )
                {
                    cell.getSurface( elements );
                }
            };
                
            template<typename CELL>
            struct GetElements<CELL,false>
            {
                static void apply( const CELL& cell, const bool inside,
                                   std::vector<typename CELL::VolIndexSimplex>&
                                   elements )
                {
                    if ( inside )
                        cell.getVolumeIn(  elements );
                    else
                        cell.getVolumeOut( elements );
                }
            };

            
        }
    }
}

//------------------------------------------------------------------------------
template<typename MESH, unsigned DIM>
struct base::cut::detail_::ExtractSimplexMesh
{
    typedef typename SimplexMesh<MESH,DIM>::Type  SimplexMesh;
    typedef base::cut::Cell<MESH::Element::shape> Cell;

    static const bool isSurface = (MESH::Element::dim == DIM+1);

    static void apply( const MESH& mesh,
                       const std::vector<Cell>& cutCells,
                       SimplexMesh& simplexMesh,
                       const bool inside = true ) 
    {
    
        // number of nodes and elements in the extracted mesh
        const std::size_t numNodes =
            detail_::countCutCellObjects( cutCells,
                                          boost::bind( &Cell::numNodes, _1 ) );
        
        const std::size_t numElements =
            ( isSurface ? 
              detail_::countCutCellObjects( cutCells,
                                            boost::bind( &Cell::numSurfaceElements, _1 ) ) :
              ( inside ?
                detail_::countCutCellObjects( cutCells,
                                              boost::bind( &Cell::numVolumeInElements, _1 ) ) :
                detail_::countCutCellObjects( cutCells,
                                              boost::bind( &Cell::numVolumeOutElements, _1 ) ) ) );
                      
        // allocate the mesh of the surface simplices
        simplexMesh.allocate( numNodes, numElements );
        
        // create nodes of the mesh
        detail_::generateNodes( mesh, cutCells, simplexMesh );
    
        // iterate over cut cells
        typedef std::vector<Cell>                     CellVector;
        typename CellVector::const_iterator cIter = cutCells.begin();
        typename CellVector::const_iterator cEnd  = cutCells.end();

        // generate elements from cut cells
        typename SimplexMesh::ElementPtrIter eIter = simplexMesh.elementsBegin();
        std::size_t elemID = 0;
        std::size_t nodeID = 0;
        for ( ; cIter != cEnd; ++cIter ) {

            if ( cIter -> isCut() ) {
                
                // get surface elements
                std::vector<typename base::cut::Simplex<DIM,unsigned>::Type > connec;
                //cIter -> getSurface( connec );
                detail_::GetElements<Cell,isSurface>::apply( *cIter, inside, connec );
                
                for ( std::size_t e = 0; e < connec.size(); e++ ) {
                    
                    // iterator for the nodes of the surface element
                    typename SimplexMesh::Element::NodePtrIter nIter =
                        (*eIter) -> nodesBegin();
                    
                    for ( unsigned n = 0; n < connec[e].size(); n++, ++nIter ) {
                        // index of that node
                        const std::size_t nodeIndex = nodeID + connec[e][n];
                        // set node pointer
                        (*nIter) = simplexMesh.nodePtr( nodeIndex );
                    }
                    
                    // set ID of element
                    (*eIter) -> setID( elemID );

                    // increment
                    elemID++;
                    ++eIter;
                }
                // increment node counter by number of nodes per cut-cell
                nodeID += cIter -> numNodes();
            }
        }

        return;
    }
    
};

//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        
        template<typename MESH>
        void extractSurfaceMeshFromCutCells(
            const MESH& mesh,
            const std::vector<base::cut::Cell<MESH::Element::shape> >& cutCells,
            typename SurfaceSimplexMesh<MESH>::Type& simplexMesh )
        {
            detail_::ExtractSimplexMesh<MESH,MESH::Element::dim-1>::apply(
                mesh, cutCells, simplexMesh );
        }

        template<typename MESH>
        void extractVolumeMeshFromCutCells(
            const MESH& mesh,
            const std::vector<base::cut::Cell<MESH::Element::shape> >& cutCells,
            typename VolumeSimplexMesh<MESH>::Type& simplexMesh,
            const bool inside )
        {
            detail_::ExtractSimplexMesh<MESH,MESH::Element::dim>::apply(
                mesh, cutCells, simplexMesh, inside );
        }

    }
}

#endif
