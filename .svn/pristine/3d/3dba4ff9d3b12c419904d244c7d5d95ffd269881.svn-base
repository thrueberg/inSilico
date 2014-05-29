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
#include <boost/bind.hpp>
// base  includes
#include <base/shape.hpp>
#include <base/geometry.hpp>
#include <base/LagrangeShapeFun.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
// base/cut includes
#include <base/cut/Cell.hpp>
#include <base/cut/SimplexMesh.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        //----------------------------------------------------------------------
        namespace detail_{

            // Helper to extract surface or volume elements (in/out)
            template<typename CELL, bool ISSURF> struct GetElements;
            
            template<typename CELL>
            struct GetElements<CELL,true>
            {
                typedef typename CELL::SurfIndexSimplex Simplex;
                
                static void apply( const CELL& cell, const bool dummy,
                                   std::vector<Simplex>& elements )
                {
                    cell.getSurface( elements );
                }
            };
                
            template<typename CELL>
            struct GetElements<CELL,false>
            {
                typedef typename CELL::VolIndexSimplex Simplex;
                
                static void apply( const CELL& cell, const bool inside,
                                   std::vector<Simplex>& elements )
                {
                    if ( inside )
                        cell.getVolumeIn(  elements );
                    else
                        cell.getVolumeOut( elements );
                }
            };

            //------------------------------------------------------------------
            template<typename MESH, unsigned SHAPEDIM, typename CELL>
            struct ExtractSimplexMesh;
            
        } // detail_

        //! Mesh of linear surface simplices
        template<typename MESH,unsigned DEGREE=1>
        struct SurfaceSimplexMesh
            : public SimplexMesh<MESH::Node::dim,MESH::Element::dim-1,DEGREE>
        { };

        //! Mesh of linear volume simplices
        template<typename MESH,unsigned DEGREE=1>
        struct VolumeSimplexMesh
            : public SimplexMesh<MESH::Node::dim,MESH::Element::dim,DEGREE>
        { };

    }
}

//------------------------------------------------------------------------------
/** Helper Object to extract a simplex mesh from mesh and cut-cell structure
 *  \tparam MESH      Volume mesh in which a surface is immersed
 *  \tparam SHAPEDIM  Dimension of the mesh to extract (volume or surface)
 *  \tparam CELL      Type of cell to extract from
*/
template<typename MESH,unsigned SHAPEDIM,typename CELL>
struct base::cut::detail_::ExtractSimplexMesh
{
    //! @name Template parameter
    //@{
    typedef MESH Mesh;
    static const unsigned shapeDim = SHAPEDIM;
    typedef CELL Cell;
    //@}

    //! Type of simplex mesh
    typedef base::cut::SimplexMesh<Mesh::Node::dim,shapeDim,Cell::geomDegree>
    SimplexMesh;

    
    typedef typename SimplexMesh::VecDim VecDim;
    typedef typename SimplexMesh::Simplex Simplex;

    //! If a surface mesh is asked for
    static const bool isSurface = (MESH::Element::dim == SHAPEDIM+1);

    //! Apply to bulk mesh and cut-cells
    static void apply( const MESH& mesh,
                       const std::vector<Cell>& cutCells,
                       SimplexMesh& simplexMesh,
                       const bool inside = true ) 
    {
        // Collect nodes and connectivity
        std::vector<VecDim>  globalNodes;
        std::vector<Simplex> globalElements;

        // begin and end of cut cells
        typename std::vector<Cell>::const_iterator cIter = cutCells.begin();
        typename std::vector<Cell>::const_iterator cEnd  = cutCells.end();

        // Access to geometry of the domain mesh needed
        typename MESH::ElementPtrConstIter eIter = mesh.elementsBegin();
        std::size_t nodeCtr = 0;

        // go through cut cells
        for ( ; cIter != cEnd; ++cIter, ++eIter ) {

            if ( cIter -> isCut() ){

                // get nodes
                std::vector<typename MESH::Element::GeomFun::VecDim> localNodes;
                cIter -> getNodes( localNodes);

                // pass all nodes to mesh
                for ( std::size_t n = 0; n < localNodes.size(); n++ ) {

                    // compute physical coordinate of local node
                    const VecDim x =
                        base::Geometry<typename MESH::Element>()( *eIter,
                                                                  localNodes[n] );

                    globalNodes.push_back( x );
                }

                // get elements (surface or volume)
                typedef detail_::GetElements<Cell,isSurface> GetElements;
                std::vector<Simplex> connec;
                GetElements::apply( *cIter, inside, connec );

                // go through these elements
                for ( std::size_t e = 0; e < connec.size(); e++ ) {

                    Simplex globalElement;
                    for ( unsigned n = 0; n < connec[e].size(); n++ ) {

                        // index of that node
                        const std::size_t nodeID = nodeCtr + connec[e][n];

                        globalElement[n] = static_cast<unsigned>( nodeID );
                    }

                    // store element
                    globalElements.push_back( globalElement );
                    
                }
                
                // increment node counter by number of nodes per cut-cell
                nodeCtr += localNodes.size();
                
            }
            
        }

        // generate mesh from nodes and connectivity
        simplexMesh.create( globalNodes, globalElements );
        
        return;

    }
    
};

//------------------------------------------------------------------------------
// Convenience functions
namespace base{
    namespace cut{

        //! Extract the surface simplex mesh from the mesh and cut-cells
        template<typename MESH, typename CELL>
        void extractSurfaceMeshFromCutCells(
            const MESH& mesh,
            const std::vector<CELL>& cutCells,
            base::cut::SurfaceSimplexMesh<MESH,CELL::geomDegree>& simplexMesh )
        {
            detail_::ExtractSimplexMesh<MESH,MESH::Element::dim-1,CELL>::apply(
                mesh, cutCells, simplexMesh );
        }

        //! Extract a volume simplex mesh from the mesh and cut-cells
        template<typename MESH, typename CELL>
        void extractVolumeMeshFromCutCells(
            const MESH& mesh,
            const std::vector<CELL>& cutCells,
            base::cut::VolumeSimplexMesh<MESH,CELL::geomDegree>& simplexMesh,
            const bool inside )
        {
            detail_::ExtractSimplexMesh<MESH,MESH::Element::dim,CELL>::apply(
                mesh, cutCells, simplexMesh, inside );
        }

    }
}

#endif
