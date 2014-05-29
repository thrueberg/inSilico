//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SimplexMesh.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_cut_simplexmesh_hpp
#define base_cut_simplexmesh_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
// boost includes
#include <boost/array.hpp>
// base  includes
#include <base/linearAlgebra.hpp>
#include <base/Unstructured.hpp>


//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        template<unsigned DIM,unsigned SHAPEDIM=DIM,unsigned DEGREE=1>
        class SimplexMesh;
    }
}

//------------------------------------------------------------------------------
/** A simplex Mesh object.
 *  Given a list of nodes and of simplex connectivities, create a mesh
 *  object of type base::Unstructured
 *  \tparam DIM      Spatial dimension
 *  \tparam SHAPEDIM Manifold dimension
 *  \tparam DEGREE   Polynomial degree 
 */
template<unsigned DIM,unsigned SHAPEDIM,unsigned DEGREE>
class base::cut::SimplexMesh
    : public base::Unstructured<base::SimplexShape<SHAPEDIM>::value,DEGREE,DIM>
{
public:
    //! Type of the mesh from which this object inherits
    typedef base::Unstructured<base::SimplexShape<SHAPEDIM>::value,DEGREE,DIM> Base;
    //! Coordinate type
    typedef typename base::Vector<DIM>::Type VecDim;
    //! Storage of linear simplex connectivity
    typedef boost::array<unsigned,Base::Element::GeomFun::numFun>  Simplex;

    //!
    SimplexMesh() { }

    //! Constructor creates the mesh object
    SimplexMesh( const std::vector<VecDim>& nodes,
                 const std::vector<Simplex>& elements )
    {
        this -> create( nodes, elements );
    }

    void create( const std::vector<VecDim>& nodes,
                 const std::vector<Simplex>& elements )
    {
        // allocate space for nodes and elements
        this -> allocate( nodes.size(), elements.size() );

        // copy node coordinates and set node IDs
        typename Base::NodePtrIter nIter = Base::nodesBegin();
        typename Base::NodePtrIter nEnd  = Base::nodesEnd();
        for ( std::size_t n = 0; nIter != nEnd; ++nIter, n++ ) {
            (*nIter) -> setX( &(nodes[n][0]) );
            (*nIter) -> setID( n );
        }

        // Copy element connectivites and set element IDs
        typename Base::ElementPtrIter eIter = Base::elementsBegin();
        typename Base::ElementPtrIter eEnd  = Base::elementsEnd();
        for ( std::size_t e = 0; eIter != eEnd; ++eIter, e++ ) {
            (*eIter) -> setID( e );
            typename Base::Element::NodePtrIter enIter =
                (*eIter) -> nodesBegin();
            typename Base::Element::NodePtrIter enEnd =
                (*eIter) -> nodesEnd();
            for ( unsigned v = 0; enIter != enEnd; v++, ++enIter )
            {
                const unsigned n = elements[e][v];
                nIter = this->nodesBegin();
                std::advance( nIter, n );
                (*enIter) = (*nIter);
            }
        }
    }
            
};

#endif
