//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Size.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_size_hpp
#define base_mesh_size_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/shape.hpp>
#include <base/geometry.hpp>
#include <base/mesh/ElementFaces.hpp>
//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename ELEMENT>
        struct Size;

        template<typename ELEMENT>
        struct MinimalEdgeLength;

        //----------------------------------------------------------------------
        // Convenience functions
        template<typename ELEMENT>
        double elementSize( const ELEMENT* ep ) {
            return Size<ELEMENT>::apply( ep );
        }

        template<typename ELEMENT>
        double minimalEdgeLength( const ELEMENT* ep ) {
            return MinimalEdgeLength<ELEMENT>::apply( ep );
        }

    }
}

//------------------------------------------------------------------------------
/** Compute an estimate of the size measure of an element.
 *  Ideally, a proper mesh size measure would be of the form
 *  \f[
 *        h_\tau = \max_{x_1, x_2 \in \tau} \| x_1 - x_2 \|
 *  \f]
 *  where \f$ \tau \f$ refers to the element under consideration. Since a true
 *  computation of this value is practically impossible for non-linear geometry
 *  representations, a very good approximation would be based on using just
 *  the element's vertices for the choices of \f$ x_1 \f$ and \f$ x_2 \f$ in
 *  above formula. Unfortunately, for 3D elements this yields to a large number
 *  of comparisons. Moreover, B-Spline based geometry representations do  not
 *  have nodal coordinates and it would be required for every \f$ x_i \f$ to
 *  evaluate the geometry representation at some \f$ \xi_i \f$.
 *  For these reasons, it is preferred to make use of the <em> aspect ratio
 *  </em> as given in PM Knupp, SIAM J SCI COMPUT 23, 2001
 *  \f[
 *       \Delta_i = \sqrt{\lambda_{ii}} = \sqrt{J_i \cdot J_i}
 *  \f]
 *  where \f$ \lambda_{ii} \f$ is a digonal entry of the metric tensor
 *  \f$ \Lambda = J^T J \f$, and \f$ J_i \f$ refers to the \f$ i \f$-th
 *  column of the Jacobi matrix \f$ J \f$.
 *  As pointed out in the reference, these quantities depend on the location
 *  of evaluation, as \f$ J = J(\xi) \f$ is a function of the local coordinate
 *  \f$ \xi \f$. In order to further reduce the effort, the centroid
 *  \f$ \xi^c \f$ of the reference shape is used.
 *
 *  Using this notation, the estimated shape size measure is given as
 *  \f[
 *       \max_{1 \leq i \leq d} \Delta_i( \xi^c )
 *  \f]
 *  with the manifold dimentions \f$ d \f$.
 *  \tparam ELEMENT Type of element to consider.
 */
template<typename ELEMENT>
struct base::mesh::Size
{
    static const unsigned localDim = base::GeomTraits<ELEMENT>::localDim;
    
    static double apply( const ELEMENT* ep )
    {
        // evaluation point is the elements centroid
        const typename base::GeomTraits<ELEMENT>::LocalVecDim centroid =
            base::ShapeCentroid<ELEMENT::shape>::apply();

        // get the elements jacobi matrix
        const typename base::JacobiMatrix<ELEMENT>::result_type J =
            base::JacobiMatrix<ELEMENT>()( ep, centroid );

        // array of length metrics
        boost::array<double,localDim> lengths;
        for ( unsigned d = 0; d < localDim; d++ )
            lengths[d] = base::norm( J.col(d) );

        return *( std::max_element( lengths.begin(), lengths.end() ) );
    }

};

//------------------------------------------------------------------------------
/** Return the minimal edge length of an element.
 *  By going through all edges of the element the distances between begin and
 *  end vertex of the edge are compared and the shortest is returned.
 *  \tparam ELEMENT Type of element to analyse.
 */
template<typename ELEMENT>
struct base::mesh::MinimalEdgeLength
{
    //! Face incidence
    typedef base::mesh::ElementFaces<ELEMENT::shape,base::EDGE> EdgeList;
    
    
    static double apply( const ELEMENT* ep )
    {
        // initialise return value with nonsense value
        double minEdgeLength = base::invalidNumber();
        // go through all edges of the element
        for ( unsigned e = 0; e < EdgeList::numFaces; e++ ) {
            // get local vertex IDs such that [v1,v2] is an edge
            const unsigned v1 = EdgeList::index( e, 0 );
            const unsigned v2 = EdgeList::index( e, 1 );
            // get node pointers
            typename ELEMENT::Node* n1 = ep -> nodePtr( v1 );
            typename ELEMENT::Node* n2 = ep -> nodePtr( v2 );

            // get coordinates of the vertices
            typename ELEMENT::Node::VecDim x1, x2;
            n1 -> getX( &(x1[0]) );
            n2 -> getX( &(x2[0]) );

            // squared distance between vertices is comparison criterion
            const double candidate = (x1-x2).squaredNorm();

            // update the storage of the minimal length
            if ( candidate < minEdgeLength ) minEdgeLength = candidate;
        }

        // return the square-root of the shortest squared distance
        return std::sqrt( minEdgeLength );
    }
};
#endif
