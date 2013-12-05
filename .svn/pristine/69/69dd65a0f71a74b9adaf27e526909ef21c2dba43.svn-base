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

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename ELEMENT>
        struct Size;
    }
}

//------------------------------------------------------------------------------
/** Compute an estimate of the size measure of an element.
 *  Ideally, a proper mesh size measure would be of the form
 *  \f[
 *        h_\tau = \max_{x_1, x_2 \in \tau} \| x_1 - x_2 \|
 *  \f]
 *  where \f$ \tau \f$ refers to the element under consideration. Since a true
 *  compuation of this value is practically impossible for non-linear geometry
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
 *  \f$ \xi^c \f$of the reference shape is used.
 *
 *  Using this notation, the estimate shape size measure is given as
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

#endif
