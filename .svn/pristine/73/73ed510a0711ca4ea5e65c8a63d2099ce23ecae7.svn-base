//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   findLocation.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_post_findlocation_hpp
#define base_post_findlocation_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <algorithm>
// base includes
#include <base/linearAlgebra.hpp>
#include <base/geometry.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace post{

        template<typename MESH>
        std::pair<typename MESH::Element*,
                  typename MESH::Element::GeomFun::VecDim>
        findLocationInMesh( const MESH& mesh,
                            const typename MESH::Node::VecDim& x,
                            const double tolerance, const unsigned maxIter );

        template<typename ELEMENT>
        std::pair<typename ELEMENT::GeomFun::VecDim,bool>
        findLocationInElement( const ELEMENT* ep,
                               const typename ELEMENT::Node::VecDim& x,
                               const double tolerance, const unsigned maxIter );

        //----------------------------------------------------------------------
        namespace detail_{
            
            //! Helper to sort the vector of pairs of distance and element pointer
            template<typename ELEMENT>
            struct CompareDistancePairs
                : public boost::function<bool (const std::pair<double,ELEMENT*>&,
                                               const std::pair<double,ELEMENT*>& )>
            {
                bool operator()( const std::pair<double,ELEMENT*>& a,
                                 const std::pair<double,ELEMENT*>& b ) const
                {
                    return (a.first < b.first);
                }
            };

        } // end detail_
    }
}

//------------------------------------------------------------------------------
/** Find a given physical coordinate in the (unstructured!) mesh.
 *  The task is to determine the pair of element \f$ e \f$ and local coordinate
 *  \f$ \xi \f$ such that 
 *  \f[
 *       | x - x_e( \xi ) | < \epsilon
 *  \f]
 *  That is the geometry representation of element \f$ e \f$ evaluated at the
 *  local coordinate \f$ \xi \f$ is within a given tolerance \f$ \epsilon \f$
 *  close to the given coordinate \f$ x \f$.
 *  
 *  In order to avoid to thoroughly check every element of the mesh, the task
 *  is prepared by
 *
 *  1.  Collection of pairs (d,e) in a vector, where \b d refers to the distance
 *      of the centroid of an element to the given coordinate \f$ x \f$ and \b e
 *      is the pointer to that element
 *  2.  Sorting of that vector according to the values of \b d
 *
 *  After this preconditioning of the problem, the sorted vector is checked
 *  beginning with the element whose centroid is closest to the point \f$ x \f$.
 *  The individual elements are then handled by the function
 *  findLocationInElement which carries out a Newton method to solve the
 *  distance problem. 
 *
 *  \note If the point \f$ x \f$ is not part of the mesh \e all elements have
 *        to be handled by the Newton method and this will become expensive.
 *
 *  \tparam MESH  Type of mesh in which a coordinate is sought
 *  \param[in] mesh      Reference to the mesh
 *  \param[in] x         Physical coordinate to find
 *  \param[in] tolerance Tolerance for coordinate comparison in Newton method
 *  \param[in] maxIter   Maximal Newton iterations 
 */
template<typename MESH>
std::pair<typename MESH::Element*,
          typename MESH::Element::GeomFun::VecDim>
base::post::findLocationInMesh( const MESH& mesh,
                    const typename MESH::Node::VecDim& x,
                    const double tolerance, const unsigned maxIter )
{
    typedef typename MESH::Element Element;
    
    typedef typename base::GeomTraits<Element> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;

    // collect elements with the distance to centroid
    typedef std::vector< std::pair<double,Element*> > DistanceElementPairVector;
    DistanceElementPairVector distanceElementPairs;
    typename MESH::ElementPtrConstIter eIter = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter eLast = mesh.elementsEnd();
    for ( ; eIter != eLast; ++eIter ) {

        const GlobalVecDim centroid =
            base::Geometry<Element>()( *eIter,
                                       base::ShapeCentroid<Element::shape>::apply() );

        const double distance = (centroid - x).norm();

        distanceElementPairs.push_back( std::make_pair( distance, *eIter ) );
    }

    // sort this vector
    detail_::CompareDistancePairs<Element> cdp;
    std::sort( distanceElementPairs.begin(), distanceElementPairs.end(), cdp );

    // go from beginning to end through this vector
    typename DistanceElementPairVector::iterator fIter = distanceElementPairs.begin();
    typename DistanceElementPairVector::iterator fLast = distanceElementPairs.end();
    for ( ; fIter != fLast; ++fIter ) {

        // check this element
        const std::pair<LocalVecDim,bool> trial =
            base::post::findLocationInElement( fIter -> second, x,
                                               tolerance, maxIter );

        // in case of success return element pointer and local coordinate
        if ( trial.second ) {
            return std::make_pair( fIter -> second, trial.first );
        }
    }

    // if this point is reached, the point has not been found in the entire mesh
    return std::make_pair( static_cast<Element*>(NULL),
                           base::invalidVector<GT::localDim>() );
}

//------------------------------------------------------------------------------
/** Find the local coordinate of an element which correpsonds to a given
 *  physical coordinate.
 *  Given an element \f$ e \f$ of a mesh and a physical coordinage \f$ x \f$,
 *  the following nonlinear problem has to be solved
 *  \f[
 *       x_e(\xi) - x = 0
 *  \f]
 *  where \f$ x_e(\xi) \f$ is the geometry representation of the considered
 *  element in function of the local coordinate \f$ \xi \f$. Using a Newton
 *  method, each iteration reads
 *  \f[
 *     \sum_{\alpha} G_\alpha(\xi^n) \Delta \xi^\alpha = x - x_e(\xi^n)
 *  \f]
 *  with the local tangent vectors \f$ G_\alpha \f$ evaluated at the current
 *  iterate \f$ \xi^n \f$. \f$ \Delta \xi^\alpha \f$ is the \f$ \alpha \f$-th
 *  component of the increment \f$ \Delta \xi \f$ and the update rule is
 *  \f[
 *      \xi^{n+1} = \xi^n + \Delta \xi
 *  \f]
 *  The initial guess of the iterate \f$ \xi^0 \f$ will be the element's
 *  centroid. A prescribed number of iterations will be carried out, unless
 *  -  the residual \f$ |x - x_e(\xi^n)| \f$ is below a given tolerance
 *  -  after the second iteration the residual increases in size
 *  -  the new iterate \f$ \xi^{n+1} \f$ would be outside of the element
 *
 *  \tparam ELEMENT      Type of element for geometry representation
 *  \param[in] ep        Pointer to element
 *  \param[in] x         Physical coordinate which is sought
 *  \param[in] tolerance Prescribed tolerance for coordinate comparison
 *  \param[in] maxIter   Maximal number of iterations
 *  \return              Pair of local coordinate and success flag
 */
template<typename ELEMENT>
std::pair<typename ELEMENT::GeomFun::VecDim,bool>
base::post::findLocationInElement( const ELEMENT* ep,
                                   const typename ELEMENT::Node::VecDim& x,
                                   const double tolerance,
                                   const unsigned maxIter )
{
    typedef typename base::GeomTraits<ELEMENT> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;

    // initial guess
    LocalVecDim xi = base::ShapeCentroid<ELEMENT::shape>::apply();

    // flag if point has been found
    bool success = false;

    // store previous residual to check if it increases
    double prevResidual = std::numeric_limits<double>::max();

    // Newton iteration
    for ( unsigned iter = 0; iter < maxIter; iter++ ) {

        // Right hand side
        const GlobalVecDim rhs = x - base::Geometry<ELEMENT>()( ep, xi );

        // already close enough, quit
        const double residual = rhs.norm();
        if ( residual < tolerance )
        {
            success = true;
            break;
        }

        // after two iterations, the residual shall not grow, otherwise quit
        if ( ( iter > 1 ) and ( residual > prevResidual ) )
        {
            success = false;
            break;
        }

        // store new residual
        prevResidual = residual;

        // Get element Jacobi matrix
        typename base::ContraVariantBasis<ELEMENT>::MatDimLDim G;
        base::ContraVariantBasis<ELEMENT>()( ep, xi, G );

        // solve system
        LocalVecDim dXi;
        dXi.noalias() = G.transpose() * rhs;

        // check if inside the reference element, if not quit
        if ( not base::InsideShape<ELEMENT::shape>::apply( xi + dXi ) )
        {
            success = false;
            break;
        }

        // update the coordinate
        xi += dXi;

    }

    // return the latest local coordinate and a success flag
    return std::make_pair( xi, success );
}

#endif
