//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   distanceToElement.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_distanceToElement_hpp
#define base_cut_distanceToElement_hpp
//------------------------------------------------------------------------------
// boost inclues
#include <boost/array.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
#include <base/fe/LagrangeElement.hpp>
#include <base/cut/LevelSet.hpp>
#include <base/auxi/compareNumbers.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename SELEMENT>
        void distanceToElement( const SELEMENT* surfEp,
                                const bool isSigned, 
                                base::cut::LevelSet<SELEMENT::Node::dim>& ls,
                                const double pointIdentityTolerance );

        //----------------------------------------------------------------------
        template<typename VECDIM>
        VECDIM closestPointToLine( const VECDIM& x1, const VECDIM& x2,
                                   const VECDIM& x,
                                   base::Vector<1,double>::Type& xi );

        base::Vector<3,double>::Type
        closestPointToTriangle( const base::Vector<3,double>::Type& x1,
                                const base::Vector<3,double>::Type& x2,
                                const base::Vector<3,double>::Type& x3,
                                const base::Vector<3,double>::Type& x,
                                base::Vector<2,double>::Type& xi );
        
        base::Vector<3,double>::Type
        closestPointToQuadrilateral( const base::Vector<3,double>::Type& x1,
                                     const base::Vector<3,double>::Type& x2,
                                     const base::Vector<3,double>::Type& x3,
                                     const base::Vector<3,double>::Type& x4,
                                     const base::Vector<3,double>::Type& x,
                                     base::Vector<2,double>::Type& xi   );

        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            template<unsigned DIM, base::Shape SHAPE> struct ComputeClosestPoint;

            template<unsigned DIM>
            struct ComputeClosestPoint<DIM,base::LINE>
            {
                typedef typename base::Vector<DIM,double>::Type   VecDim;
                typedef typename base::Vector<1,double>::Type     VecLDim;
                    
                template<typename ARRAY>
                static VecDim apply( const ARRAY& array, const VecDim& x, VecLDim& xi )
                {
                    return base::cut::closestPointToLine( array[0], array[1], x, xi );
                }
            };

            template<unsigned DIM>
            struct ComputeClosestPoint<DIM,base::TRI>
            {
                typedef typename base::Vector<DIM,double>::Type VecDim;
                typedef typename base::Vector<2,double>::Type   VecLDim;
                    
                template<typename ARRAY>
                static VecDim apply( const ARRAY& array, const VecDim& x, VecLDim& xi )
                {
                    return base::cut::closestPointToTriangle( array[0], array[1],
                                                              array[2], x, xi );
                }
            };

            template<unsigned DIM>
            struct ComputeClosestPoint<DIM,base::QUAD>
            {
                typedef typename base::Vector<DIM,double>::Type VecDim;
                typedef typename base::Vector<2,double>::Type   VecLDim;
                    
                template<typename ARRAY>
                static VecDim apply( const ARRAY& array, const VecDim& x, VecLDim& xi )
                {
                    return base::cut::closestPointToQuadrilateral( array[0], array[1],
                                                                   array[2], array[3],
                                                                   x, xi );
                }
            };

            //------------------------------------------------------------------
            template<unsigned DIM> struct SignedDistanceToPlane;
            template<>             struct SignedDistanceToPlane<2>;
            template<>             struct SignedDistanceToPlane<3>;
            
        } // namespace detail_

        
    }
}

//------------------------------------------------------------------------------
/** Compute the distance between a point and an element.
 *  A given surface element is sampled by a Lagrangian element with the same
 *  shape and the same degree of geometry representation. Then, in function of
 *  the shape, the closest point is computed based on a linear geometry
 *  approximation (using the vertices of the Lagrangian element only).
 *  If the new distance is closer than the previously stored one, the passed
 *  level set datum is updated, otherwise untouched.
 *  \tparam SELEMENT  Type of surface element (not necessarily interpolatory)
 *  \param[in]     surfEp   Pointer to surface element
 *  \param[in]     isSigned Flag if distance function shall be signed
 *  \param[in,out] ls       Old and new level set datum
 *  \param[in]     pointIdentityTolerance Tolerance for point comparison
 */
template<typename SELEMENT>
void base::cut::distanceToElement( const SELEMENT* surfEp,
                                   const bool isSigned, 
                                   base::cut::LevelSet<SELEMENT::Node::dim>& ls,
                                   const double pointIdentityTolerance )
{
    typedef base::cut::LevelSet<SELEMENT::Node::dim> LevelSet;
        
    // type of an equivalent lagrange element
    typedef base::fe::LagrangeElement<SELEMENT::shape,
                                      SELEMENT::GeomFun::degree> LagrangeElement;

    // get interpolation support points of the lagrange element
    static const unsigned numPoints = LagrangeElement::ShapeFun::numFun;
    boost::array<typename LagrangeElement::ShapeFun::VecDim,
                 numPoints> supportPoints;
    LagrangeElement::ShapeFun::supportPoints( supportPoints );

    // compute the corresponding global coordinates
    boost::array<typename SELEMENT::Node::VecDim, numPoints> globalX;
    for ( unsigned p = 0; p < numPoints; p++ ) {

        globalX[p] =
            base::Geometry<SELEMENT>()( surfEp, supportPoints[p] );
    }

    // compute the closest point
    const typename LevelSet::VecDim x = ls.getX();
    typename base::Vector<SELEMENT::dim>::Type xi;
    const typename LevelSet::VecDim newClosestPoint =
        detail_::ComputeClosestPoint<LevelSet::dim,
                                     SELEMENT::shape>::apply( globalX, x, xi );
                                                                        
    // distances
    const double oldDistance = ls.getUnsignedDistance();
    const double newDistance = base::norm(x - newClosestPoint);

    // if new distances is shorter, update the level set data
    if ( ( not isSigned ) and ( newDistance < oldDistance ) ) {

        // set closest point and element
        ls.setClosestPoint(           newClosestPoint );
        ls.setClosestElement(         surfEp -> getID() );
        ls.setClosestLocalCoordinate( xi );

    }
    else if ( isSigned ) {

        // if closest points are numerically equal, decide by distance to plane
        const bool decideByDistancesToPlane =
            base::auxi::almostEqualVectors<SELEMENT::Node::dim>(
                ls.getClosestPoint(), newClosestPoint,
                pointIdentityTolerance );

        if ( decideByDistancesToPlane ) {

            // get old and new signed spanned volumina
            const double oldDistanceToPlane = ls.getDistanceToPlane();
            const double newDistanceToPlane =
                detail_::SignedDistanceToPlane<LevelSet::dim>::apply( globalX, x );

            // if the distance to the plane of the element is larger, update
            if ( std::abs( newDistanceToPlane ) > std::abs( oldDistanceToPlane ) ) {

                // set new level set data
                ls.setClosestPoint(           newClosestPoint );
                ls.setClosestElement(         surfEp -> getID() );
                ls.setClosestLocalCoordinate( xi );
                ls.setDistanceToPlane(        newDistanceToPlane );

            }
            // else leave old datum unchanged
        }
        else if ( newDistance < oldDistance ) {

            const double newDistanceToPlane =
                detail_::SignedDistanceToPlane<LevelSet::dim>::apply( globalX, x );

            // set new level set data
            ls.setClosestPoint(           newClosestPoint );
            ls.setClosestElement(         surfEp -> getID() );
            ls.setClosestLocalCoordinate( xi );
            ls.setDistanceToPlane(        newDistanceToPlane );

        }

    }
    
    
    return;
}

//------------------------------------------------------------------------------
/** Compute the point on a line element closest to a given point.
 *  With the line element's parameterisation, any point in space can be
 *  expressed as 
 *  \f[
 *         x(\xi,\eta) = a + (b - a) \xi + n \eta
 *  \f]
 *  with \f$ n \f$ a vector orthogonal to the line element. Now a dot-product
 *  of this representation with a the tangent \f$ t = x_2 - x_1 \f$ gives
 *  in case of a give point \f$ x \f$
 *  \f[
 *         (t, t) \xi = (p - a, t)
 *  \f]
 *  The code makes use of Ericson's (Real-time collision detection, 5.1.2)
 *  optimised version which defers the division until really necessary.
 *  The code is copied from the book and adapted to this code's conventions.
 *  \tparam VECDIM Type of vector (should be 2D or 3D)
 *  \param[in] a, b  Begin and end points of the line segment
 *  \param[in] p     Point to which the closest point is sought
 *  \param[in] xi    Local coordinate of the closest point
 *  \return          Cloeset point on the segment (a,b) to p
 */
template<typename VECDIM>
VECDIM base::cut::closestPointToLine( const VECDIM& a, const VECDIM& b,
                                      const VECDIM& p,
                                      typename base::Vector<1,double>::Type& xi)
{
    // tangent vector
    const VECDIM t = b - a;
    
    // Project c onto ab, but deferring divide by Dot(ab, ab)
    const VECDIM pa = p - a;
    double s = base::dotProduct( pa, t);

    if (s <= 0.0) {
        // c projects outside the [a,b] interval, on the a side; clamp to a
        // t = 0.0;
        xi[0] = 0.0;
        return a;
    }
    else {
        const double denom = base::dotProduct(t, t);
        // Always nonnegative since denom = |ab|^2
        if (s >= denom) {
            // c projects outside the [a,b] interval, on the b side; clamp to b
            //t = 1.0;
            xi[0] = 1.0;
            return b;
        }
        else {
            // c projects inside the [a,b] interval; must do deferred divide now
            s /= denom;
        }
    }
    
    xi[0] = s;
    
    return a + s * t;
    
}

//------------------------------------------------------------------------------
/** Compute the point in a triangle closest to a given point.
 *  Given a point \f$ p \f$ and a triangle spanned by \f$ a \f$, \f$ b \f$ and
 *  \f$ c \f$, this point can be represented as
 *  \f[
 *      p = p(s,t) = a + (b-a) s + (c-a) t + n w
 *                 = a (1-s-t) + b s + c t + n w
 *  \f]
 *  with the bary-centric coordiantes \f$ u = 1 - s - t \f$, \f$ s \f$ and
 *  \f$ t \f$. \f$ s \f$ and \f$ t \f$ can be computed by taking the dot-
 *  products of this representation with the tangent vectors \f$ b-a \f$ and
 *  \f$ c-a \f$ respectively. But if the resulting projection
 *  \f[
 *      p^\prime = a + (b-a) s + (c-a) t
 *  \f]
 *  is not inside the triangle, vertices and edges have to be tested in order
 *  to find the true point inside the triangle closest to \f$ p \f$.
 *  The following code does this in an optimised fashion and is copied and
 *  adapted from Ericson's book Real-time collision detection, chapter 5.1.5.
 *  \param[in] a, b, c  The vertices of the triangle
 *  \param[in] p        The point whose closest point is sought
 *  \param[in] xi    Local coordinate of the closest point
 *  \return             Closest point inside the triangle to p
 */
base::Vector<3,double>::Type
base::cut::closestPointToTriangle( const base::Vector<3,double>::Type& a,
                                   const base::Vector<3,double>::Type& b,
                                   const base::Vector<3,double>::Type& c,
                                   const base::Vector<3,double>::Type& p,
                                   base::Vector<2,double>::Type& xi  )
{
    typedef base::Vector<3,double>::Type Vec3;
    
    // Check if P in vertex region outside A
    const Vec3 ab = b - a;
    const Vec3 ac = c - a;
    const Vec3 ap = p - a;
    const double d1 = base::dotProduct(ab, ap);
    const double d2 = base::dotProduct(ac, ap);
    if ( (d1 <= 0.0) and (d2 <= 0.0) ) {
        xi[0] = xi[1] = 0.0;
        return a; // barycentric coordinates (1,0,0)
    }
    
    // Check if P in vertex region outside B
    const Vec3 bp = p - b;
    const double d3 = base::dotProduct(ab, bp);
    const double d4 = base::dotProduct(ac, bp);
    if ( (d3 >= 0.0) and (d4 <= d3) ) {
        xi[0] = 1.0; xi[1] = 0.0;
        return b; // barycentric coordinates (0,1,0)
    }
    
    // Check if P in edge region of AB, if so return projection of P onto AB
    const double vc = d1 * d4 - d3 * d2;
    if ( (vc <= 0.0) and (d1 >= 0.0) and (d3 <= 0.0)) {
        const double v = d1 / (d1 - d3);
        xi[0] = v; xi[1] = 0.0;
        return a + v * ab; // barycentric coordinates (1-v,v,0)
    }
    
    // Check if P in vertex region outside C
    const Vec3   cp = p - c;
    const double d5 = base::dotProduct(ab, cp);
    const double d6 = base::dotProduct(ac, cp);
    if ( (d6 >= 0.0) and (d5 <= d6)) {
        xi[0] = 0.0; xi[1] = 1.0;
        return c; // barycentric coordinates (0,0,1)
    }

    // Check if P in edge region of AC, if so return projection of P onto AC
    const double vb = d5 * d2 - d1 * d6;
    if ( (vb <= 0.0) and (d2 >= 0.0) and (d6 <= 0.0) ) {
        const double w = d2 / (d2 - d6);
        xi[0] = 0.; xi[1] = w;
        return a + w * ac; // barycentric coordinates (1-w,0,w)
    }
    
    // Check if P in edge region of BC, if so return projection of P onto BC
    const double va = d3 * d6 - d5 * d4;
    if ( (va <= 0.0) and ((d4 - d3) >= 0.0) and ((d5 - d6) >= 0.0)) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        xi[0] = 1.0 - w; xi[1] = w;
        return b + w * (c - b); // barycentric coordinates (0,1-w,w)
    }
    
    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    xi[0] = v; xi[1] = w;
    return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0-v-w
}

//------------------------------------------------------------------------------
/** Compute the point closest inside a quadrilateral to a given point.
 *  For simplicity the given quadrilateral is subdivided into two triangles,
 *  the closest points w.r.t. each triangle are computed and the closer one
 *  is returned.
 *  \param[in] x1, x2, x3, x4  The four vertices of the quadrilateral
 *  \param[in] x               Point of interest
 *  \param[in] xi    Local coordinate of the closest point
 *  \return                    Cloesest point to p inside the quadrilateral
 */
base::Vector<3,double>::Type
base::cut::closestPointToQuadrilateral( const base::Vector<3,double>::Type& x1,
                                        const base::Vector<3,double>::Type& x2,
                                        const base::Vector<3,double>::Type& x3,
                                        const base::Vector<3,double>::Type& x4,
                                        const base::Vector<3,double>::Type& x,
                                        base::Vector<2,double>::Type& xi   )
{
    typedef base::Vector<3,double>::Type Vec3;

    base::Vector<2,double>::Type xi1, xi2;

    // First triangle
    const Vec3 cp1 = closestPointToTriangle( x1, x2, x3, x, xi1 );
    // Second triangle
    const Vec3 cp2 = closestPointToTriangle( x1, x3, x4, x, xi2 );

    // distances squared
    const double dist1 = base::dotProduct( cp1 - x, cp1 - x );
    const double dist2 = base::dotProduct( cp2 - x, cp2 - x );

    // return the closer point
    if ( dist1 < dist2 ) {
        xi = xi1;
        return cp1;
    }

    xi = xi2;
    return cp2;
    
}

//------------------------------------------------------------------------------
template<>
struct base::cut::detail_::SignedDistanceToPlane<2>
{
    /** Compute the signed distance to the plane (here line) that embeds the
     *  considered element. This distance is used in case of ambiguities when
     *  the closest point is shared by many elements and the 'correct' closest
     *  element has to be chose.
     *  Consider the picture
     *  \code{.txt}
     *       .......(A)_____________(B)..............
     *                      |                 |
     *                      | n               |
     *                      V                (P)
     *  \endcode
     *  where a point \f$ P\f$ is to be tested against the line segment
     *  \f$ [A,B] \f$. The vertical line between \f$ p \f$ and the infinite
     *  extension of \f$ [A,B] \f$ is of interest. The signed distance between
     *  \f$ P \f$ and that line is given by
     *  \f[
     *       d = \frac{ n \cdot (PA) }{ |n| }
     *  \f]
     *  where \f$ n \f$ is the non-normalised (!) outward normal vector to the
     *  hyperplane and \f$ PA \f$ the vector from \f$ P \f$ to \f$ A \f$.
     *  In this picture the result will be negative since the vectors \f$ n \f$
     *  and \f$ PA \f$ point in opposite directions.
     *  
     *  \tparam ARRAY Type of arry with element vertices
     *  \param[in] array Element's vertices
     *  \param[in] P     Point to check
     *  \return          Signed distance to line
     */
    template<typename ARRAY>
    static double apply( const ARRAY& array,
                         const base::Vector<2,double>::Type& P )
    {
        // construct normal vector 
        base::Vector<2,double>::Type n;
        n[0] =  array[1][1] - array[0][1];
        n[1] = -array[1][0] + array[0][0];

        // vector between point and begin of line segment
        const base::Vector<2,double>::Type PA = array[0] - P;
        
        const double numer  = base::dotProduct( n, PA );
        const double denom2 = base::dotProduct( n, n  );
        return numer / std::sqrt( denom2 );
    }
};

//------------------------------------------------------------------------------
template<>
struct base::cut::detail_::SignedDistanceToPlane<3>
{
    /** Compute the signed distance between given point and plane.
     *  The distance with sign between the given point and the plane which
     *  embeds a given triangle is of interest. Let \f$ [A,B,C] \f$ denote the
     *  oriented triangle and \f$ P \f$ the point of interest. Then
     *  \f[
     *      d = \frac{ n \cdot PA }{ |n| }
     *  \f]
     *  gives the signed distance based on the outward normal vector \f$ n \f$
     *  (not normaalised) and the vector \f$ PA \f$ pointing from \f$ P \f$ to
     *  \f$ A \f$.
     *  \param[in] array Element's vertices
     *  \param[in] P     Point to check
     *  \return          Signed distance     
     */
    template<typename ARRAY>
    static double apply( const ARRAY& array,
                         const base::Vector<3,double>::Type& P )
    {
        // normal vector
        const base::Vector<3,double>::Type n =
            base::crossProduct( array[1] - array[0], array[2] - array[0] );

        const base::Vector<3,double>::Type PA = array[0] - P;
        
        const double numer  = base::dotProduct( n, PA );
        const double denom2 = base::dotProduct( n, n  );
        return numer / std::sqrt( denom2 );
    }
    
};


            


#endif
