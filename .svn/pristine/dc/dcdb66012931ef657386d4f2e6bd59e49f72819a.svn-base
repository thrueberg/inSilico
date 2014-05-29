//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   shape.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_shape_hpp
#define base_shape_hpp

//------------------------------------------------------------------------------
// std includes
#include <string>
// base includes
#include <base/meta.hpp>
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{

    //--------------------------------------------------------------------------
    //! Enumerator for description of the geometry of an element.
    enum Shape {
        POINT = 0,       //!< point
        LINE,            //!< line element
        TRI,             //!< triangle element
        QUAD,            //!< quadrilateral element
        TET,             //!< tetrehedral element
        HEX              //!< hexahedral element
    };

    //--------------------------------------------------------------------------
    /** Enumerator for description of a shape's faces.
     *  http://en.wikipedia.org/wiki/Face_%28geometry%29
     *  \note The numeric values of this enum are used in the code.
     *        Therefore, they <em> must not </em> be changed.
     */
    enum NFace {
        VERTEX = 0, //!< Vertex = 0-face
        EDGE,       //!< Edge   = 1-face
        FACE,       //!< Face   = 2-face
        CELL        //!< Cell   = 3-face
    };


    //--------------------------------------------------------------------------
    /** \class base::ShapeDim
     * The manifold dimension a shape implies.
     * Every base::Shape implies a spatial dimension, e.g. base::TRI is
     * two-dimensional. This is not necessarily the spatial dimension of the
     * considered physical problem. Consider, e.g., shells for an example where
     * the shape is two- but the problem is three-dimensional.
     * \tparam SHAPE  The shape which implies a dimension
     */
    template<Shape SHAPE> struct ShapeDim;

    //! \cond SKIPDOX
    template<> struct ShapeDim<POINT>{ static const unsigned value = 0; };
    template<> struct ShapeDim<LINE>{  static const unsigned value = 1; };
    template<> struct ShapeDim<TRI>{   static const unsigned value = 2; };
    template<> struct ShapeDim<QUAD>{  static const unsigned value = 2; };
    template<> struct ShapeDim<TET>{   static const unsigned value = 3; };
    template<> struct ShapeDim<HEX>{   static const unsigned value = 3; };
    //! \endcond 
    
    //--------------------------------------------------------------------------
    /** \class base::ShapeName
     *  A human-readable name used for I/O of shapes.
     *  Associates with every base::Shape a unique string which can be used
     *  as a descriptor in in- and out-put methods.
     *  \tparam SHAPE The shape which is paired with a string
     */
    template<Shape SHAPE> struct ShapeName;

    //! \cond SKIPDOX
    template<> struct ShapeName<POINT>
    { static std::string apply() { return std::string( "point" ); } };

    template<> struct ShapeName<LINE>
    { static std::string apply() { return std::string( "line" ); } };

    template<> struct ShapeName<TRI>
    { static std::string apply() { return std::string( "triangle" ); } };

    template<> struct ShapeName<QUAD>
    { static std::string apply() { return std::string( "quadrilateral" ); } };

    template<> struct ShapeName<TET>
    { static std::string apply() { return std::string( "tetrahedron" ); } };

    template<> struct ShapeName<HEX>
    { static std::string apply() { return std::string( "hexahedron" ); } };
    //! \endcond

    //--------------------------------------------------------------------------
    /** \class base::HyperCubeShape
     *  Define the shape which forms a hyper-cube in the given dimension.
     *  \sa http://en.wikipedia.org/wiki/Hypercube
     *  \tparam DIM Spatial dimension for which a hyper-cube is of interest
     */
    template<unsigned DIM> struct HyperCubeShape;

    //! \cond SKIPDOX
    template<> struct HyperCubeShape<0> { static const Shape value = POINT; };
    template<> struct HyperCubeShape<1> { static const Shape value = LINE; };
    template<> struct HyperCubeShape<2> { static const Shape value = QUAD; };
    template<> struct HyperCubeShape<3> { static const Shape value = HEX;  };
    //! \endcond

    /** \class base::SimplexShape
     *  Define the shape which would form a simplex in the given dimension.
     *  \sa http://en.wikipedia.org/wiki/Simplex
     *  \tparam DIM Spatial dimension for which a simplex is of interest
     */
    template<unsigned DIM> struct SimplexShape;

    //! \cond SKIPDOX
    template<> struct SimplexShape<0> { static const Shape value = POINT; };
    template<> struct SimplexShape<1> { static const Shape value = LINE; };
    template<> struct SimplexShape<2> { static const Shape value = TRI; };
    template<> struct SimplexShape<3> { static const Shape value = TET;  };
    //! \endcond 

    namespace detail_{
        //----------------------------------------------------------------------
        /** Compute the number of N-faces of a DIM-Hypercube.
         *  Remember: 0-face = vertex, 1-face = edge, etc.
         *  The number of N-faces for dimension D is computed as
         *  \f[
         *      E_{D,M} = 2^{M-D} C(D,M)
         *  \f]
         *  with the binomial coefficient \f$ C(D,M) \f$.
         *  See: http://en.wikipedia.org/wiki/Hypercube
         *  \tparam DIM    Spatial dimension of the hypercube
         *  \tparam NFACE  Type of M-face to be counted
         */
        template<unsigned DIM, base::NFace NFACE>
        struct NumHyperCubeFaces
        {
            static const unsigned diff  = ( NFACE > DIM ? 0 : DIM-NFACE );
            static const unsigned value =
                ( NFACE > DIM ? 0 :
                  base::MToTheN<2,diff>::value * base::Binomial<DIM,NFACE>::value );
        };

        /** Compute the number of N-faces of a DIM-Simplex.
         *  Remember: 0-face = vertex, 1-face = edge, etc.
         *  The number of N-faces for dimension D is computed as
         *  \f[
         *      E_{D,M} = C(D+1,M+1)
         *  \f]
         *  with the binomial coefficient \f$ C(D,M) \f$.
         *  See http://en.wikipedia.org/wiki/Simplex
         *  \tparam DIM   Spatial dimension of simplex
         *  \tparam NFACE Type of M-face to be counted
         */
        template<unsigned DIM, base::NFace NFACE>
        struct NumSimplexFaces
        {
            static const unsigned aux = NFACE + 1;
            static const unsigned value =
                ( NFACE > DIM ? 0 : base::Binomial<DIM+1,aux>::value );
        };
    }

    //--------------------------------------------------------------------------
    /** Provide the number of N-faces for a given shape.
     *  \tparam SHAPE   Type of geometric shape
     *  \tparam NFACE   Type of N-face to be counted
     */
    template<base::Shape SHAPE,base::NFace NFACE>
    struct NumNFaces
        : base::IfElse< SHAPE == SimplexShape<ShapeDim<SHAPE>::value>::value,
                        detail_::NumSimplexFaces<  ShapeDim<SHAPE>::value,NFACE>,
                        detail_::NumHyperCubeFaces<ShapeDim<SHAPE>::value,NFACE> >::Type
    {
        // empty
    };

    //--------------------------------------------------------------------------
    /** Deduce the shape of the (dim-codim)-face of a shape.
     *  Here the notion of co-dimension is used
     *  (http://en.wikipedia.org/wiki/Codimension), which represents the
     *  difference in spatial dimensions between the shape and its face.
     *  As examples, the co-dimension of an EDGE of a QUAD is 1 (=2-1) and the
     *  the co-dimension of an EDGE of a TET is 2 (=3-1).
     *  \tparam SHAPE  Geometric shape of which the surface shape is queried
     *  \tparam CODIM  Co-dimension of the face
     */
    template<base::Shape SHAPE, unsigned CODIM=1>
    struct FaceShape
        : base::IfElse< SHAPE == SimplexShape<ShapeDim<SHAPE>::value>::value,
                        SimplexShape<  ShapeDim<SHAPE>::value-CODIM>,
                        HyperCubeShape<ShapeDim<SHAPE>::value-CODIM> >::Type
    {
        STATIC_ASSERT_MSG( (CODIM <= ShapeDim<SHAPE>::value),
                           "Co-dimension larger than shape dimension");
    };

    //--------------------------------------------------------------------------
    /** Define the NFace which is the surface of a given shape.
     *  E.g. for a HEX we would have a FACE, for a TRI an EDGE, etc.
     *  The action of this meta-function could be put directly into the
     *  template defintions but the required casting makes it more legible
     *  this way.
     *  \tparam SHAPE  Shape of which the surface is sought.
     */
    template<base::Shape SHAPE>
    struct ShapeSurface
    {
        STATIC_ASSERT_MSG( (SHAPE != POINT), "You are out of your senses!" );
        
        static const base::NFace value =
            static_cast<base::NFace>( base::ShapeDim<SHAPE>::value - 1 );
    };

    //--------------------------------------------------------------------------
    namespace detail_{
        /** Centroid of a simplex.
         *  The controid of any simplex is located at a point with equal
         *  coordinates with value \f$ 1 / f(d) \f$, \f$ d \f$ being the spatial
         *  dimension and \f$ f(d) = 2, 3, 6 \f$, respectively.
         *  \tparam DIM Spatial dimension
         */
        template<unsigned DIM>
        struct SimplexCentroid
        {
            static typename base::Vector<DIM,double>::Type apply()
            {
                const unsigned aux = static_cast<unsigned>( DIM==3 );
                const double denom =
                    static_cast<double>( (DIM+1) + 2 * aux );
                return base::constantVector<DIM>( 1. / denom );
            }
        };

        /** Centroid of a hypercube.
         *  \tparam DIM Spatial dimension
         */
        template<unsigned DIM>
        struct HyperCubeCentroid
        {
            static typename base::Vector<DIM,double>::Type apply()
            {
                return base::constantVector<DIM>( 0.5 );
            }
        };
    }

    /** Centroid of a shape. \tparam SHAPE  The shape
     */
    template<base::Shape SHAPE>
    struct ShapeCentroid
        :base::IfElse< SHAPE == SimplexShape<ShapeDim<SHAPE>::value>::value,
                       detail_::SimplexCentroid<   ShapeDim<SHAPE>::value>,
                       detail_::HyperCubeCentroid< ShapeDim<SHAPE>::value> >::Type
    {};

    //--------------------------------------------------------------------------
    namespace detail_{

        /** Size of the DIM-simplex reference element is \f$ 1/ d! \f$.
         *  That is 1, 1/2, 1/6 in the dimensions 1, 2, and 3.
         *  \tparam DIM Spatial dimension of the element
         */
        template<unsigned DIM>
        struct SimplexRefSize
        {
            static double apply()
            {
                return 1. /static_cast<double>( base::Factorial<DIM>::value );
            }
        };

        /** Size of the DIM-hypercube reference element is always 1.
         *  \tparam DIM Spatial dimension of the element
         */
        template<unsigned DIM>
        struct HyperCubeRefSize { static double apply() { return 1.0; } };
    }

    /** Return the size of the reference element of the given shape.
     *  \tparam SHAPE Shape to consider
     */
    template<base::Shape SHAPE>
    struct RefSize
        : base::IfElse< SHAPE == SimplexShape<ShapeDim<SHAPE>::value>::value,
                        detail_::SimplexRefSize<  ShapeDim<SHAPE>::value>,
                        detail_::HyperCubeRefSize<ShapeDim<SHAPE>::value> >::Type
    {};


    //--------------------------------------------------------------------------
    namespace detail_{

        /** Check if coordinate is inside the reference Simplex.
         *  The bounds are that all coordinate components have to be
         *  positive
         *  \f[
         *        \xi_i \geq 0
         *  \f]
         *  and that upper bounds can be expressed as
         *  \f[
         *        \xi_i \leq 1 - \sum_{j < i} \xi_j
         *  \f]
         *  \tparam DIM Spatial dimension of the simplex
         */
        template<unsigned DIM>
        struct InsideSimplex
        {
            static bool apply( const typename base::Vector<DIM>::Type& xi,
                               const double tol )
            {
                // lower bounds
                for ( unsigned d = 0; d < DIM; d++ ) {
                    if ( xi[d] < -tol ) return false;
                }

                // upper bounds
                double rhs = 1.;
                for ( unsigned d = 0; d < DIM; d++ ) {
                    if ( xi[d] > rhs+tol ) return false;
                    rhs -= xi[d];
                }
                
                return true;
            }
        };

        /** Check if given coordinate is inside the hypercube reference domain.
         *  The test is simply
         *  \f[
         *        0 \leq \xi_i \leq 1
         *  \f]
         *  for all coordinate components.
         *  \tparam DIM Spatial dimension of the simplex
         */
        template<unsigned DIM>
        struct InsideHyperCube
        {
            static bool apply( const typename base::Vector<DIM>::Type& xi,
                               const double tol )
            {
                // lower bounds
                for ( unsigned d = 0; d < DIM; d++ ) {
                    if ( xi[d] <   -tol ) return false;
                    if ( xi[d] > 1.+tol ) return false;
                }

                return true;
            }
        };
        
    }

    /** Check if a given coordinate is inside the reference domain of a shape.
     *  \tparam SHAPE  Shape of reference domain
     */
    template<base::Shape SHAPE>
    struct InsideShape
        : base::IfElse< SHAPE == SimplexShape<ShapeDim<SHAPE>::value>::value,
                        detail_::InsideSimplex<  ShapeDim<SHAPE>::value>,
                        detail_::InsideHyperCube<ShapeDim<SHAPE>::value> >::Type
    {};


    //--------------------------------------------------------------------------
    namespace detail_{

        template<unsigned DIM>
        struct SnapToHyperCube
        {
            typedef typename base::Vector<DIM>::Type VecDim;

            static VecDim apply( const VecDim& xi )
            {
                VecDim result = xi;
                for ( unsigned d = 0; d < DIM; d++ ) {
                    if      ( xi[d] < 0. ) result[d] = 0.;
                    else if ( xi[d] > 1. ) result[d] = 1.;
                }
                
                return result;
            }
        };

        //----------------------------------------------------------------------
        //! Dummy not to be called
        template<unsigned DIM>
        struct SnapToSimplex
        {
            typedef typename base::Vector<DIM>::Type VecDim;

            static VecDim apply( const VecDim& xi )
            {
                STATIC_ASSERT_MSG( DIM == 0, "Not to be called" );
            }
        };
        
        //----------------------------------------------------------------------
        /** Find point inside standard triangle closest to the given point.
         *  The standard triangle has vertices A=(0,0), B=(1,0) and C=(0,1). A
         *  given point \f$ \xi \f$ can have the following relative locations
         *
         *  1. closest to vertices A, B, or C
         *  2. closest to edges A-B, A-C, or B-C
         *  3. inside the triangle A-B-C
         *
         *  In this object these locations are checked in the given order. In
         *  case the point falls in one of the categories, the closest point
         *  (i.e. the detected vertex, the projection onto an edge or simply
         *  the point itself) is computed and returned.
         */
        template<>
        struct SnapToSimplex<2>
        {
            typedef base::Vector<2>::Type VecDim;

            static VecDim apply( const VecDim& xi )
            {
                VecDim result; 

                //--------------------------------------------------------------
                // Vertex regions
                if ( ( xi[0] <= 0. ) and ( xi[1] <= 0. ) ) {
                    result[0] = 0.;
                    result[1] = 0.;
                    return result; // Vertex A=(0,0)
                }

                if ( ( xi[0] >= 1. ) and ( xi[0] >= 1. + xi[1] ) ) {
                    result[0] = 1.;
                    result[1] = 0.;
                    return result; // Vertex B=(1,0)
                }
                
                if ( ( xi[1] >= 1. ) and ( xi[1] >= 1. + xi[0] ) ) {
                    result[0] = 0.;
                    result[1] = 1.;
                    return result; // Vertex C=(0,1)
                }

                //--------------------------------------------------------------
                // Edge regions
                if ( ( xi[0] >= 0. ) and ( xi[0] <= 1. ) and ( xi[1] <= 0. ) ) {
                    result[0] = xi[0];
                    result[1] = 0.;
                    return result; // Edge A-B
                }

                if ( ( xi[1] >= 0. ) and ( xi[1] <= 1. ) and ( xi[0] <= 0. ) ) {
                    result[0] = 0.;
                    result[1] = xi[1];
                    return result; // Edge A-C
                }

                const double s = (1. - xi[0] + xi[1]) / 2.;
                if ( ( s >= 0. ) and ( s <= 1. ) and (xi[0] + xi[1] >= 1.) ) {
                    result[0] = 1. - s;
                    result[1] = s;
                    return result; // Edge B-C
                }

                // Point is inside the triangle
                return xi;
            }
        };

        //----------------------------------------------------------------------
        /** Find point inside standard tetrahedron closest to the given point.
         *  The standard tetrahedron has vertices A=(0,0,0), B=(1,0,0),
         *  C=(0,1,0) and D=(0,0,1). A given point \f$ \xi \f$ can have the
         *  following relative locations
         *
         *  1. closest to vertices A, B, C, or D
         *  2. closest to edges A-B, A-C, A-D, B-C, C-D, or D-B
         *  3. closest to faces A-B-C, A-B-D, A-C-D, or B-C-D
         *  4. inside the tetrahedron A-B-C-D
         *
         *  In this object these locations are checked in the given order. In
         *  case the point falls in one of the categories, the closest point
         *  (i.e. the detected vertex, the projection onto an edge or face, or
         *  simply the point itself) is computed and returned.
         */
        template<>
        struct SnapToSimplex<3>
        {
            typedef base::Vector<3>::Type VecDim;

            static VecDim apply( const VecDim& xi )
            {
                VecDim result; 

                //--------------------------------------------------------------
                // Vertex regions
                if ( ( xi[0] <= 0. ) and ( xi[1] <= 0. ) and ( xi[2] <= 0. ) ){
                    result[0] = 0.;
                    result[1] = 0.;
                    result[2] = 0.;
                    return result; // Vertex A=(0,0,0)
                }

                if ( ( xi[0] >= 1. ) and
                     ( xi[0] >= 1. + xi[1] ) and ( xi[0] >= 1. + xi[2] ) ) {
                    result[0] = 1.;
                    result[1] = 0.;
                    result[2] = 0.;
                    return result; // Vertex B=(1,0,0)
                }
                
                if ( ( xi[1] >= 1. ) and
                     ( xi[1] >= 1. + xi[0] ) and ( xi[1] >= 1. + xi[2] ) ) {
                    result[0] = 0.;
                    result[1] = 1.;
                    result[2] = 0.;
                    return result; // Vertex C=(0,1,0)
                }
                
                if ( ( xi[2] >= 1. ) and
                     ( xi[2] >= 1. + xi[0] ) and ( xi[2] >= 1. + xi[1] ) ) {
                    result[0] = 0.;
                    result[1] = 0.;
                    result[2] = 1.;
                    return result; // Vertex D=(0,0,1)
                }

                //--------------------------------------------------------------
                // Edge regions

                // edges along the main coordinate axes
                if ( ( xi[0] >= 0. ) and ( xi[0] <= 1. ) and
                     ( xi[1] <= 0. ) and ( xi[2] <= 0. ) ) {
                    result[0] = xi[0];
                    result[1] = 0.;
                    result[2] = 0.;
                    return result; // Edge A-B
                }
                
                if ( ( xi[1] >= 0. ) and ( xi[1] <= 1. ) and
                     ( xi[0] <= 0. ) and ( xi[2] <= 0. ) ) {
                    result[0] = 0.;
                    result[1] = xi[1];
                    result[2] = 0.;
                    return result; // Edge A-C
                }

                if ( ( xi[2] >= 0. ) and ( xi[2] <= 1. ) and
                     ( xi[0] <= 0. ) and ( xi[1] <= 0. ) ) {
                    result[0] = 0.;
                    result[1] = 0.;
                    result[2] = xi[2];
                    return result; // Edge A-D
                }

                // diagonal edges
                const double sum = xi[0] + xi[1] + xi[2];
                if ( sum >= 1. ) {
                
                    const double s = (1. - xi[0] + xi[1]) / 2.;
                    if ( ( s >= 0. ) and ( s <= 1. ) and ( xi[2] <= 0. ) ) {
                        result[0] = 1. - s;
                        result[1] = s;
                        result[2] = 0.;
                        return result; // Edge B-C
                    }

                    const double t = (1. - xi[1] + xi[2]) / 2.;
                    if ( ( t >= 0. ) and ( t <= 1. ) and ( xi[0] <= 0. ) ) {
                        result[0] = 0.;
                        result[1] = 1. - t;
                        result[2] = t;
                        return result; // Edge C-D
                    }
                    
                    const double u = (1. - xi[2] + xi[1]) / 2.;
                    if ( ( u >= 0. ) and ( u <= 1. ) and ( xi[1] <= 0. ) ) {
                        result[0] = u;
                        result[1] = 0.;
                        result[2] = 1. - u;
                        return result; // Edge D-B
                    }

                }

                //--------------------------------------------------------------
                // Face regions
                typedef detail_::InsideSimplex<2> CheckTri;
                typedef base::Vector<2>::Type     Vec2;

                Vec2 abc; abc[0] = xi[0]; abc[1] = xi[1];
                if ( CheckTri::apply( abc, 0. ) and ( xi[2] <= 0. ) ) {
                    result[0] = xi[0];
                    result[1] = xi[1];
                    result[2] = 0.;
                    return result; // Face A-B-C
                }
                
                Vec2 abd; abd[0] = xi[0]; abd[1] = xi[2];
                if ( CheckTri::apply( abd, 0. ) and ( xi[1] <= 0. ) ) {
                    result[0] = xi[0];
                    result[1] = 0.;
                    result[2] = xi[2];
                    return result; // Face A-B-D
                }

                Vec2 acd; acd[0] = xi[1]; acd[1] = xi[2];
                if ( CheckTri::apply( acd, 0. ) and ( xi[0] <= 0. ) ) {
                    result[0] = 0.;
                    result[1] = xi[1];
                    result[2] = xi[2];
                    return result; // Face A-C-D
                }

                const double w = sum - 1.;
                Vec2 bcd; bcd[0] = xi[1] - w/3.; bcd[1] = xi[2] - w/3.;
                if ( CheckTri::apply( bcd, 0. ) and ( w >= 0. ) ) {
                    result[0] = xi[0] - w/3.;
                    result[1] = xi[1] - w/3.;
                    result[2] = xi[2] - w/3.;
                    return result; // Face B-C-D
                }

                // reaching this point -> point is inside the TET
                return xi;
            }
        };
    }

    /**  Snap to a reference element's polytope.
     *   Given a parametric coordinate \f$ \xi \f$, these objects calculate the
     *   parametric coordinate \f$ \eta \f$ inside of the reference shape which
     *   is closest to the given \f$ \xi \f$. Let \f$ \hat{\tau} \f$ denote the
     *   reference polytope, then the task is to find \f$ \eta \f$ such that
     *   \f[
     *       \min_{\eta \in \hat{\tau}} | \eta - \xi |
     *   \f]
     *   is achieved.
     *   If the given \f$ \xi \f$ is already inside \f$ \hat{\tau} \f$ we
     *   obviously have \f$ \eta = \xi \f$. Otherwise the vertex, edge or face
     *   of \f$ \hat{\tau} \f$ has to be found to which \f$ \xi \f$ is closest
     *   and its orthogonal projection onto that entity is returned.
     *   Whereas the computations for the unit hypercube shapes is trivially
     *   the check against \f$ [0,1] \f$ in all coordinate directions, the
     *   computations for the triangle and tetrahedron are quite tedious.
     *   \tparam SHAPE  Type of shape to find closest point in
     */
    template<base::Shape SHAPE>
    struct SnapToShape
        : base::IfElse< SHAPE == HyperCubeShape<ShapeDim<SHAPE>::value>::value,
                        detail_::SnapToHyperCube<ShapeDim<SHAPE>::value>,
                        detail_::SnapToSimplex<  ShapeDim<SHAPE>::value> >::Type
    { };
}

#endif
