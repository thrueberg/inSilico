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
         *  coordinates with value \f$  1 / d! \f$, \f$ d \f$ being the spatial
         *  dimension.
         *  \tparam DIM Spatial dimension
         */
        template<unsigned DIM>
        struct SimplexCentroid
        {
            static typename base::Vector<DIM,double>::Type apply()
            {
                const double value =
                    1. / static_cast<double>( base::Factorial<DIM>::value );
                return base::constantVector<DIM>( value );
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

    //--------------------------------------------------------------------------
    /** Centroid of a shape. \tparam SHAPE  The shape
     */
    template<base::Shape SHAPE>
    struct ShapeCentroid
        :base::IfElse< SHAPE == SimplexShape<ShapeDim<SHAPE>::value>::value,
                       detail_::SimplexCentroid<   ShapeDim<SHAPE>::value>,
                       detail_::HyperCubeCentroid< ShapeDim<SHAPE>::value> >::Type
    {};


}

#endif
