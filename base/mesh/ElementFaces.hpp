//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ElementFaces.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_elementfaces_hpp
#define base_mesh_elementfaces_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/array.hpp>
// baes includes
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        //----------------------------------------------------------------------
        //! Define the basic features of an element's n-Face
        template<base::Shape ELEMSHAPE, base::NFace NFACE,
                 typename IDXTYPE = unsigned>
        struct ElementFaceTraits
        {
            //! Shape of the requested face
            static const base::Shape shape =
                base::FaceShape<ELEMSHAPE,
                                base::ShapeDim<ELEMSHAPE>::value - NFACE>::value;

            //! Number of vertices on the face
            static const unsigned numVertices =
                base::NumNFaces<shape, base::VERTEX>::value;

            //! Storage of surface element connectivity
            typedef boost::array<IDXTYPE,numVertices> Type;

            //! Assign an index
            static void assignIndex( Type& face, const unsigned num,
                                     const IDXTYPE index )
            {
                face[num] = index;
            }

        };

        //! \cond SKIPDOX
        //----------------------------------------------------------------------
        //! Define the basic features of an element's n-Face
        template<base::Shape SHAPE, typename IDXTYPE>
        struct ElementFaceTraits<SHAPE,base::VERTEX,IDXTYPE>
        {
            static const base::Shape shape = base::POINT;
            
            static const unsigned numVertices = 1;
            typedef IDXTYPE Type;

            //! Assign an index
            static void assignIndex( Type& vertex, const unsigned dummy,
                                     const IDXTYPE index )
            {
                vertex = index;
            }

        };
        //! \endcond
        
        //----------------------------------------------------------------------
        template<base::Shape SHAPE, base::NFace NFACE>
        struct ElementFaces;

        //----------------------------------------------------------------------
        template<base::Shape SHAPE>
        struct FaceEdges;

    }
}

//------------------------------------------------------------------------------
/** Provides a connectivity look-up table for faces of elements.
 *  Given an array of vertex IDs of an element (i.e. its connectivity), the
 *  look-up table of this object provides the connectivities of the element's
 *  faces, i.e. edges and faces (as 2D objects). Note that a hierarchical
 *  ordering (as opposed to the lexicographic ordering of the shape functions)
 *  is assumed here. The hierarchical ordering is assumed for the input
 *  connectivity AND for the output connectivities!!
 *
 *  \note Only the functionality of edges of all shapes and faces of
 *        TET and HEX elements is considered. 
 *
 *  \param SHAPE  Type of geometric shape of an element
 *  \param NFACE  Type of face, that is of interest
 */
template<base::Shape SHAPE, base::NFace NFACE>
struct base::mesh::ElementFaces
{
    //! Number of n-faces of element
    static const unsigned numFaces = base::NumNFaces<SHAPE, NFACE>::value;

    //! Traits class for access of properties
    typedef base::mesh::ElementFaceTraits<SHAPE,NFACE> ElementFaceTraits;

    //! Type of Face
    typedef typename ElementFaceTraits::Type Face;

    //! Storage of all surface arrays
    typedef boost::array<Face, numFaces> FaceTable;

    //! Return j-th vertex index of i-th face
    static unsigned index( const unsigned i, const unsigned j )
    {
        return faceTable_[i][j];
    }

private:
    /** Lookup table:
     *   faceTable_[ i ][ j ]
     *  gives the j-th vertex index of the i-th face
     */
    static const FaceTable faceTable_;

    //------------------------------------------------------------------
    //! Sanity checks
    STATIC_ASSERT_MSG( (NFACE <= base::ShapeDim<SHAPE>::value),
                       "These cases are not implemented" );

    //! Note that the cases vertex and cell would be have to be implemented
    //! for every shape value even though they look the same. Hence,
    //! These partial specialisations are implemented at the bottom of this
    //! file.
    STATIC_ASSERT_MSG( (NFACE != base::VERTEX),
                       "Not implemented here for technical reasons" );
    //! See above.
    STATIC_ASSERT_MSG( (NFACE != base::CELL),
                       "Not implemented here for technical reasons" );
};

//==============================================================================
namespace base{
    namespace mesh{
        //----------------------------------------------------------------------
        // HERE follow the partial specialisations of ElementFaces for the cases
        // base::VERTEX and base::CELL.
        
        //! \cond SKIPDOX
        //----------------------------------------------------------------------
        // VERTEX case
        template<base::Shape SHAPE>
        struct ElementFaces<SHAPE,base::VERTEX>
        {
            //! Number of n-faces of element
            static const unsigned numFaces = base::NumNFaces<SHAPE,
                                                             base::VERTEX>::value;

            //! Type of face, just an index
            typedef unsigned Face;

            //! Storage of all surface arrays
            typedef boost::array<Face, numFaces> FaceTable;

            //! Return vertex index of i-th face=vertex
            static unsigned index( const unsigned i, const unsigned j )
            {
                return i; //!< trivial result: the i-th face is vertex no. i
            }

        };
        
        //----------------------------------------------------------------------
        // CELL case
        template<base::Shape SHAPE>
        struct ElementFaces<SHAPE,base::CELL>
        {
            //! Number of n-faces of element
            static const unsigned numFaces = 1;

            static const unsigned numVertices = base::NumNFaces<SHAPE,
                                                                base::VERTEX>::value;

            //! Type of face, just an index
            typedef boost::array<unsigned,numVertices> Face;

            //! Storage of all surface arrays
            typedef boost::array<Face, numFaces> FaceTable;

            //! Return j-th vertex index 
            static unsigned index( const unsigned i, const unsigned j )
            {
                return j; //!< Trivial result: the j-th vertex of the CELL is j
            }

        };
        //! \endcond
    }
}

//==============================================================================
// Initialisations of the static member faceTable_ for all full specialisations
namespace base{
    namespace mesh{

        namespace detail_{

            //! Generic definition, never implemented
            template<base::Shape SHAPE,base::NFace NFACE>
            typename ElementFaces<SHAPE,NFACE>::FaceTable
            initialiseFaceTable();
                
            //==================================================================
            // EDGES
            //==================================================================

            //------------------------------------------------------------------
            /**  base::LINE
             *   <pre>
             *
             *                 0--->--1
             *   </pre>
             */  
            
            template<> inline
            ElementFaces<base::LINE,base::EDGE>::FaceTable
            initialiseFaceTable<base::LINE,base::EDGE>()
            {
                ElementFaces<base::LINE,base::EDGE>::FaceTable aux = {{
                        {{ 0, 1 }}
                    }};
                return aux;
            }

            //------------------------------------------------------------------
            /**  base::TRIANGLE
             *   <pre>
             *
             *   2     2         2
             *   V     |`\        |             \
             *   |     |  `\        `           \
             *   |[1]  |    `\    [1] `         \
             *   |     |      `\        `       \
             *   V     |        `\        |     \
             *   0     0----------1         1
             *     
             *             [0]
             *         0>-------->1
             *   </pre>
             */
            template<> inline
            ElementFaces<base::TRI,base::EDGE>::FaceTable
            initialiseFaceTable<base::TRI,base::EDGE>()
            {

                ElementFaces<base::TRI,base::EDGE>::FaceTable aux = {{
                        {{ 0, 1 }}, {{ 1, 2 }}, {{ 2, 0 }}
                    }};

                return aux;
            }

            //------------------------------------------------------------------
            /** base::QUAD
             *  <pre>
             *              3<---------<2
             *                   [2] (xi_1=1)
             *
             *      3       3-----------2       2
             *      V       |           |       A
             *      |       |           |       |
             *      |[3]    |(xi_0=0)   |    [1]| (xi_0=1)
             *      |       |           |       |
             *      V       |           |       A
             *      0       0-----------1       1
             *  
             *                   [0] (xi_1=0)
             *              0>--------->1
             *  </pre>
             */
            template<> inline
            ElementFaces<base::QUAD,base::EDGE>::FaceTable
            initialiseFaceTable<base::QUAD,base::EDGE>()
            {
                ElementFaces<base::QUAD,base::EDGE>::FaceTable aux = {{
                        {{ 0, 1 }}, {{ 1, 2 }}, {{ 2, 3 }}, {{ 3, 0 }}
                    }};
                return aux;
            }

            //------------------------------------------------------------------
            /** base::TET (view on faces from outside, i.e. against positive normal)
             *  <pre>
             *
             *                3              [0]:  (0) ---- (1)
             *              ,/|`\                  
             *            ,/  |  `\          [1]:  (1) ---- (2)      
             *          ,/    '.   `\              
             *        ,/       |     `\      [2]:  (2) ---- (0)      
             *      ,/         |       `\  
             *     0-----------'.--------2 
             *      `\.         |      ,/    [3]:  (3) ---- (0)
             *         `\.      |    ,/            
             *            `\.   '. ,/        [4]:  (3) ---- (1)
             *               `\. |/                
             *                  `1           [5]:  (3) ---- (2)
             *  </pre>
             */
            template<> inline
            ElementFaces<base::TET,base::EDGE>::FaceTable
            initialiseFaceTable<base::TET,base::EDGE>()
            {
                ElementFaces<base::TET,base::EDGE>::FaceTable aux = {{
                        {{ 0, 1 }}, {{ 1, 2 }}, {{ 2, 0 }},
                        {{ 3, 0 }}, {{ 3, 1 }}, {{ 3, 2 }}
                    }};
                return aux;
            }

            //------------------------------------------------------------------
            /** base::HEX (view on faces from outside, i.e. against positive normal)
             *  <pre>
             *                         {bottom}
             *                         [0]: (0) ---- (1)     [1]: (1) ---- (2)
             *     4----------7                                               
             *     |\         |\       [2]: (2) ---- (3)     [3]: (3) ---- (0)
             *     | \        | \      
             *     |  \       |  \     {top}                                  
             *     |   5----------6    [4]: (4) ---- (5)     [5]: (5) ---- (6)
             *     |   |      |   |                                           
             *     0---|------3   |    [6]: (6) ---- (7)     [7]: (7) ---- (4)
             *      \  |       \  |    
             *       \ |        \ |    {sides}
             *        \|         \|    [8]: (0) ---- (4)     [9]: (1) ---- (5)
             *         1----------2                                           
             *  </pre>                [10]: (2) ---- (6)    [11]: (3) ---- (7)
             */
            template<> inline
            ElementFaces<base::HEX,base::EDGE>::FaceTable
            initialiseFaceTable<base::HEX,base::EDGE>()
            {
                const ElementFaces<base::HEX,base::EDGE>::FaceTable aux = {{
                        {{ 0, 1 }}, {{ 1, 2 }}, {{ 2, 3 }}, {{ 3, 0 }},
                        {{ 4, 5 }}, {{ 5, 6 }}, {{ 6, 7 }}, {{ 7, 4 }},
                        {{ 0, 4 }}, {{ 1, 5 }}, {{ 2, 6 }}, {{ 3, 7 }}
                    }};
                return aux;
            }

            //==================================================================
            // FACES
            //==================================================================

            //! Dummies for TRI and QUAD
            template<> inline
            ElementFaces<base::TRI,base::FACE>::FaceTable
            initialiseFaceTable<base::TRI,base::FACE>()
            {
                ElementFaces<base::TRI,base::FACE>::FaceTable aux = {{
                        {{ 0, 1, 2 }} }};
                return aux;
            }

            template<> inline
            ElementFaces<base::QUAD,base::FACE>::FaceTable
            initialiseFaceTable<base::QUAD,base::FACE>()
            {
                ElementFaces<base::QUAD,base::FACE>::FaceTable aux = {{
                        {{ 0, 1, 2, 3 }} }};
                return aux;
            }

            //------------------------------------------------------------------
            /** base::TET (view on faces from outside, i.e. against positive normal)
             *  <pre>
             *
             *                3              [0]:  1               [1]:  3     
             *              ,/|`\                  | \                   | \   
             *            ,/  |  `\                |  \                  |  \  
             *          ,/    '.   `\              |   \                 |   \ 
             *        ,/       |     `\            0----2                0----1
             *      ,/         |       `\  
             *     0-----------'.--------2 
             *      `\.         |      ,/    [2]:  3               [3]:  3     
             *         `\.      |    ,/            | \                   | \   
             *            `\.   '. ,/              |  \                  |  \  
             *               `\. |/                |   \                 |   \ 
             *                  `1                 1----2                2----0
             *  </pre>
             */
            template<> inline
            ElementFaces<base::TET,base::FACE>::FaceTable
            initialiseFaceTable<base::TET,base::FACE>()
            {
                ElementFaces<base::TET,base::FACE>::FaceTable aux = {{
                        {{ 0, 2, 1 }},
                        {{ 0, 1, 3 }},
                        {{ 1, 2, 3 }},
                        {{ 2, 0, 3 }}
                    }};
                return aux;
            }

            //----------------------------------------------------------------------
            /** base::HEX (view on faces from outside, i.e. against positive normal)
             *  <pre>
             *
             *                               
             *     4----------7        [0]:  1----2    [1]:  7----6    [2]:  4----5
             *     |\         |\             |    |          |    |          |    |
             *     | \        | \            |    |          |    |          |    |
             *     |  \       |  \           0----3          4----5          0----1
             *     |   5----------6                 
             *     |   |      |   |        
             *     0---|------3   |    [3]:  5----6    [4]:  6----7    [5]:  7----4
             *      \  |       \  |          |    |          |    |          |    |             
             *       \ |        \ |          |    |          |    |          |    |             
             *        \|         \|          1----2          2----3          3----0             
             *         1----------2
             *  </pre>
             */
            template<> inline
            ElementFaces<base::HEX,base::FACE>::FaceTable
            initialiseFaceTable<base::HEX,base::FACE>()
            {
                ElementFaces<base::HEX,base::FACE>::FaceTable aux = {{
                        {{ 0, 3, 2, 1 }},
                        {{ 4, 5, 6, 7 }},
                        {{ 0, 1, 5, 4 }},
                        {{ 1, 2, 6, 5 }},
                        {{ 2, 3, 7, 6 }},
                        {{ 3, 0, 4, 7 }}
                    }};
                return aux;
            }
                
        } // detail_
        
        // Call specialised initialiser function
        template<base::Shape SHAPE, base::NFace NFACE>
        const typename ElementFaces<SHAPE,NFACE>::FaceTable
        ElementFaces<SHAPE,NFACE>::faceTable_ =
            detail_::initialiseFaceTable<SHAPE,NFACE>();


    }
}

//==============================================================================

//------------------------------------------------------------------------------
/** Provides a connectivity look-up table for edges of the faces of an element.
 *
 *  \param SHAPE  Type of geometric shape of an element
 */
template<base::Shape SHAPE>
struct base::mesh::FaceEdges
{
    //! Number of n-faces of element
    static const unsigned numFaces = base::NumNFaces<SHAPE,base::FACE>::value;

    //! Dimension of element shape
    static const unsigned shapeDim = base::ShapeDim<SHAPE>::value;

    //! Co-dimension of the face w.r.t to shape
    static const unsigned coDimension = shapeDim - 2;

    //! Shape of a face
    static const base::Shape faceShape =
        base::FaceShape<SHAPE,coDimension>::value;

    //! Number of edges on a face
    static const unsigned numFaceEdges =
        base::NumNFaces<faceShape,base::EDGE>::value;

    //! Storage of edge numbers per face
    typedef boost::array<int,numFaceEdges>  EdgesPerFace;

    //! Storage of all surface arrays
    typedef boost::array<EdgesPerFace, numFaces> EdgeTable;

    //! Return shape's edge number corresponding to face,edgeIndex
    static int index( const unsigned edgeIndex,
                      const unsigned faceNum  )
    {
        return edgeTable_[faceNum][edgeIndex];
    }

    //! Return the orientation of ....
    static int sign( const unsigned edgeIndex,
                     const unsigned faceNum )
    {
        return signTable_[faceNum][edgeIndex];
    }

private:
    /** Lookup table:
     *   edgeTable_[ i ][ j ]
     *  gives the j-th edge number of the i-th face
     */
    static const EdgeTable edgeTable_;
    static const EdgeTable signTable_;
            
    //------------------------------------------------------------------
    //! Sanity check
    STATIC_ASSERT_MSG( (shapeDim > 1), "Cannot collect edges from no face" );
};

namespace base{
    namespace mesh{

        namespace detail_{
        
            //==================================================================
            // The numbers below result from combining the edge numbers of
            // ElementFaces<SHAPE,base::EDGE> with the face numbers of
            // ElementFaces<SHAPE,base::FACE> in the pictures and implementations
            // above. Hence, if you touch any of the above numbers you MUST also
            // adapt the following numbers!!
            
            //! Generic definition, never implemented
            template<base::Shape SHAPE>
            typename FaceEdges<SHAPE>::EdgeTable
            initialiseEdgeTable();

            template<base::Shape SHAPE>
            typename FaceEdges<SHAPE>::EdgeTable
            initialiseSignTable();

            // TRI
            template<> inline
            FaceEdges<base::TRI>::EdgeTable initialiseEdgeTable<base::TRI>()
            {
                FaceEdges<base::TRI>::EdgeTable aux = {{ {{ 0, 1, 2 }} }};
                return aux;
            }

            template<> inline
            FaceEdges<base::TRI>::EdgeTable initialiseSignTable<base::TRI>()
            {
                FaceEdges<base::TRI>::EdgeTable aux = {{ {{ +1, +1, +1 }} }};
                return aux;
            }

            // QUAD
            template<> inline
            FaceEdges<base::QUAD>::EdgeTable initialiseEdgeTable<base::QUAD>()
            {
                FaceEdges<base::QUAD>::EdgeTable aux = {{ {{ 0, 1, 2, 3 }} }};
                return aux;
            }

            template<> inline
            FaceEdges<base::QUAD>::EdgeTable initialiseSignTable<base::QUAD>()
            {
                FaceEdges<base::QUAD>::EdgeTable aux = {{ {{ +1, +1, +1, +1 }} }};
                return aux;
            }

            // TET
            template<> inline
            FaceEdges<base::TET>::EdgeTable initialiseEdgeTable<base::TET>()
            {
                FaceEdges<base::TET>::EdgeTable aux = {{
                        {{ 2, 1, 0 }},
                        {{ 0, 4, 3 }},
                        {{ 1, 5, 4 }},
                        {{ 2, 3, 5 }}
                    }};
                return aux;
            }

            template<> inline
            FaceEdges<base::TET>::EdgeTable initialiseSignTable<base::TET>()
            {
                FaceEdges<base::TET>::EdgeTable aux = {{
                        {{ -1, -1, -1 }},
                        {{ +1, +1, -1 }},
                        {{ +1, +1, -1 }},
                        {{ +1, +1, -1 }}
                    }};
                return aux;
            }

            // HEX
            template<> inline
            FaceEdges<base::HEX>::EdgeTable initialiseEdgeTable<base::HEX>()
            {
                FaceEdges<base::HEX>::EdgeTable aux = {{
                        {{  3,  2,  1,  0 }},
                        {{  4,  5,  6,  7 }},
                        {{  0,  9,  4,  8 }},
                        {{  1, 10,  5,  9 }},
                        {{  2, 11,  6, 10 }},
                        {{  3,  8,  7, 11 }}
                    }};
                return aux;
            }

            template<> inline
            FaceEdges<base::HEX>::EdgeTable initialiseSignTable<base::HEX>()
            {
                FaceEdges<base::HEX>::EdgeTable aux = {{
                        {{ -1, -1, -1, -1 }},
                        {{ +1, +1, +1, +1 }},
                        {{ +1, +1, -1, -1 }},
                        {{ +1, +1, -1, -1 }},
                        {{ +1, +1, -1, -1 }},
                        {{ +1, +1, -1, -1 }}
                    }};
                return aux;
            }
        } // namespace detail_

        // Call specialised initialiser function
        template<base::Shape SHAPE>
        const typename FaceEdges<SHAPE>::EdgeTable
        FaceEdges<SHAPE>::edgeTable_ =
            detail_::initialiseEdgeTable<SHAPE>();

        template<base::Shape SHAPE>
        const typename FaceEdges<SHAPE>::EdgeTable
        FaceEdges<SHAPE>::signTable_ =
            detail_::initialiseSignTable<SHAPE>();

    }
}

#endif
