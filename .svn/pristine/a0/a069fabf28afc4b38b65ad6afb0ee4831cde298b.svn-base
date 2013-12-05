//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ParameterSurface.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_parametersurface_hpp
#define base_mesh_parametersurface_hpp

//------------------------------------------------------------------------------
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<base::Shape SHAPE>
        struct ParameterSurface
        {
            //! Shape of volume element
            static const base::Shape shape = SHAPE;

            //! Dimension of the shape
            static const unsigned dim = base::ShapeDim<shape>::value;

            //! Number of faces of volume element
            static const unsigned numSurfaces = base::NumMFaces<shape,dim-1>::value;

            //! Shape of surface of volume element
            static const base::Shape faceShape =
                base::FaceShape<shape>::value;

            //! Number of vertices of surface element
            static const unsigned numFaceVertices =
                base::NumMFaces<faceShape,base::VERTEX>::value;

            //! Storage of surface element connectivity
            typedef boost::array<unsigned, numFaceVertices> Surface;

            //! Storage of all surface arrays
            typedef boost::array<Surface, numSurfaces> SurfaceTable;

            /** Lookup table:
             *    surfaceTable[ i ][ j ]
             *  gives the j-th vertex index of the i-th surface
             */
            static const SurfaceTable surfaceTable;
        };

        //----------------------------------------------------------------------
        /**  base::LINE
         *   <pre>
         *
         *                 [0] 0-----1 [1]
         *   </pre>
         */  
        template<>
        const typename ParameterSurface<base::LINE>::SurfaceTable
        ParameterSurface<base::LINE>::surfaceTable = {{
                {{ 0 }},
                {{ 1 }}
            }};

        //----------------------------------------------------------------------
        /**  base::TRIANGLE
         *   <pre>
         *
         *   2     2         2
         *   V     |`\        |\
         *   |     |  `\        `\
         *   |[0]  |    `\    [2] `\
         *   |     |      `\        `\
         *   V     |        `\        |\
         *   0     0----------1         1
         *     
         *             [1]
         *         0>-------->1
         *   </pre>
         */
        template<>
        const typename ParameterSurface<base::TRI>::SurfaceTable
        ParameterSurface<base::TRI>::surfaceTable = {{
                {{ 2, 0 }},
                {{ 0, 1 }},
                {{ 1, 2 }}
            }};

        //----------------------------------------------------------------------
        /** base::QUAD
         *  <pre>
         *              2<---------<3
         *                   [3] (xi_1=1)
         *
         *      2       2-----------3       3
         *      V       |           |       A
         *      |       |           |       |
         *      |[0]    |(xi_0=0)   |    [2]| (xi_0=1)
         *      |       |           |       |
         *      V       |           |       A
         *      0       0-----------1       1
         *  
         *                   [1] (xi_1=0)
         *              0>--------->1
         *  </pre>
         */
        template<>
        const typename ParameterSurface<base::QUAD>::SurfaceTable
        ParameterSurface<base::QUAD>::surfaceTable = {{
                {{ 2, 0 }},
                {{ 0, 1 }},
                {{ 1, 3 }},
                {{ 3, 2 }}
            }};

        //----------------------------------------------------------------------
        /** base::TET (view on faces from outside, i.e. against positive normal)
         *  <pre>
         *
         *                3              [0]:  2               [1]:  3     
         *              ,/|`\                  | \                   | \   
         *            ,/  |  `\                |  \                  |  \  
         *          ,/    '.   `\              |   \                 |   \ 
         *        ,/       |     `\            0----3                0----1
         *      ,/         |       `\  
         *     0-----------'.--------2 
         *      `\.         |      ,/    [2]:  1               [3]:  3     
         *         `\.      |    ,/            | \                   | \   
         *            `\.   '. ,/              |  \                  |  \  
         *               `\. |/                |   \                 |   \ 
         *                  `1                 0----2                1----2
         *  </pre>
         */
        template<>
        const typename ParameterSurface<base::TET>::SurfaceTable
        ParameterSurface<base::TET>::surfaceTable = {{
                {{ 0, 3, 2 }},
                {{ 0, 1, 3 }},
                {{ 0, 2, 1 }},
                {{ 1, 2, 3 }}
            }};

        //----------------------------------------------------------------------
        /** base::HEX (view on faces from outside, i.e. against positive normal)
         *  <pre>
         *
         *                               xi_0=0          xi_1=0          xi_=0
         *     4----------6        [0]:  6----4    [1]:  4----5    [2]:  1----3
         *     |\         |\             |    |          |    |          |    |
         *     | \        | \            |    |          |    |          |    |
         *     |  \       |  \           2----0          0----1          0----2
         *     |   5----------7                 
         *     |   |      |   |          xi_0=1          xi_1=1          xi_2=1
         *     0---|------2   |    [3]:  5----7    [4]:  7----6    [5]:  6----7
         *      \  |       \  |          |    |          |    |          |    |             
         *       \ |        \ |          |    |          |    |          |    |             
         *        \|         \|          1----3          3----2          4----5             
         *         1----------3
         *  </pre>
         */
        template<>
        const typename ParameterSurface<base::HEX>::SurfaceTable
        ParameterSurface<base::HEX>::surfaceTable = {{
                {{ 2, 0, 6, 4 }},
                {{ 0, 1, 4, 5 }},
                {{ 0, 2, 1, 3 }},
                {{ 1, 3, 5, 7 }},
                {{ 3, 2, 7, 6 }},
                {{ 4, 5, 6, 7 }}
            }};


    }
}


#endif
