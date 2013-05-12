//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   elementTypes.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_converter_gmsh2smf_elementtypes_hpp
#define tools_converter_gmsh2smf_elementtypes_hpp

// std includes
#include <string>
#include <algorithm>
// boost includes
#include <boost/array.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace gmsh2smf{
            //! Type of node: just an array of 3 numbers
            typedef boost::array<double,3> Node;

            //! get number of nodes depending on element type
            unsigned numNodesPerElement( const unsigned type );

            //! get smf element type name  depending on element type
            std::string smfNameOfElementType( const unsigned type );

            //! Reorder GMSH index logic to SMF index logic
            void reorderConnectivity( const unsigned type,
                                      std::vector<std::size_t>& connec );

        }
    }
}
    
//------------------------------------------------------------------------------
/** Definition of the elment types in GMSH
 *  Copied from: http://www.geuz.org/gmsh/doc/texinfo/gmsh.html
 *
 *  GMSH element ID | Description
 *  --------------: | :----------------------------------------------
 *          1       |    2-node line                               
 *          2       |    3-node triangle                          
 *          3       |    4-node quadrangle                        
 *          4       |    4-node tetrahedron                       
 *          5       |    8-node hexahedron                        
 *          8       |    3-node second order line                 
 *          9       |    6-node second order triangle             
 *         10       |    9-node second order quadrangle           
 *         11       |    10-node second order tetrahedron         
 *         12       |    27-node second order hexahedron          
 *         15       |    1-node point                             
 *         16       |    8-node second order quadrangle           
 *         17       |    20-node second order hexahedron          
 *         21       |    10-node third order triangle             
 *         23       |    15-node fourth order triangle            
 *         25       |    21-node fifth order complete triangle    
 *         26       |    4-node third order edge                  
 *         27       |    5-node fourth order edge                 
 *         28       |    6-node fifth order edge                  
 *         29       |    20-node third order tetrahedron          
 *         30       |    35-node fourth order tetrahedron         
 *         31       |    56-node fifth order tetrahedron          
 *  
 *
 * \code{.txt}
 *  
 *  Line:                   Line3:           Line4:
 *       
 *  0----------1 --> u      0-----2----1     0----2----3----1
 *  
 *  
 *  Triangle:               Triangle6:          Triangle10:          Triangle15:
 *       
 *  v
 *  ^                                                                   2
 *  |                                                                   | \
 *  2                       2                    2                      9   8
 *  |`\                     |`\                  | \                    |     \
 *  |  `\                   |  `\                7   6                 10  14   7
 *  |    `\                 5    `4              |     \                |         \
 *  |      `\               |      `\            8   9   5             11  12   13  6
 *  |        `\             |        `\          |         \            |             \
 *  0----------1 --> u      0-----3----1         0---3---4---1          0---3---4---5---1
 *  
 *  
 *  
 *  Quadrangle:            Quadrangle8:            Quadrangle9:
 *       
 *  v
 *  ^
 *  |
 *  3-----------2          3-----6-----2           3-----6-----2
 *  |     |     |          |           |           |           |
 *  |     |     |          |           |           |           |
 *  |     +---- | --> u    7           5           7     8     5
 *  |           |          |           |           |           |
 *  |           |          |           |           |           |
 *  0-----------1          0-----4-----1           0-----4-----1
 *  
 *  
 *
 *     Tetrahedron:                          Tetrahedron10:
 *     
 *                        v
 *                      .
 *                    ,/
 *                   /
 *                2                                     2
 *              ,/|`\                                 ,/|`\
 *            ,/  |  `\                             ,/  |  `\
 *          ,/    '.   `\                         ,6    '.   `5
 *        ,/       |     `\                     ,/       8     `\
 *      ,/         |       `\                 ,/         |       `\
 *     0-----------'.--------1 --> u         0--------4--'.--------1
 *      `\.         |      ,/                 `\.         |      ,/
 *         `\.      |    ,/                      `\.      |    ,9
 *            `\.   '. ,/                           `7.   '. ,/
 *               `\. |/                                `\. |/
 *                  `3                                    `3
 *                     `\.
 *                        ` w
 *
 *
 *     Hexahedron:             Hexahedron20:          Hexahedron27:
 *     
 *            v
 *     3----------2            3----13----2           3----13----2
 *     |\     ^   |\           |\         |\          |\         |\
 *     | \    |   | \          | 15       | 14        |15    24  | 14
 *     |  \   |   |  \         9  \       11 \        9  \ 20    11 \
 *     |   7------+---6        |   7----19+---6       |   7----19+---6
 *     |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
 *     0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
 *      \  |    \  \  |         \  17      \  18       \ 17    25 \  18
 *       \ |     \  \ |         10 |        12|        10 |  21    12|
 *        \|      w  \|           \|         \|          \|         \|
 *         4----------5            4----16----5           4----16----5
 *
 *  \endcode
 *
 *  \param[in] type   Type of element (GMSH element ID)
 *  \returns          Number of nodes of the element
 */
unsigned tools::converter::gmsh2smf::numNodesPerElement( const unsigned type )
{
    unsigned result = static_cast<unsigned>( -1 );

    switch (type)
    {
    case  1: result =  2; break; //  2-node line
    case  2: result =  3; break; //  3-node triangle
    case  3: result =  4; break; //  4-node quad
    case  4: result =  4; break; //  4-node tet
    case  5: result =  8; break; //  8-node hex

    case  8: result =  3; break; //  3-node line
    case  9: result =  6; break; //  6-node triangle
    case 10: result =  9; break; //  9-node quad
    case 11: result = 10; break; // 10-node tet
    case 12: result = 27; break; // 27-node hex

    case 15: result =  1; break; //  a point
    case 16: result =  8; break; // 8-node quad (serendipity)
    case 17: result = 20; break; // 20-node hex (serendipity)

    case 21: result = 10; break; // 10-node triangle
        
    case 23: result = 15; break; // 15-node triangle
        
    case 25: result = 15; break; // 21-node triangle
    case 26: result =  4; break; //  4-node line
    case 27: result =  5; break; //  5-node line
    case 28: result =  6; break; //  6-node line
    case 29: result = 20; break; // 20-node tet
    case 30: result = 35; break; // 35-node tet
    case 31: result = 56; break; // 56-node tet
    default: ;
    }

    // user warnings
    switch( type )
    {
    case 16: case 17:
        std::cerr << "(WW) Serendipity elements not yet supported \n"
                  << "(WW) Shape functions missing and reorientation\n";
        break;

    case 23: case 25:
        std::cerr << "(WW) Triangle of fourth and fifth order not yet implemented\n";
        break;

    case 29: case 30: case 31:
        std::cerr << "(WW) Tets of orders 3, 4 and 5 currently not implemented\n";
        break;
    default: ;
    }

    return result;
}

//------------------------------------------------------------------------------
/** For details see numNodesPerElement()
 *  \param[in] type  GMSH ID of element type
 *  \returns         Human-readable aname of element (a la smf specifications)
 */
std::string tools::converter::gmsh2smf::smfNameOfElementType( const unsigned type )
{
    std::string result = "invalid";

    switch (type)
    {
    case 15:                                     result = "point";         break;
    case  1: case  8: case 26: case 27: case 28: result = "line";          break;
    case  2: case  9: case 21: case 23: case 25: result = "triangle";      break;
    case  3: case 10: case 16:                   result = "quadrilateral"; break;
    case  4: case 11: case 29: case 30: case 31: result = "tetrahedron";   break;
    case  5: case 12: case 17:                   result = "hexahedron";    break;
    default: ;
    }
    return result;
}

//------------------------------------------------------------------------------
/** Reorder the connectivity of certain element types.
 *  \param[in]     type   GMSH ID of the element type.
 *  \param[in,out] connec  Connectivity before and after reordering
 */
void tools::converter::gmsh2smf::reorderConnectivity( const unsigned type,
                                                      std::vector<std::size_t>& connec )
{
    switch( type )
    {
    case 11: // Tet-10 needs reordering
        std::swap( connec[8], connec[9] );
        break;

    case 27:
    {
        std::vector<std::size_t> tmp = connec;
        const boost::array<unsigned,18> eightPlusIdx =
            {{ 8, 11, 13,  9, 16, 18, 19, 17, 10, 12, 14, 15, // edges
               24, 25, 21, 23, 24, 22 }};                      // faces

        for ( unsigned i = 0; i < 17; i++ )
            tmp[ 8+i ] = connec[ eightPlusIdx[i] ];
        
        connec = tmp;
    }
    break;
    default: ; // do not reorder
    }

    return;
}

#endif
