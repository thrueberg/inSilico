//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ParametricSurface.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_meshgeneration_parametricsurface_h
#define tools_meshgeneration_parametricsurface_h

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <sstream>
// tools include
#include <tools/meshGeneration/meshGeneration.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace meshGeneration{
        namespace boundaries3D{
            
            template<typename SURF>
            struct ParametricSurface;
            
        }
    }
}
//------------------------------------------------------------------------------
/** Generate mesh of a parametric surface given the coordinate map.
 *  This function generates a mesh from the map
 *  \f[
 *      [0,1) x [0,1) \to  R^3:  x = S(\xi_1,\xi_2)
 *  \f]
 *  where \f$ \xi \f$ is two-dimensional parameteric coordinate and \f$ S \f$
 *  a three-dimensional point on the surface.
 *  \tparam SURF Type of surface which provides the map
 */
template<typename SURF>
struct tools::meshGeneration::boundaries3D::ParametricSurface
{
    typedef tools::meshGeneration::Point Point;

    //--------------------------------------------------------------------------
    /** Generate points
     */
    static void makePoints( const SURF& surf,
                            std::vector<Point>& points )
    {
        const unsigned n1 =
            surf.numIntervals1() + (surf.periodic1() ? 0 : 1 );
        const unsigned n2 =
            surf.numIntervals2() + (surf.periodic2() ? 0 : 1 );
        
        for ( unsigned j = 0; j < n2; j++ ) {
            for ( unsigned i = 0; i < n1; i++ ) {
                points.push_back( surf.makePoint( i, j ) );
            }
        }

    }

    typedef tools::meshGeneration::Element Element;

    //--------------------------------------------------------------------------
    /** Generate elements
     */
    static void makeElements( const SURF& surf,
                              std::vector<Element>& elements )
    {
        const unsigned n1 = surf.numIntervals1();
        const unsigned n2 = surf.numIntervals2();
        
        for ( unsigned j = 0; j < n2; j++ ) {
            for ( unsigned i = 0; i < n1; i++ ) {

                // following indices with boundary periodicity
                const std::size_t iNext =
                    ( i < n1 - 1 ? i+1 :
                      (surf.periodic1() ? 0 : i+1) );
                const std::size_t jNext = 
                    ( j < n2 - 1 ? j+1 :
                      (surf.periodic2() ? 0 : j+1) );

                // index shift to next line
                const unsigned stride =
                    n1 + static_cast<unsigned>(not surf.periodic1());
                
                // the four vertices of a quadrilateral
                const std::size_t i1 = j     * stride + i;
                const std::size_t i2 = j     * stride + iNext;
                const std::size_t i4 = jNext * stride + i;
                const std::size_t i3 = jNext * stride + iNext;

                if ( surf.withTriangles() ) {
                    // generate two triangles
                    Element tri1, tri2;
                    tri1.push_back( i1 );
                    tri1.push_back( i2 );
                    tri1.push_back( i3 );
                    tri2.push_back( i3 );
                    tri2.push_back( i4 );
                    tri2.push_back( i1 );
                    elements.push_back( tri1 );
                    elements.push_back( tri2 );
                
                }
                else {
                    // generate a quadrilateral
                    Element quad;
                    quad.push_back( i1 );
                    quad.push_back( i2 );
                    quad.push_back( i3 );
                    quad.push_back( i4 );
                    elements.push_back( quad );
                }
            }
        }
        return;
    }

    //--------------------------------------------------------------------------
    /** Write to stream
     */
    static std::ostream& write( const SURF& surf, 
                                const std::vector<Point>& points,
                                const std::vector<Element>& elements, 
                                std::ostream& out )
    {
        // write output to stream
        const bool makeTriangles = surf.withTriangles();
        const std::string elementShape = ( makeTriangles ? "triangle" : "quadrilateral" );
        const unsigned    elementNumPoints = ( makeTriangles ? 3 : 4 );

        std::stringstream buffer;
        surf.write( buffer );
        
        tools::meshGeneration::writeSMFComment( buffer.str(), out );
        tools::meshGeneration::writeSMFHeader( elementShape, elementNumPoints, 
                                               points.size(), elements.size(), out );
        tools::meshGeneration::writePoints(   points,   out );
        tools::meshGeneration::writeElements( elements, out );

        return out;
    }

    //--------------------------------------------------------------------------
    /**
     */
    static void apply( const SURF& surf, std::ostream& out )
    {
        std::vector<Point> points;
        makePoints( surf, points );

        std::vector<Element> elements;
        makeElements( surf, elements );

        write( surf, points, elements, out );
    }

};


#endif
