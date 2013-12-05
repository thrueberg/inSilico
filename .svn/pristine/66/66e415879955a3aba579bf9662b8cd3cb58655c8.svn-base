//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   boundaries2D.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_meshgeneration_boundaries2D_boundaries2D_h
#define tools_meshgeneration_boundaries2D_boundaries2D_h


#include <tools/meshGeneration/meshGeneration.hpp>


//------------------------------------------------------------------------------
namespace tools{
    namespace meshGeneration{
        namespace boundaries2D{

            //------------------------------------------------------------------
            // generate chain  of elements
            void generateElements( const std::size_t nNodes, 
                                   std::vector<tools::meshGeneration::Element> & elements, 
                                   const bool closed = true ) 
            {
                const std::size_t nElements = nNodes - ( closed ? 0 : 1 );
                for ( std::size_t e = 0; e < nElements; e ++ ) {
                    const std::size_t v1 = e;
                    const std::size_t v2 = (closed and (e == nElements - 1) ) ? 0 : v1 + 1 ;
                    tools::meshGeneration::Element el;
                    el.push_back( v1 );
                    el.push_back( v2 );
                    elements.push_back( el );
                }
            }

            
        }
    }
}

#if 0
//------------------------------------------------------------------------------
// Functor for a coordinate rotation around (centerX,centerY) with angle theta
class PointRotator 
    : public std::unary_function<Point &,void>
{
public:
    PointRotator( double centerX, double centerY, double angle ) 
        : centerX_(centerX), 
          centerY_(centerY), 
          cosTheta_( cos( angle ) ), 
          sinTheta_( sin( angle ) )
    { }
                    
    void operator()( Point & p )
    {
        const double tx = p.x; 
        const double ty = p.y;
        p.x = cosTheta_ * (tx - centerX_) - sinTheta_ * (ty - centerY_) + centerX_;
        p.y = sinTheta_ * (tx - centerX_) + cosTheta_ * (ty - centerY_) + centerY_;
        return;
    }
                    
private:
    double centerX_, centerY_;
    double cosTheta_, sinTheta_;
};

// rotate points
PointRotator pointRotation( std::vector<Point> & points )
{
    std::cout << "Do you want to rotate the object? (0/1)";
    bool yes;
    std::cin >> yes;
    if ( yes ) {
        double centerX, centerY, angle;
        std::cout << " Enter center of rotation (x/y)" << std::endl
                  << " x = ";
        std::cin >> centerX;
        std::cout << " y = ";                   
        std::cin >> centerY;
        std::cout << " Enter angle of rotation in degrees : ";
        std::cin  >> angle;
        PointRotator pointRotator( centerX, centerY, M_PI * angle / 180. );
        pointRotator = std::for_each( points.begin(), points.end(), pointRotator );
        return pointRotator;
    }
    PointRotator dummy( 0., 0., 0. );
    return dummy;
}

// write the mesh to stream
std::ostream & writeMesh( std::ostream & out, 
                          const std::vector<Point>   & points, 
                          const std::vector<Element> & elements ) 
{
    // write header
    corlib::SmfHead smfHead;
    smfHead.setElementShape( corlib::LINE );
    smfHead.setElementNumPoints( 2 );
    smfHead.write( out );
    // write smf data
    out << points.size() << "  " << elements.size() << std::endl;
    std::for_each(    points.begin(),  points.end(), boost::bind2nd( PointWriter(),   out ) );
    std::for_each( elements.begin(), elements.end(), boost::bind2nd( ElementWriter(), out ) );

    return out;
}
#endif


//------------------------------------------------------------------------------
#endif
