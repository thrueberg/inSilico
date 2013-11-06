//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   cylinder.cpp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <vector>
// boost includes
#include <boost/lexical_cast.hpp>
// tools include
#include <tools/meshGeneration/meshGeneration.hpp>
#include <tools/meshGeneration/boundaries3D/ParametricSurface.hpp>

//------------------------------------------------------------------------------
/** Parametric representation of a cylinder surface.  This function
 *  reads all the relevant parameters and generates the map
 *  \f[
 *  [0,1) x [0,1) \to R^3: x = C(\xi_1,\xi_2)
 *  \f] where \f$ \xi \f$ is
 *  two-dimensional parameteric coordinate and \f$ C \f$ a
 *  three-dimensional point on the surface of a cylinder.  The
 *  cylinder is perpendicular to the 1-2 plane, is defined by the
 *  centre point, its radius and height. For the parameterisation, the
 *  intervals in 1- and 2-direction are provided and, if wanted, the
 *  parameter space is triangulated.
 */
class Cylinder
{
public:
    typedef tools::meshGeneration::Point Point;

    //--------------------------------------------------------------------------
    //! Read all needed values from the input arguments
    Cylinder( const int argc, char* const argv[] )
        : radius_( boost::lexical_cast<double>( argv[1] ) ), 
          height_( boost::lexical_cast<double>( argv[2] ) ), 
          centre_(  tools::meshGeneration::makePoint(
                        boost::lexical_cast<double>( argv[3] ),
                        boost::lexical_cast<double>( argv[4] ),
                        boost::lexical_cast<double>( argv[5] ) ) ), 
          nIntv1_( boost::lexical_cast<unsigned>( argv[6] ) ), 
          nIntv2_( boost::lexical_cast<unsigned>( argv[7] ) ), 
        withTriangles_( ( argc == 8 ? false :
                          boost::lexical_cast<bool>( argv[8] ) ) )
    {

    }

    //--------------------------------------------------------------------------
    //! Generate point coordinates given the indices of a grid point 
    Point makePoint( const unsigned i, const unsigned j ) const
    {
        Point p = centre_;
        const double dPhi   = 2. * M_PI / static_cast<double>( nIntv1_ );
        const double dZ     = height_   / static_cast<double>( nIntv2_ );
        const double phi    = static_cast<double>( i ) * dPhi;
        const double z      = -height_/2. + static_cast<double>( j ) * dZ;
        p[0] += radius_ * std::cos( phi );
        p[1] += radius_ * std::sin( phi );
        p[2] += z;         
        return p;
    }

    //! Number of intervals in 1-direction
    unsigned numIntervals1() const { return nIntv1_; }
    //! Number of intervals in 2-direction
    unsigned numIntervals2() const { return nIntv2_; }

    //! Periodicity in 1-direction
    bool periodic1() const { return true; }
    //! Periodicity in 2-direction
    bool periodic2() const { return false; }

    //! If space is trinagulated
    bool withTriangles() const { return withTriangles_; }

    //! Write parameters to a stream
    std::ostream& write( std::ostream& out ) const
    {
        out << "Centre=("
            << centre_[0] << ", " << centre_[1] << ", "<< centre_[2]
            << "), r=" << radius_ << ", h=" << height_
            << ", n1="  << nIntv1_   << ", n2=" << nIntv2_ 
            << ", using "
            << (withTriangles_ ? "triangle" : "quadrilateral")
            << " elements ";
        
        return out;
    }
    
private:
    const double   radius_;        //!< Radius
    const double   height_;        //!< Height
    const Point    centre_;        //!< Centre of cylinder
    const unsigned nIntv1_;        //!< Number of intervals in 1-direction
    const unsigned nIntv2_;        //!< Number of intervals in 2-direction
    const bool     withTriangles_; //!< True for a triangulation of (0,1)^2
};
    
//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    if ( ( argc != 8 ) and ( argc != 9 ) ){
        std::cerr << "Usage: " << argv[0]
                  << " r h cx cy cz n1 n2 [mt=0]\n\n"
                  << "To create a cylinder with radius r and height h  around "
                  << "centre (cx,cy,cz) with a grid of n1 x n2 elements.\n"
                  << "The optional flag mt (default=false=0) forces the "
                  << "generation of triangles instead of quadrilaterals\n\n";
        return false;
    }

    Cylinder cylinder( argc, argv );
    tools::meshGeneration::ParametricSurface<Cylinder>::apply( cylinder, std::cout );
    
    return 0;
}
