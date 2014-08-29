//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   torus.cpp
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
namespace tools{
    namespace meshGeneration{
        namespace boundaries3D{
            class Torus;
        }
    }
}

//------------------------------------------------------------------------------
/** Parametric representation of a torus.
 *  This function reads all the relevant parameters and generates the map
 *  \f[
 *      [0,1) x [0,1) \to  R^3:  x = T(\xi_1,\xi_2)
 *  \f]
 *  where \f$ \xi \f$ is two-dimensional parameteric coordinate and \f$ T \f$
 *  a three-dimensional point on the surface of a torus.
 *  The torus lies in the 1-2 plane, is defined by the centre point and its
 *  two radii. For the parameterisation, the intervals in 1- and 2-direction
 *  are provided and, if wanted, the parameter space is triangulated. 
 */
class tools::meshGeneration::boundaries3D::Torus
{
public:
    typedef tools::meshGeneration::Point Point;

    //--------------------------------------------------------------------------
    //! Read all needed values from the input arguments
    Torus( const int argc, char* const argv[] )
        : radius1_( boost::lexical_cast<double>( argv[1] ) ), 
          radius2_( boost::lexical_cast<double>( argv[2] ) ), 
          centre_(  tools::meshGeneration::makePoint(
                        boost::lexical_cast<double>( argv[3] ),
                        boost::lexical_cast<double>( argv[4] ),
                        boost::lexical_cast<double>( argv[5] ) ) ), 
          nIntv1_( boost::lexical_cast<unsigned>( argv[6] ) ), 
          nIntv2_( boost::lexical_cast<unsigned>( argv[7] ) ), 
        withTriangles_( ( argc == 8 ? false :
                          boost::lexical_cast<bool>( argv[8] ) ) )
    {
        // a word to the user
        if ( radius1_ <= radius2_)
        {
            std::cerr << "(WW) The major radius is NOT greater than the minor radius.\n"
                      << "(WW) The torus will be degenerate. \n\n";
        }

    }

    //--------------------------------------------------------------------------
    //! Generate point coordinates given the indices of a grid point 
    Point makePoint( const unsigned i, const unsigned j ) const
    {
        Point p = centre_;
        const double dPhi   = 2. * M_PI / static_cast<double>( nIntv1_ );
        const double dTheta = 2. * M_PI / static_cast<double>( nIntv2_ );
        const double phi    = static_cast<double>( i ) * dPhi;
        const double theta  = static_cast<double>( j ) * dTheta;
        p[0] += ( radius1_ + radius2_ * std::cos( theta ) )* std::cos( phi );
        p[1] += ( radius1_ + radius2_ * std::cos( theta ) )* std::sin( phi );
        p[2] +=              radius2_ * std::sin( theta );
        return p;
    }

    //! Number of intervals in 1-direction
    unsigned numIntervals1() const { return nIntv1_; }
    //! Number of intervals in 2-direction
    unsigned numIntervals2() const { return nIntv2_; }

    //! Periodicity in 1-direction
    bool periodic1() const { return true; }
    //! Periodicity in 2-direction
    bool periodic2() const { return true; }

    //! If space is trinagulated
    bool withTriangles() const { return withTriangles_; }

    //! Write parameters to a stream
    std::ostream& write( std::ostream& out ) const
    {
        out << "Centre=("
            << centre_[0] << ", " << centre_[1] << ", "<< centre_[2]
            << "), r1=" << radius1_ << ", r2=" << radius2_
            << ", n1="  << nIntv1_   << ", n2=" << nIntv2_ 
            << ", using "
            << (withTriangles_ ? "triangle" : "quadrilateral")
            << " elements ";
        
        return out;
    }
    
private:
    const double   radius1_;       //!< Major radius
    const double   radius2_;       //!< Minor radius
    const Point    centre_;        //!< Centre of torus
    const unsigned nIntv1_;        //!< Number of intervals in 1-direction
    const unsigned nIntv2_;        //!< Number of intervals in 2-direction
    const bool     withTriangles_; //!< True for a triangulation of (0,1)^2
};
    
//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    if ( ( argc != 8 ) and ( argc != 9 ) ){
        std::cerr << "Usage: " << argv[0]
                  << " r1 r2 cx cy cz n1 n2 [mt=0]\n\n"
                  << "To create a torus with radii r1 and r2 (r1 > r2) around "
                  << "centre (cx,cy,cz) with a grid of n1 x n2 elements.\n"
                  << "The optional flag mt (default=false=0) forces the "
                  << "generation of triangles instead of quadrilaterals\n\n";
        return false;
    }

    typedef tools::meshGeneration::boundaries3D::Torus Torus;
    Torus torus( argc, argv );
    tools::meshGeneration::boundaries3D::ParametricSurface<Torus>::apply( torus,
                                                                          std::cout );
    
    return 0;
}
