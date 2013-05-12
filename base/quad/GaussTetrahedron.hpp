//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   GaussTetrahedron.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_quad_gausstetrahedron_hpp
#define base_quad_gausstetrahedron_hpp

//------------------------------------------------------------------------------
// std   includes
#include <utility>
// boost includes
#include <boost/array.hpp>
#include <boost/utility.hpp>
// base  includes
#include <base/verify.hpp>
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace quad{
        
        template<unsigned DEGREE>
        class GaussTetrahedron;

        // Declare a lookup table for the require number of points
        namespace detail_{

            template<unsigned DEGREE> struct NumTetPoints;
            template<> struct NumTetPoints<1>{ static const unsigned value =  1; };
            template<> struct NumTetPoints<2>{ static const unsigned value =  4; };
            template<> struct NumTetPoints<3>{ static const unsigned value =  5; };
            template<> struct NumTetPoints<4>{ static const unsigned value = 11; };
            template<> struct NumTetPoints<5>{ static const unsigned value = 15; };
           
        }
    }
}

//------------------------------------------------------------------------------
/** \brief Gauss quadrature rules for tetrahedra with given polynomial degree
 *  \tparam DEGREE Polynomial degree to be integrated exactly
 */
template<unsigned DEGREE> 
class base::quad::GaussTetrahedron
    : public boost::noncopyable
{
    //! Sanity check NUMPOINTS > 0
    STATIC_ASSERT_MSG( (DEGREE > 0),
                       "Degree has to be larger than zero" );

public:
    //! Template parameter: polynomial degree to be integrated exactly
    static const unsigned degree = DEGREE;

    //! Deduce number of quadrature points from degree and lookup table
    static const unsigned numPoints =
        base::quad::detail_::NumTetPoints<degree>::value;

    //! The local coordinate dimension
    static const unsigned dim = 3;

    //! Type of coordinate-vector
    typedef typename base::Vector<dim>::Type VecDim;

    //! Type of iterator for external access
    typedef typename boost::array< std::pair<double, VecDim>,
                                   numPoints>::const_iterator Iter;

    //! Specialized constructor for different polynomial degrees
    GaussTetrahedron();

    //! Begin of array iterator
    Iter begin() const  { return weightsAndPoints_.begin(); }
    //! End of array iterator
    Iter end()   const  { return weightsAndPoints_.end();   }
    
private:
    //! Pair of weights and  points representing the quadrature rule
    boost::array<std::pair<double,VecDim>, numPoints> weightsAndPoints_;
};

//------------------------------------------------------------------------------
namespace base{
    namespace quad{
        
        //----------------------------------------------------------------------
        //! Mid-point rule, exact for linears
        template<> GaussTetrahedron<1>::GaussTetrahedron()
        {
            const VecDim point( 0.25, 0.25, 0.25 );
            weightsAndPoints_[0] = std::make_pair( 0.166666666666666, point );
        }

        //----------------------------------------------------------------------
        //! 4-point rule, exact for quadratic polynomials
        template<> GaussTetrahedron<2>::GaussTetrahedron()
        {
            const double a = 0.13819660112501051518; // (5-sqrt(5))/20
            const double b = 0.58541019662496845446; // 1 - 3*a
            const double w = 0.04166666666666666667; // (1/4)*(1/6)
            
            const VecDim p1( a, a, a ); weightsAndPoints_[0] = std::make_pair( w, p1 );
            const VecDim p2( b, a, a ); weightsAndPoints_[1] = std::make_pair( w, p2 );
            const VecDim p3( a, b, a ); weightsAndPoints_[2] = std::make_pair( w, p3 );
            const VecDim p4( a, a, b ); weightsAndPoints_[3] = std::make_pair( w, p4 );
        }

        //----------------------------------------------------------------------
        //! 5-point rule, exact for cubic polynomials
        template<> GaussTetrahedron<3>::GaussTetrahedron()
        {
            const double x  =   0.1666666666666667; // 1/6
            // weights: -4/30, 3/40
            const double w1 = -0.13333333333333333;
            
            const VecDim p1(0.25, 0.25, 0.25); weightsAndPoints_[0] = std::make_pair( w1,    p1 );
            const VecDim p2(   x,    x,    x); weightsAndPoints_[1] = std::make_pair( 0.075, p2 );
            const VecDim p3( 0.5,    x,    x); weightsAndPoints_[2] = std::make_pair( 0.075, p3 );
            const VecDim p4(   x,  0.5,    x); weightsAndPoints_[3] = std::make_pair( 0.075, p4 );
            const VecDim p5(   x,    x,  0.5); weightsAndPoints_[4] = std::make_pair( 0.075, p5 );
        }

        //----------------------------------------------------------------------
        //! 11-point rule, exact for fourth-order polynomials
        template<> GaussTetrahedron<4>::GaussTetrahedron()
        {
            const VecDim p1( 0.25, 0.25, 0.25);
            weightsAndPoints_[0] = std::make_pair( -0.01315555555555555556, p1 );

            double x = 0.071428571428571, y = 0.785714285714286, w = 0.0076222222222222;
            const VecDim p2( x, x, x); weightsAndPoints_[1] = std::make_pair( w, p2 );
            const VecDim p3( y, x, x); weightsAndPoints_[2] = std::make_pair( w, p3 );
            const VecDim p4( x, y, x); weightsAndPoints_[3] = std::make_pair( w, p4 );
            const VecDim p5( x, x, y); weightsAndPoints_[4] = std::make_pair( w, p5 );

            x = 0.100596423833201, y = 0.399403576166799, w = 0.024888888888889;
            const VecDim p6(  y, y, x); weightsAndPoints_[ 5] = std::make_pair( w, p6  );
            const VecDim p7(  y, x, x); weightsAndPoints_[ 6] = std::make_pair( w, p7  );
            const VecDim p8(  x, y, x); weightsAndPoints_[ 7] = std::make_pair( w, p8  );
            const VecDim p9(  x, x, y); weightsAndPoints_[ 8] = std::make_pair( w, p9  );
            const VecDim p10( y, x, y); weightsAndPoints_[ 9] = std::make_pair( w, p10 );
            const VecDim p11( x, y, y); weightsAndPoints_[10] = std::make_pair( w, p11 );
        }

        //----------------------------------------------------------------------
        //! 15-point rule, exact for fifth-order polynomials
        template<> GaussTetrahedron<5>::GaussTetrahedron()
        {
            const VecDim p1( 0.25, 0.25, 0.25);
            weightsAndPoints_[0] = std::make_pair( 0.030283678097089, p1 );

            double x = 0.333333333333333, y = 0.0, w = 0.006026785714286;
            const VecDim p2( x, x, x ); weightsAndPoints_[1] = std::make_pair( w, p2 );
            const VecDim p3( y, x, x ); weightsAndPoints_[2] = std::make_pair( w, p3 );
            const VecDim p4( x, y, x ); weightsAndPoints_[3] = std::make_pair( w, p4 );
            const VecDim p5( x, x, y ); weightsAndPoints_[4] = std::make_pair( w, p5 );

            x = 0.090909090909091, y = 0.727272727272727, w = 0.011645249086029;
            const VecDim p6( x, x, x ); weightsAndPoints_[5] = std::make_pair( w, p6 );
            const VecDim p7( y, x, x ); weightsAndPoints_[6] = std::make_pair( w, p7 );
            const VecDim p8( x, y, x ); weightsAndPoints_[7] = std::make_pair( w, p8 );
            const VecDim p9( x, x, y ); weightsAndPoints_[8] = std::make_pair( w, p9 );

            x = 0.066550153573664, y = 0.433449846426336, w = 0.010949141561386;
            const VecDim p10( y, y, x ); weightsAndPoints_[ 9] = std::make_pair( w, p10 );
            const VecDim p11( y, x, x ); weightsAndPoints_[10] = std::make_pair( w, p11 );
            const VecDim p12( x, y, x ); weightsAndPoints_[11] = std::make_pair( w, p12 );
            const VecDim p13( x, x, y ); weightsAndPoints_[12] = std::make_pair( w, p13 );
            const VecDim p14( y, x, y ); weightsAndPoints_[13] = std::make_pair( w, p14 );
            const VecDim p15( x, y, y ); weightsAndPoints_[14] = std::make_pair( w, p15 );
        }

    }
}
#endif
