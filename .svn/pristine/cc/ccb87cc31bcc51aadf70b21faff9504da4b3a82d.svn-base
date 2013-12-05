//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   GaussTriangle.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_quad_gausstriangle_hpp
#define base_quad_gausstriangle_hpp

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
        class GaussTriangle;

        // Declare a lookup table for the require number of points
        namespace detail_{

            template<unsigned DEGREE> struct NumTriaPoints;
            template<> struct NumTriaPoints<1>{ static const unsigned value =  1; };
            template<> struct NumTriaPoints<2>{ static const unsigned value =  3; };
            template<> struct NumTriaPoints<3>{ static const unsigned value =  4; };
            template<> struct NumTriaPoints<4>{ static const unsigned value =  6; };
            template<> struct NumTriaPoints<5>{ static const unsigned value =  7; };
            template<> struct NumTriaPoints<6>{ static const unsigned value = 13; };
            template<> struct NumTriaPoints<7>{ static const unsigned value = 13; };
            template<> struct NumTriaPoints<8>{ static const unsigned value = 19; };
            template<> struct NumTriaPoints<9>{ static const unsigned value = 19; };
        }

    }
}

//------------------------------------------------------------------------------
/** \brief Gauss quadrature rules for triangles with given polynomial degree
 *  \tparam DEGREE Polynomial degree to be integrated exactly
 */
template<unsigned DEGREE> 
class base::quad::GaussTriangle
    : public boost::noncopyable
{
    //! Sanity check DEGREE > 0
    STATIC_ASSERT_MSG( (DEGREE > 0),
                       "Degree has to be larger than zero" );

public:
    //! Template parameter: polynomial degree to be integrated exactly
    static const unsigned degree = DEGREE;

    //! Deduce number of quadrature points from degree and lookup table
    static const unsigned numPoints =
        base::quad::detail_::NumTriaPoints<degree>::value;

    //! The local coordinate dimension
    static const unsigned dim = 2;
    
    //! Type of coordinate-vector
    typedef typename base::Vector<dim>::Type VecDim;

    //! Type of iterator for external access
    typedef typename boost::array< std::pair<double, VecDim>,
                                   numPoints>::const_iterator Iter;
    
    //! Specialized constructor for different polynomial degrees
    GaussTriangle();
    
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
        template<> GaussTriangle<1>::GaussTriangle()
        {
            const VecDim point( 0.333333333333333, 0.333333333333333 );
            weightsAndPoints_[0] = std::make_pair( 0.5, point );
        }

        //----------------------------------------------------------------------
        //! 3-point rule, exact for quadratic polynomials
        template<> GaussTriangle<2>::GaussTriangle()
        {
            const VecDim point1( 0.666666666666667, 0.166666666666667 );
            weightsAndPoints_[0] = std::make_pair( 0.166666666666666, point1 );

            const VecDim point2( 0.166666666666667, 0.666666666666667 );
            weightsAndPoints_[1] = std::make_pair( 0.166666666666666, point2 );

            const VecDim point3( 0.166666666666667, 0.166666666666667 );
            weightsAndPoints_[2] = std::make_pair( 0.166666666666666, point3 );
        }

        //----------------------------------------------------------------------
        //! 4-point rule, exact for cubic polynomials
        template<> GaussTriangle<3>::GaussTriangle()
        {
            const VecDim point1( 0.333333333333333, 0.333333333333333 );
            weightsAndPoints_[0] = std::make_pair( -0.28125, point1 );

            const VecDim point2( 0.6, 0.2 ); 
            weightsAndPoints_[1] = std::make_pair( 0.260416666666667, point2 );

            const VecDim point3( 0.2, 0.6 ); 
            weightsAndPoints_[2] = std::make_pair( 0.260416666666667, point3 );

            const VecDim point4( 0.2, 0.2 ); 
            weightsAndPoints_[3] = std::make_pair( 0.260416666666667, point4 );
        }

        //----------------------------------------------------------------------
        //! 6-point rule, exact for fourth-order polynomials
        template<> GaussTriangle<4>::GaussTriangle()
        {
            const VecDim point1( 0.10810301816807, 0.445948490915965 ); 
            weightsAndPoints_[0] = std::make_pair( 0.111690794839005, point1 );

            const VecDim point2( 0.816847572980459, 0.091576213509771 ); 
            weightsAndPoints_[1] = std::make_pair( 0.054975871827661, point2 );

            const VecDim point3( 0.445948490915965, 0.10810301816807 );
            weightsAndPoints_[2] = std::make_pair( 0.111690794839005, point3 );

            const VecDim point4( 0.445948490915965, 0.445948490915965 ); 
            weightsAndPoints_[3] = std::make_pair( 0.111690794839005, point4 );

            const VecDim point5( 0.091576213509771, 0.816847572980459 );
            weightsAndPoints_[4] = std::make_pair( 0.054975871827661, point5 );

            const VecDim point6( 0.091576213509771, 0.091576213509771 );
            weightsAndPoints_[5] = std::make_pair( 0.054975871827661, point6 );
        }

        //----------------------------------------------------------------------
        //! 7-point rule, exact for fifth-order polynomials
        template<> GaussTriangle<5>::GaussTriangle()
        {
            const VecDim point1( 0.333333333333333, 0.333333333333333 ); 
            weightsAndPoints_[0] = std::make_pair( 0.1125,             point1 );

            const VecDim point2( 0.05971587178977,  0.470142064105115 ); 
            weightsAndPoints_[1] = std::make_pair( 0.066197076394253,  point2 );

            const VecDim point3( 0.797426985353087, 0.101286507323456 ); 
            weightsAndPoints_[2] = std::make_pair( 0.0629695902724135, point3 );

            const VecDim point4( 0.470142064105115, 0.05971587178977 ); 
            weightsAndPoints_[3] = std::make_pair( 0.066197076394253,  point4 );

            const VecDim point5( 0.470142064105115, 0.470142064105115 ); 
            weightsAndPoints_[4] = std::make_pair( 0.066197076394253,  point5 );

            const VecDim point6( 0.101286507323456, 0.797426985353087 ); 
            weightsAndPoints_[5] = std::make_pair( 0.0629695902724135, point6 );

            const VecDim point7( 0.101286507323456, 0.101286507323456 ); 
            weightsAndPoints_[6] = std::make_pair( 0.0629695902724135, point7 );

        }

        //----------------------------------------------------------------------
        //! 13-point rule, exact for seventh-order polynomials
        template<> GaussTriangle<7>::GaussTriangle()
        {
            const double a0 =  0.333333333333333;
            const double w0 = -0.074785022233835;

            const VecDim point1( a0, a0 );
            weightsAndPoints_[0] = std::make_pair( w0, point1 );

            const double a1 = 0.260345966079038;
            const double b1 = 0.479308067841923;
            const double w1 = 0.087807628716602;

            const VecDim point2( a1, a1 ); 
            weightsAndPoints_[1] = std::make_pair( w1, point2 );

            const VecDim point3( a1, b1 ); 
            weightsAndPoints_[2] = std::make_pair( w1, point3 );

            const VecDim point4( b1, a1 );
            weightsAndPoints_[3] = std::make_pair( w1, point4 );

            const double a2 = 0.065130102902216;
            const double b2 = 0.869739794195568;
            const double w2 = 0.02667361780442;

            const VecDim point5( a2, a2 ); 
            weightsAndPoints_[4] = std::make_pair( w2, point5 );

            const VecDim point6( a2, b2 ); 
            weightsAndPoints_[5] = std::make_pair( w2, point6 );

            const VecDim point7( b2, a2 ); 
            weightsAndPoints_[6] = std::make_pair( w2, point7 );

            const double a3 = 0.048690315425316;
            const double b3 = 0.312865496004875;
            const double c3 = 0.638444188569809;
            const double w3 = 0.038556880445128;

            const VecDim point8( a3, b3 ); 
            weightsAndPoints_[7] = std::make_pair( w3, point8 );

            const VecDim point9( b3, a3 ); 
            weightsAndPoints_[8] = std::make_pair( w3, point9 );

            const VecDim point10( a3, c3 ); 
            weightsAndPoints_[9] = std::make_pair( w3, point10 );

            const VecDim point11( b3, c3 );
            weightsAndPoints_[10] = std::make_pair( w3, point11 );

            const VecDim point12( c3, a3 );
            weightsAndPoints_[11] = std::make_pair( w3, point12 );

            const VecDim point13( c3, b3 );
            weightsAndPoints_[12] = std::make_pair( w3, point13 );
        }

        //----------------------------------------------------------------------
        //! No sixth-order rule is implemented copy from seventh-order rule
        template<> GaussTriangle<6>::GaussTriangle()
        {
            GaussTriangle<7> gaussSeven;
            std::copy( gaussSeven.begin(), gaussSeven.end(),
                       weightsAndPoints_.begin() );
        }


        //----------------------------------------------------------------------
        //! 19-point rule, exact for ninth-order polynomials
        template<> GaussTriangle<9>::GaussTriangle()
        {
            const double xi[19] = {
                .3333333333333329 , .0206349616025250 , .1258208170141270 ,
                .6235929287619349 , .9105409732110950 , .0368384120547360 ,
                .4896825191987380 , .4896825191987380 , .4370895914929371 ,
                .4370895914929371 , .1882035356190330 , .1882035356190330 ,     
                .0447295133944530 , .0447295133944530 , .0368384120547360 ,
                .2219629891607660 , .7411985987844981 , .2219629891607660 ,
                .7411985987844981  };

            const double eta[19] = {
                .3333333333333329 , .4896825191987380 , .4370895914929371 ,
                .1882035356190330 , .0447295133944530 , .2219629891607660 ,
                .0206349616025250 , .4896825191987380 , .1258208170141270 ,
                .4370895914929371 , .6235929287619350 , .1882035356190330 ,  
                .9105409732110950 , .0447295133944530 , .7411985987844981 ,
                .7411985987844981 , .0368384120547360 , .0368384120547360 ,
                .2219629891607660  };

            const double omega[19] = {
                .0971357962827990 , .0313347002271390 , .0778275410047740 ,
                .0796477389272100 , .0255776756586980 , .0432835393772890 ,
                .0313347002271390 , .0313347002271390 , .0778275410047740 ,
                .0778275410047740 , .0796477389272100 , .0796477389272100 , 
                .0255776756586980 , .0255776756586980 , .0432835393772890 ,
                .0432835393772890 , .0432835393772890 , .0432835393772890 ,
                .0432835393772890  };

            for ( unsigned i = 0; i < 19; i ++ ) {
                const VecDim point( xi[i], eta[i] );
                weightsAndPoints_[i] = std::make_pair( 0.5 * omega[i], point );
            }
            
        }
        
        //----------------------------------------------------------------------
        //! No eigth-order rule is implemented copy from ninth-order rule
        template<> GaussTriangle<8>::GaussTriangle()
        {
            GaussTriangle<9> gaussNine;
            std::copy( gaussNine.begin(), gaussNine.end(),
                       weightsAndPoints_.begin() );
        }


    }
}
#endif
