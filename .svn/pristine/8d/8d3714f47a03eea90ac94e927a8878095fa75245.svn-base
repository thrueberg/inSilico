//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   GaussLegendre.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_quad_gausslegendre_hpp
#define base_quad_gausslegendre_hpp

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
        
        template<unsigned DEGREE> class GaussLegendre;

        namespace detail_{

            //! Provide the number of Gauss points, given the polynomial degree
            template<unsigned DEGREE>
            struct NumGaussPoints
            {
                //! Note '+2' is necessary for correct rounding
                static const unsigned value = (DEGREE + 2) / 2;
            };

            //! Provide the N-th point rule in an extra structure
            template<unsigned NUMPOINTS>
            struct GaussLegendreRule
            {
                boost::array< std::pair<double,double>, NUMPOINTS > data;
                GaussLegendreRule();
            };
            
        }
        
    }
}

//------------------------------------------------------------------------------
/** \brief Gauss-Legendre quadrature rules for given polynomial degree
 *  \tparam DEGREE Polynomial degree to be integrated exactly
 */
template<unsigned DEGREE> 
class base::quad::GaussLegendre
    : public boost::noncopyable
{
    //! Sanity check DEGREE > 0
    STATIC_ASSERT_MSG( (DEGREE > 0), "Degree has to be larger than zero" );
    
public:
    //! Template parameter: polynomial degree to be integrated exactly
    static const unsigned degree = DEGREE;

    //! Number of points of the rule
    static const unsigned numPoints =
        base::quad::detail_::NumGaussPoints<degree>::value;

    //! The local coordinate dimension
    static const unsigned dim = 1;
    
    //! Type of coordinate-vector
    typedef typename base::Vector<dim>::Type VecDim;

    //! Type of iterator for external access
    typedef typename boost::array< std::pair<double, VecDim>,
                                   numPoints>::const_iterator Iter;

    //! Specialized constructor for different number of Points
    GaussLegendre( )
    {
        //! Collect from data base
        detail_::GaussLegendreRule<numPoints> glr;
        for ( unsigned g = 0; g < numPoints; g ++ ) {
            const double weight = glr.data[g].first;
            const VecDim   point  = base::constantVector<1>( glr.data[g].second );
            
            weightsAndPoints_[g] = std::make_pair( weight, point );
        }
    }

    //! Begin of array iterator
    Iter begin() const  { return weightsAndPoints_.begin(); }
    //! End of array iterator
    Iter end()   const  { return weightsAndPoints_.end();   }

private:
    //! Pair of weights and points representing the NPTS-point quadrature rule
    boost::array<std::pair<double,VecDim>, numPoints> weightsAndPoints_;
};


//------------------------------------------------------------------------------
namespace base{
    namespace quad{
        namespace detail_{
            //------------------------------------------------------------------
            /* The implementation of Gauss-Legendre quadrature rules for
             *        the reference interval (0,1)
             */

            
            //! 1-point rule
            template<> GaussLegendreRule<1>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 1, 0.5);
            }

            //! 2-point rule
            template<> GaussLegendreRule<2>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.5, 0.788675134594813);
                data[1] =  std::make_pair( 0.5, 0.211324865405187);
            }

            //! 3-point rule
            template<> GaussLegendreRule<3>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.277777777777777, 0.887298334620741);
                data[1] =  std::make_pair( 0.444444444444444, 0.5);
                data[2] =  std::make_pair( 0.277777777777777, 0.112701665379259);
            }

            //! 4-point rule
            template<> GaussLegendreRule<4>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.173927422568727, 0.930568155797026);
                data[1] =  std::make_pair( 0.326072577431273, 0.669990521792428);
                data[2] =  std::make_pair( 0.326072577431273, 0.330009478207572);
                data[3] =  std::make_pair( 0.173927422568727, 0.069431844202974);
            }

            //! 5-point rule
            template<> GaussLegendreRule<5>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.118463442528095, 0.953089922969332);
                data[1] =  std::make_pair( 0.239314335249683, 0.769234655052841);
                data[2] =  std::make_pair( 0.284444444444444, 0.5);
                data[3] =  std::make_pair( 0.239314335249683, 0.230765344947159);
                data[4] =  std::make_pair( 0.118463442528095, 0.046910077030668);
            }

            //! 6-point rule
            template<> GaussLegendreRule<6>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.085662246189585, 0.966234757101576);
                data[1] =  std::make_pair( 0.180380786524069, 0.830604693233132);
                data[2] =  std::make_pair( 0.233956967286345, 0.619309593041598);
                data[3] =  std::make_pair( 0.233956967286345, 0.380690406958402);
                data[4] =  std::make_pair( 0.180380786524069, 0.169395306766868);
                data[5] =  std::make_pair( 0.085662246189585, 0.033765242898424);
            }

            //! 7-point rule
            template<> GaussLegendreRule<7>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.064742483084435, 0.974553956171379);
                data[1] =  std::make_pair( 0.139852695744638, 0.870765592799697);
                data[2] =  std::make_pair( 0.190915025252559, 0.702922575688699);
                data[3] =  std::make_pair( 0.208979591836735, 0.5);
                data[4] =  std::make_pair( 0.190915025252559, 0.297077424311301);
                data[5] =  std::make_pair( 0.139852695744638, 0.129234407200303);
                data[6] =  std::make_pair( 0.064742483084435, 0.025446043828621);
            }

            //! 8-point rule
            template<> GaussLegendreRule<8>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.0506142681451885, 0.980144928248768);
                data[1] =  std::make_pair( 0.111190517226687, 0.898333238706813);
                data[2] =  std::make_pair( 0.156853322938943, 0.762766204958164);
                data[3] =  std::make_pair( 0.181341891689181, 0.591717321247825);
                data[4] =  std::make_pair( 0.181341891689181, 0.408282678752175);
                data[5] =  std::make_pair( 0.156853322938943, 0.237233795041836);
                data[6] =  std::make_pair( 0.111190517226687, 0.101666761293187);
                data[7] =  std::make_pair( 0.0506142681451885, 0.019855071751232);
            }

            //! 9-point rule
            template<> GaussLegendreRule<9>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.0406371941807875, 0.984080119753813);
                data[1] =  std::make_pair( 0.0903240803474285, 0.918015553663318);
                data[2] =  std::make_pair( 0.130305348201468, 0.806685716350295);
                data[3] =  std::make_pair( 0.156173538520001, 0.662126711701905);
                data[4] =  std::make_pair( 0.16511967750063, 0.5);
                data[5] =  std::make_pair( 0.156173538520001, 0.337873288298095);
                data[6] =  std::make_pair( 0.130305348201468, 0.193314283649705);
                data[7] =  std::make_pair( 0.0903240803474285, 0.081984446336682);
                data[8] =  std::make_pair( 0.0406371941807875, 0.015919880246187);
            }

            //! 10-point rule
            template<> GaussLegendreRule<10>::GaussLegendreRule()
            {
                data[0] =  std::make_pair( 0.033335672154344, 0.986953264258586);
                data[1] =  std::make_pair( 0.0747256745752905, 0.932531683344493);
                data[2] =  std::make_pair( 0.109543181257991, 0.839704784149512);
                data[3] =  std::make_pair( 0.134633359654998, 0.716697697064623);
                data[4] =  std::make_pair( 0.147762112357376, 0.574437169490815);
                data[5] =  std::make_pair( 0.147762112357376, 0.425562830509185);
                data[6] =  std::make_pair( 0.134633359654998, 0.283302302935377);
                data[7] =  std::make_pair( 0.109543181257991, 0.160295215850488);
                data[8] =  std::make_pair( 0.0747256745752905, 0.0674683166555075);
                data[9] =  std::make_pair( 0.033335672154344, 0.013046735741414);
            }
            
        } // end namespace detail_
    } // end namespace quad
} // end namespace base
//------------------------------------------------------------------------------
#endif
