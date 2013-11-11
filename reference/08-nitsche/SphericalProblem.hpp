//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SphericalProblem.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef apps_nitsche_sphericalproblem_hpp
#define apps_nitsche_sphericalproblem_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace apps{
    namespace nitsche{

        template<unsigned DIM>
        struct SphericalProblem;
    }
}


    
//------------------------------------------------------------------------------
/** Square domain with a circular inclusion with different material parameter.
 *  In spherical coordinates the problem reduces to
 *  \f[
 *         - \alpha_i u_i^{\prime\prime}(r) = -2d
 *  \f]
 *  and the material parameter \f$ \alpha_i \f$ are potentially different for
 *  the domains \f$ \Omega_1 = \{ |x| < R \} \f$ and
 *  \f$ \Omega_2 = [0,1]^2 \setminus \Omega_1 \f$.
 *  The analytic solution becomes
 *  \f[
 *     u_i(r) = \frac{r^2}{\alpha_i}
 *                    - \frac{R^2}{\alpha_i} + \frac{R^2}{\alpha_1}
 *  \f]
 *  
 *
 *  Idea taken from Hansbo & Hansbo, CMAME, 2002.
 */
template<unsigned DIM> 
struct apps::nitsche::SphericalProblem
{
    static const unsigned dim = DIM;

    typedef typename base::Vector<dim>::Type VecDim;

    static VecDim makeCentre() { return base::constantVector<dim>( 0.5 ); }

    // at x=delta lies the interface
    static bool interface( const VecDim& x, const double radius, VecDim& xClosest )
    {
        const VecDim centre = base::constantVector<dim>( 0.5 );
        const VecDim diff   = x - centre;
        const double r      = diff.norm();

        // catch the centre of the sphere
        if ( r    > 1.e-10 )
            xClosest    = centre + (radius / r) * diff;
        else 
            xClosest    = centre + (radius /centre.norm()) * centre;
        
        return ( r < radius );
    }

    // description of the boundary
    static bool dirichletBoundary( const VecDim& x )
    {
        const double eps = 1.e-8;

        for ( unsigned d = 0; d < dim; d++ ) {
            
            if ( ( std::abs( x[d] - 0.0 ) < eps ) or
                 ( std::abs( x[d] - 1.0 ) < eps ) ) return true;
            
        }

        return false;
    }

    // Forcing function:  f=1  everywhere
    static base::Vector<1>::Type forceFun( const VecDim& x )
    {
        base::Vector<1>::Type result;
        result[0] = -2.0 * static_cast<double>( dim );
        return result;
    }

    // Solution in Omega_1: u_1(x)
    static base::Vector<1>::Type sol1( const VecDim& x,
                                       const double radius, const double alpha1,
                                       const double alpha2 )
    {
        const VecDim centre = base::constantVector<dim>( 0.5 );
        const VecDim diff   = x - centre;
        const double r      = diff.norm();

        base::Vector<1>::Type result;
        result[0] = r * r / alpha1;
        
        return result;
    }

    // Solution in Omega_2: u_2(x)
    static base::Vector<1>::Type sol2( const VecDim& x,
                                       const double radius, const double alpha1,
                                       const double alpha2 )
    {
        const VecDim centre = base::constantVector<dim>( 0.5 );
        const VecDim diff   = x - centre;
        const double r      = diff.norm();

        base::Vector<1>::Type result;
        result[0] = r*r / alpha2 - radius*radius/alpha2 + radius*radius/alpha1;
        return result;
    }

    static base::Vector<1>::Type sol( const VecDim& x,
                                      const double radius, const double alpha1,
                                      const double alpha2 )
    {
        VecDim dummy;
        return (interface( x, radius, dummy ) ?
                sol1( x, radius, alpha1, alpha2 ) :
                sol2( x, radius, alpha1, alpha2 ) );
    }


    // Dirichlet function on boundary (zero anyway)
    static base::Vector<1>::Type dirichlet( const VecDim& x, const double radius,
                                            const double alpha1, const double alpha2 )
    {
        return sol( x, radius, alpha1, alpha2 );
    }


    // Derivative of solution in Omega_1: u_1'(x)
    static typename base::Matrix<dim,1>::Type sol1Deriv( const VecDim& x,
                                                         const double radius,
                                                         const double alpha1,
                                                         const double alpha2 )
    {
        const VecDim diff   = x - base::constantVector<dim>( 0.5 );
        
        typename base::Matrix<dim,1>::Type result = base::constantMatrix<dim,1>(0.);
        for ( unsigned d = 0; d < dim; d++ )
            result(d,0) = 2.0 / alpha1 * diff[d];

        return result;
    }

    // Derivative of solution in Omega_2: u_2'(x)
    static typename base::Matrix<dim,1>::Type sol2Deriv( const VecDim& x,
                                                         const double radius,
                                                         const double alpha1,
                                                         const double alpha2 )
    {
        const VecDim diff   = x - base::constantVector<dim>( 0.5 );
        
        typename base::Matrix<dim,1>::Type result = base::constantMatrix<dim,1>(0.);
        for ( unsigned d = 0; d < dim; d++ )
            result(d,0) = 2.0 / alpha2 * diff[d];

        return result;
    }

    static typename base::Matrix<dim,1>::Type solDeriv( const VecDim& x,
                                                        const double radius,
                                                        const double alpha1,
                                                        const double alpha2 )
    {
        VecDim dummy;
        return (interface( x, radius, dummy ) ?
                sol1Deriv( x, radius, alpha1, alpha2 ) :
                sol2Deriv( x, radius, alpha1, alpha2 ) );
    }

};

#endif
