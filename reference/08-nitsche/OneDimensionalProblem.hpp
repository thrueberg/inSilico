//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   OneDimensionalProblem.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef apps_nitsche_onedimensionalproblem_hpp
#define apps_nitsche_onedimensionalproblem_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace apps{
    namespace nitsche{

        template<unsigned DIM>
        struct OneDimensionalProblem;
    }
}

//------------------------------------------------------------------------------
/**  One-dimensional problem with plane interface.
 *   Sketch:
 *   <pre>
 * 
 *
 *          x=0            x=delta       x=1
 *           |----------------X===========|
 *               Omega_1         Omega_2
 *
 *                 d(x) > 0      d(x) < 0
 *
 *   </pre>
 *
 *  Signed distance function:  d(x)
 *
 *  1D Interface BVP:
 *  \f{eqnarray*}{
 *
 *         -\alpha_i d^2 (u_i) / d x^2 &=& 1   \quad   0   < x_1 < 1   \\
 *                                 u_1 &=& 0   \quad   x_1 = 0         \\
 *                                 u_2 &=& 0   \quad   x_1 = 1         \\
 *  \f}
 *
 *  Interface:
 *                    [u] = u_1 - u_2 = 0      x = delta
 *         alpha_1 u_1' - alpha_2 u_2'= 0      x = delta
 *
 *
 *  Analytical solution:
 *                        u_i = a_i x^2 + b_i x + c_i
 *
 *  Omega_1:  a_1 = -1/2/alpha_1, b_1=(see below), c_1 = 0
 *  Omega_2:  a_2 = -1/2/alpha_2, b_2=(see below), c_2 = -a_2 - b_2
 *
 *  Idea taken from Hansbo & Hansbo, CMAME, 2002.
 */
template<unsigned DIM> 
struct apps::nitsche::OneDimensionalProblem
{
    static const unsigned dim = DIM;

    typedef typename base::Vector<dim>::Type VecDim;

    // at x=delta lies the interface
    static bool interface( const VecDim& x, const double delta, VecDim& xClosest )
    {
        xClosest    = x;
        xClosest[0] = delta;

        //// inclined interface, makes only sense for alpha1==alpha2
        //const double alpha = 0.75;
        //const double s = std::sin(alpha) * (x[0] - delta) + std::cos(alpha) * (x[1]-0.5);
        //xClosest[0] = delta + std::sin(alpha) * s;
        //xClosest[1] = 0.5   + std::cos(alpha) * s;

        //// curved interface, makes only sense for alpha1==alpha2
        //xClosest    = x;
        //xClosest[0] = delta + 0.2 * std::sin( 2. * (x[1] - 0.5) );

        return (x[0] < xClosest[0]);
    }

    // Forcing function:  f=1  everywhere
    static base::Vector<1>::Type forceFun( const VecDim& x )
    {
        base::Vector<1>::Type result;
        result[0] = 1.0;
        return result;
    }

    // Coefficient b_1
    static double bOne( const double delta, const double alpha1, const double alpha2 )
    {
        const double aux1 = 2. * (alpha2 * delta - alpha1 * delta + alpha1);
        const double aux2 = alpha2/alpha1 * delta*delta - delta*delta + 1.0;
        return aux2 / aux1;
    }

    // Coefficient b_2
    static double bTwo( const double delta, const double alpha1, const double alpha2 )
    {
        return alpha1 * bOne( delta, alpha1, alpha2 ) / alpha2;
    }

    // Solution in Omega_1: u_1(x)
    static base::Vector<1>::Type sol1( const VecDim& x,
                                       const double delta, const double alpha1,
                                       const double alpha2 )
    {
        const double a1 = -0.5 / alpha1;
        const double b1 = bOne( delta, alpha1, alpha2 );
        const double c1 = 0.0;
        
        base::Vector<1>::Type result;
        VecDim dummy;
        result[0] = a1 * x[0] * x[0] + b1 * x[0] + c1;
        return result;
    }

    // Solution in Omega_2: u_2(x)
    static base::Vector<1>::Type sol2( const VecDim& x,
                                       const double delta, const double alpha1,
                                       const double alpha2 )
    {
        const double a2 = -0.5 / alpha2;
        const double b2 = bTwo( delta, alpha1, alpha2 );
        const double c2 = -a2 - b2;

        base::Vector<1>::Type result;
        VecDim dummy;
        result[0] = a2 * x[0] * x[0] + b2 * x[0] + c2;
        return result;
    }

    // Combined solution
    static base::Vector<1>::Type sol( const VecDim& x,
                                      const double delta, const double alpha1,
                                      const double alpha2 )
    {
        VecDim dummy;
        return (interface(x, delta, dummy) ?
                sol1( x, delta, alpha1, alpha2 ) :
                sol2( x, delta, alpha1, alpha2 ) );
    }

    // Dirichlet function on boundary (zero anyway)
    static base::Vector<1>::Type dirichlet( const VecDim& x,
                                            const double delta,
                                            const double alpha1, const double alpha2 )
    {
        VecDim dummy;
        base::Vector<1>::Type result = ( interface( x, delta, dummy ) ?
                                         sol1( x, delta, alpha1, alpha2 ) :
                                         sol2( x, delta, alpha1, alpha2 ) );
        return result;
    }

    // description of the boundary
    static bool dirichletBoundary( const VecDim& x )
    {
        const double eps = 1.e-8;

        if ( ( std::abs( x[0] - 0.0 ) < eps ) or
             ( std::abs( x[0] - 1.0 ) < eps ) ) return true;

        return false;
    }

    // Derivative of solution in Omega_1: u_1'(x)
    static typename base::Matrix<dim,1>::Type sol1Deriv( const VecDim& x,
                                                         const double delta,
                                                         const double alpha1,
                                                         const double alpha2 )
    {
        const double a1 = -0.5 / alpha1;
        const double b1 = bOne( delta, alpha1, alpha2 );
        
        typename base::Matrix<dim,1>::Type result = base::constantMatrix<dim,1>(0.);
        result(0,0) = 2. * a1 * x[0] + b1;
        return result;
    }

    // Derivative of solution in Omega_2: u_2'(x)
    static typename base::Matrix<dim,1>::Type sol2Deriv( const VecDim& x,
                                                         const double delta,
                                                         const double alpha1,
                                                         const double alpha2 )
    {
        const double a2 = -0.5 / alpha2;
        const double b2 = bTwo( delta, alpha1, alpha2 );

        typename base::Matrix<dim,1>::Type result = base::constantMatrix<dim,1>(0.);
        result(0,0) = 2. * a2 * x[0] + b2;
        return result;
    }

    // combined derivative
    static typename base::Matrix<dim,1>::Type solDeriv( const VecDim& x,
                                                        const double delta,
                                                        const double alpha1,
                                                        const double alpha2 )
    {
        VecDim dummy;
        return (interface(x,delta,dummy) ?
                sol1Deriv( x, delta, alpha1, alpha2 ) :
                sol2Deriv( x, delta, alpha1, alpha2 ) );
    }


};

#endif
