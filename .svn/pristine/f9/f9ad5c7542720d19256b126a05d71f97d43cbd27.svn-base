//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FundamentalSolution.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_auxi_fundamentalsolution_hpp
#define base_auxi_fundamentalsolution_hpp

//------------------------------------------------------------------------------
// base  includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace auxi{

        template<unsigned DIM> class FundSolLaplace;
        template<unsigned DIM> class FundSolElastoStatic;

        namespace detail_{

            template<unsigned DIM> struct LaplaceKernel;

            template<>
            struct LaplaceKernel<2>
            {
                static double apply( const double dist )
                {
                    return std::log( 1./dist );
                }
            };

            template<>
            struct LaplaceKernel<3>
            {
                static double apply( const double dist )
                {
                    return 1./dist;
                }
            };

        }

    }
}

//------------------------------------------------------------------------------
/** Fundamental solution of the Laplace operator.
 *  The d-dimensional Laplace operator reads
 *  \f[
 *       L u = - \Delta u = - \sum_{i=1}^d \frac{\partial^2 u}{\partial x_i^2}
 *  \f]
 *  Its fundamental solution, \f$ u^*(x,y) \f$ is defined such that
 *  \f[
 *        L_y u^*(x,y) = \delta(y-x)
 *  \f]
 *  That is the application of the operator acting on the \f$ y \f$ coordinates
 *  gives a point source of unit strength at the location \f$ x \f$.
 *
 *  The fundamental solution \f$ u^* \f$ can be written in the concise notation
 *  \f[
 *        u^*(x,y) = \frac{1}{2(d-1) \pi} \gamma(x,y)
 *  \f]
 *  with the dimension dependent kernel function \f$ \gamma \f$ defined as
 *  \f[
 *        \gamma_{2D} (x,y) = \log \frac{1}{|y-x|}
 *  \f]
 *  in 2D and
 *  \f[
 *        \gamma_{3D} (x,y) = \frac{1}{|y-x|}
 *  \f]
 *  in 3D.
 *
 *  \tparam DIM  Spatial dimension of the problem.
 */
template<unsigned DIM>
class base::auxi::FundSolLaplace
{
public:
    static const unsigned dim     = DIM;
    static const unsigned doFSize = 1;

    typedef typename base::Vector<dim,    double>::Type     VecDim;
    typedef typename base::Vector<doFSize,double>::Type     VecDoF;
    typedef typename base::Matrix<dim,doFSize,double>::Type Grad;

    //! Evaluate the fundamental solution for arguments x and y
    VecDoF fun( const VecDim& x, const VecDim& y ) const
    {
        const double dist = base::norm(y - x);
        const double factor = 1./(2.*(dim-1.)*M_PI);
        VecDoF result;
        result[0] = factor * detail_::LaplaceKernel<dim>::apply( dist );
        return result;
    }

    //! Evaluate the gradient of the fundamental solution for given x and y
    Grad grad( const VecDim& x, const VecDim& y ) const
    {
        const double dist = base::norm(y - x);
        const double factor =
            1./(2.*(dim-1.)*M_PI) / (base::Power<dim>::apply( dist ) );

        Grad result;
        for ( unsigned d = 0; d < dim; d++ )
            result(d,0) = factor * (y[d] - x[d]);

        return result;
    }

    //! Evaluate the co-normal derivative for given x, y and normal vector
    VecDoF coNormal( const VecDim& x, const VecDim& y,
                     const VecDim& normal ) const
    {
        const Grad grad = this -> grad( x, y );
        return grad.transpose() * normal;
    }
    
};

//==============================================================================

//------------------------------------------------------------------------------
/** Fundamental solution of the elastostatic operator.
 *  Story to be told.
 */
template<unsigned DIM>
class base::auxi::FundSolElastoStatic
{
public:
    static const unsigned dim     = DIM;
    static const unsigned doFSize = dim;

    typedef typename base::Vector<dim,    double>::Type     VecDim;
    typedef typename base::Vector<doFSize,double>::Type     VecDoF;
    typedef typename base::Matrix<dim,doFSize,double>::Type Grad;

    FundSolElastoStatic( const double lambda, const double mu )
        : lambda_( lambda ), mu_( mu ) { }
    

    VecDoF fun( const VecDim& x, const VecDim& y, const VecDim& dir ) const
    {
        const double factor =
            1. /(4. * (dim-1) * M_PI ) *
            (lambda_ + mu_) / (lambda_ + 2.* mu_) / mu_;

        const double dist   = base::norm(y - x);

        typedef Grad MatDimDim;

        MatDimDim U;
        for ( unsigned i = 0; i < dim; i++ ) {
            for ( unsigned j = 0; j < dim; j++ ) {
                U(i,j) = factor *
                    ( (i==j ?
                       (lambda_+3.*mu_)/(lambda_+mu_) *
                       detail_::LaplaceKernel<dim>::apply( dist ) :
                       0.
                        ) +
                      ( (y[i]-x[i])*(y[j]-x[j]) / base::Power<dim>::apply( dist ) )
                        );
            }
        }

        VecDoF result;
        result.noalias() = U * dir;
        return result;
    }

    Grad grad( const VecDim& x, const VecDim& y, const VecDim& dir ) const
    {
        /* empty */
        VERIFY_MSG( false, "Not implemented" );
    }

    VecDoF coNormal( const VecDim& x, const VecDim& y, const VecDim& dir,
                     const VecDim& normal ) const
    {
        /* empty */
        VERIFY_MSG( false, "Not implemented" );
    }

private:
    //! @name Lame parameters
    //@{
    const double lambda_;
    const double mu_;
    //@}
    
};

#endif
