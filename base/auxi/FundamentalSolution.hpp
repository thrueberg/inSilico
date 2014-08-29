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
            struct LaplaceKernel<1>
            {
                static double apply( const double dist )
                {
                    return -dist;
                }
            };

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

            template<unsigned DIM>
            struct UnitSphere
            {
                static double surface() { return 2.*(DIM-1.)*M_PI; }
            };
            
            template<>
            struct UnitSphere<1>
            {
                static double surface() { return 2.; }
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
 *        u^*(x,y) = \frac{1}{\Gamma_d} \gamma(x,y)
 *  \f]
 *  with the surface of the d-dimensional unit sphere \f$ \Gamma_d \f$ and 
 *  dimension-dependent kernel function \f$ \gamma \f$ defined as
 *  \f[
 *        \gamma_{2D} (x,y) = \log \frac{1}{|y-x|}
 *  \f]
 *  in 2D and
 *  \f[
 *        \gamma_{3D} (x,y) = \frac{1}{|y-x|}
 *  \f]
 *  in 3D.
 *
 *  \note For the computation of the derivatives the roles of \f$ x \f$ and
 *        \f$ y \f$ are not interchangeable. Here, the notation is to have
 *        \f$ x \f$ the location of the unit source (i.e., it becomes a
 *        parameter) and \f$ y \f$ is the location where the function is
 *        evaluated.
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
    VecDoF fun( const VecDim& y, const VecDim& x ) const
    {
        const double dist = base::norm(y - x);
        const double factor = 1./ detail_::UnitSphere<dim>::surface();
        VecDoF result;
        result[0] = factor * detail_::LaplaceKernel<dim>::apply( dist );
        return result;
    }

    //! Evaluate the gradient of the fundamental solution for given x and y
    Grad grad( const VecDim& y, const VecDim& x ) const
    {
        const double dist = base::norm(y - x);
        const double factor =
            -1./detail_::UnitSphere<dim>::surface()/(base::Power<dim>::apply( dist ) );
        
        Grad result;
        for ( unsigned d = 0; d < dim; d++ )
            result(d,0) = factor * (y[d] - x[d]);

        return result;
    }

    //! Evaluate the co-normal derivative for given x, y and normal vector
    VecDoF coNormal( const VecDim& y, const VecDim& x,
                     const VecDim& normal ) const
    {
        const Grad grad = this -> grad( y, x );
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
    typedef          Grad                                   MatDimDim;

    //! Convert and store material parameters
    FundSolElastoStatic( const double lambda, const double mu )
        : G_(  mu ),
          nu_( lambda/2./(lambda+mu) ) { }
          

    //! Kelvin Tensor
    MatDimDim U( const VecDim& y, const VecDim& x ) const
    {
        const double fac1 =
            1./ (4. * detail_::UnitSphere<dim>::surface() * G_ * (1.-nu_));
        const double fac2 = 3. - 4. * nu_;
        const VecDim R      = y - x;
        const double dist   = base::norm(R);

        MatDimDim U;
        for ( unsigned i = 0; i < dim; i++ ) {
            for ( unsigned j = 0; j < dim; j++ ) {
                U(i,j) = fac1 *
                    ( (i==j ? fac2 * detail_::LaplaceKernel<dim>::apply( dist ) : 0. ) +
                      (R[i]*R[j]) / base::Power<dim>::apply( dist )
                        );
            }
        }

        return U;
    }

    //! Traction kernel
    MatDimDim T( const VecDim& y, const VecDim& x, const VecDim& normal ) const
    {
        const double fac1 =
            -1./ (2. * detail_::UnitSphere<dim>::surface() * (1.-nu_));
        const double fac2 = 1. - 2. * nu_;
        const VecDim R    = y - x;
        const double dist = base::norm( R );
        const double rXn  = R.dot( normal );

        MatDimDim T;
        for ( unsigned i = 0; i < dim; i++ ) {
            for ( unsigned j = 0; j < dim; j++ ) {

                T(i,j) = fac1 / base::Power<dim>::apply( dist ) * (
                    fac2 * (i==j ? rXn : (normal[j] * R[i] - normal[i] * R[j])  ) +
                    static_cast<double>( dim ) * (R[i] * R[j] / dist /dist ) * rXn
                    );
                
            }
        }
            
        return T;
    }
    

    VecDoF fun( const VecDim& y, const VecDim& x, const VecDim& dir ) const
    {
        VecDoF result;
        result.noalias() = U( y, x ) * dir;
        return result;
    }

    Grad grad( const VecDim& y, const VecDim& x, const VecDim& dir ) const
    {
        /* empty */
        VERIFY_MSG( false, "Not implemented" );
        return Grad();
    }

    VecDoF coNormal( const VecDim& y, const VecDim& x, const VecDim& dir,
                     const VecDim& normal ) const
    {
        VecDoF result;
        result.noalias() = (this -> T(y,x,normal)) * dir;
        return result;
    }

private:
    //! @name Material parameters
    //@{
    const double G_;
    const double nu_;
    //@}
    
};

#endif
