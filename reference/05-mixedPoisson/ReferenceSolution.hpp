//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   05-mixedPoisson/ReferenceSolution.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef mixedpoisson_referencesolution_hpp
#define mixedpoisson_referencesolution_hpp

//------------------------------------------------------------------------------
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace ref05{
    template<unsigned DIM> class ReferenceSolution;
}

//------------------------------------------------------------------------------
/** Constructed analytic reference solution.
 *  By inserting a prescribed function \f$ u(x) \f$ into the boundary value
 *  problem, the Dirichlet and Neumann data and the forcing function are
 *  obtained. Here a function of the type
 *  \f[
 *       u(x) = \prod_{d=1}^{DIM} \cos (\alpha_d x_d )
 *  \f]
 *  is used, that is a tensor product of cosine functions.
 *
 *  The image shows a two-dimensional function of this type (it is actually the
 *  FE solution on a 20x20 grid, but visually indistinguishable from the exact
 *  solution). Here the parameters are \f$ \alpha_1 = 3 \f$ and
 *  \f$ \alpha_2 = 4 \f$.
 *
 *  \image html solution.png "2D solution"
 *
 *
 *  \tparam DIM Spatial dimension of the problem
 */
template<unsigned DIM>
class ref05::ReferenceSolution
{
public:
    //! Template parameter: the spatial dimension
    static const unsigned dim = DIM;
    
    //! @name Linear Algebra types
    //@{
    typedef typename base::Vector<1  >::Type   VecDoF;
    typedef typename base::Vector<dim>::Type   VecDim;
    typedef typename base::Matrix<dim,1>::Type GradType;
    //@}

    //! Constructor with solution parameters
    ReferenceSolution( const double alpha1,
                       const double alpha2 = 0.,
                       const double alpha3 = 0. )
    {
        alpha_[0] = alpha1;
        alpha_[1] = alpha2;
        alpha_[2] = alpha3;
    }

    //! Evaluate the solution
    VecDoF evaluate( const VecDim& x ) const
    {
        VecDoF result;
        result[0] = 1;
        for ( unsigned d = 0; d < dim; d++ )
            result *= std::cos( alpha_[d] * x[d] );

        return result;
    }

    //! Evaluate the gradient
    GradType evaluateGradient( const VecDim& x ) const
    {
        GradType result;

        for ( unsigned d = 0; d < dim; d++ ) {
            result[d] = -alpha_[d] * std::sin( alpha_[d] * x[d] );

            for ( unsigned e = 0; e < dim; e++ ) {
                if ( e != d ) {
                    result[d] *= std::cos( alpha_[e] * x[e] );
                }
            }
        }
        return result; 
    }

    //! Evaluate the Laplace operator applied to the solution
    VecDoF laplacian( const VecDim& x ) const
    {
        double factor = 0.;
        for ( unsigned d = 0; d < dim; d++ )
            factor += -alpha_[d] * alpha_[d];

        return factor * evaluate( x );
    }

private:
    //! Parameters
    base::Vector<3,double>::Type alpha_;
};


#endif
