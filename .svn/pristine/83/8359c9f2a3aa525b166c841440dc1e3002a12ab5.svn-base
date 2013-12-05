//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   FenicsTest.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef mat_thermal_fenicstest_hpp
#define mat_thermal_fenicstest_hpp

#include <cmath>

//------------------------------------------------------------------------------
namespace mat{
    namespace thermal{

        class FenicsTest;
    }
}

//------------------------------------------------------------------------------
/**  Solution-dependent material bevhaviour for testing a nonlinear solver.
 *   The idea is taking from Fenics
 *   (see http://fenicsproject.org/documentation/tutorial/nonlinear.html)
 *   and has the advantage that for a uni-dimensional setup an analytical
 *   solution is available. For instance the BVP
 *   \f[
 *          -( K(u) u^\prime)^\prime = 0,   0 < x < 1
 *   \f]
 *   with BCs  \f$ u = 0 \f$ at \f$ x = 0 \f$ and \f$ u = 1 \f$ at
 *   \f$ x = 1 \f$ and the here implemented material behaviour
 *   \f[
 *        K(u) = (1 + u)^m
 *   \f]
 *   has the analytic solution
 *   \f[
 *        u = ( (2^{m+1} - 1) x - 1)^{1/(m+1)}
 *   \f]
 *   
 */
class mat::thermal::FenicsTest
{
public:
    // construct with exponent \f$ m \f$
    FenicsTest( const double m )
        : m_( m )
    { }

    typedef mat::Vector Vector;
    typedef mat::Tensor Tensor;

    void conductivity( const double u, const Vector& gradU,
                       Tensor& K ) const
    {
        const double factor = std::pow( 1. + u, m_ );
        K = factor * Tensor::Identity();
    }

    void conductivityGradient( const double u, const Vector& gradU,
                               Tensor& gradK ) const
    {
        const double factor = m_ * std::pow( 1. + u, m_ - 1.0 );
        gradK = factor * Tensor::Identity();
    }

private:
    const double m_;
};

#endif
