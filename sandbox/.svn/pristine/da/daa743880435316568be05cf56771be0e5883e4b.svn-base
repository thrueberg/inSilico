#include <base/linearAlgebra.hpp>
#include <mat/Lame.hpp>

//------------------------------------------------------------------------------
/** \ingroup thomas
 *  Deformation of an infinite, linear elastic and homogeneous solid due to
 *  a Dipole load.
 */
base::Matrix<3,3>::Type
dipoleSol( const base::Vector<3>::Type& x,
           const base::Vector<3>::Type& y,
           const double E, const double nu )
{
    // x - evaluation, y - location of dipole
    const double G = mat::Lame::mu( E, nu );
    const base::Vector<3>::Type R = x - y;

    const double r = R.norm();

    const double fac = 1./ (8. * M_PI * G * (nu -1.) * r * r * r);

    base::Matrix<3,3>::Type D;
    for ( unsigned i = 0; i < 3; i++ ) {
        for ( unsigned j = 0; j < 3; j++ ) {
            D(i,j)
                //    D(j,i) 
                = fac * (
                0.5 * (3./r/r * R[j] * R[j] - 1.) * R[i] +
                (i==j? (1. - 2.*nu) * R[j] : 0. ) );
        }
    }

    return D;
}
