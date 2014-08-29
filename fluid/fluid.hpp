#error Documentation only

//------------------------------------------------------------------------------
/** \namespace fluid
 *  Physical model of incompressible viscous fluid flow.
 *  The Navier-Stokes equations are
 *  \f{eqnarray*}{
 *   \rho u_t + u \cdot \nabla u - \nabla \cdot \sigma(u,p) &=& f \cr
 *                                  \nabla \cdot u          &=& 0
 *  \f}
 *  with the fluid velocity \f$ u \f$ and the pressure \f$ p \f$.
 *  The first term is treated by a time integration from base::time,
 *  the second is fluid::Convection and the remaining terms are a
 *  Stokes system.
 */

/** \namespace fluid::gls
 *  Stabilisation by means of a Galerkin Least Squares method.
 */
