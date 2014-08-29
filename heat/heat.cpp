#error Documentation only

//------------------------------------------------------------------------------
/** \namespace heat
 *  Physical model of heat equation or similar.
 *  Consider the scalar differential equation
 *  \f[
 *      u_t - \nabla \cdot (\kappa \nabla u) = f
 *  \f]
 *  where \f$ u \f$ is the temperature or concentration (or whatever the model
 *  stands for), \f$ \kappa \f$ is a \f$ d\times d\f$-tensor, possibly
 *  dependent on the location \f$ x \f$ and the solution \f$ u \f$.
 *  In this module, the element stiffness matrix and residual forces are
 *  implemented which are needed to solve above equation in combination with
 *  time integration method of the module base::time.
 *  Optionally, above equation can be augmented by a convection term
 *  \f[
 *        \vec{v} \cdot (\nabla u)
 *  \f]
 *  
 */
