//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Stokes.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef fluid_stokes_hpp
#define fluid_stokes_hpp

//------------------------------------------------------------------------------
// fluid includes
#include <fluid/VectorLaplace.hpp>
#include <fluid/StressDivergence.hpp>
#include <fluid/PressureGradient.hpp>
#include <fluid/VelocityDivergence.hpp>

namespace fluid{

    template<typename T>
    class Stokes; //!< Documentation class

}


//------------------------------------------------------------------------------
/** Stokes system.
 *  The Stokes' system is based on the coupled partial differential equations
 *  \f[
 *       -\nabla \cdot \sigma(u,p) = f  \quad and \quad \nabla \cdot u = 0
 *  \f]
 *  which express the equilibrium of the fluid stresses \f$ \sigma \f$ in
 *  function of the velocity \f$ u \f$ and the hydrostatic pressure \f$ p \f$,
 *  and the incompressibility of the fluid leads to the condition that
 *  \f$ u \f$ be solenoidal.
 *
 *  Based on the typical weighted residual approach, this system converts to
 *  \f[
 *      \int_\Omega \sigma(u,p) : \nabla v d x -
 *      \int_\Omega q (\nabla \cdot u) d x -
 *      \int_\Gamma \sigma \cdot n \cdot v d s = \int_\Omega f \cdot v d s
 *  \f]
 *  where \f$ v \f$ is the velocity test field and \f$ q \f$ the pressure
 *  test field. Here, no boundary conditions have been applied.
 *
 *  Moreover, the fluid is assumed to be incompressible Newtonian and thus
 *  obey the material law
 *  \f[
 *      \sigma( u, p ) = \mu (\nabla u + \nabla u^T) - p I
 *                     = \mu \varepsilon(u) - p I
 *  \f]
 *  with the dynamic fluid viscosity \f$ \mu \f$ and the identity
 *  matrix \f$ I \f$.
 *
 *  Based on the incompressibility condition (\f$ \nabla \cdot u = 0\f$),
 *  the stress equilibrium equation can also be written as
 *  \f[
 *      \nabla \cdot \sigma =
 *      \nabla \cdot (\mu \nabla u + \mu \nabla u^T - p I) =
 *      \mu \nabla^2 u + \mu \nabla (\nabla \cdot u) - \nabla \cdot (p I) =
 *      \mu \nabla^2 u - \nabla p
 *  \f]
 *  Hence, just the vector Laplacian is used.
 *
 *  The implementation of the terms of the weighted residual method are
 *  given in the classes
 *  fluid::VectorLaplace or fluid::StressDivergence,
 *  fluid::PressureGradient and fluid::VelocityDivergence.
 *  
 */
template<typename T>
class fluid::Stokes
{
    STATIC_ASSERT_MSG( (sizeof(T) == 0), 
                       "This class exists for documentation purpose only" );
};


#endif

