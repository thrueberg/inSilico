//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   05-mixedPoisson/GivenData.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef mixedpoisson_givendata_hpp
#define mixedpoisson_givendata_hpp

//------------------------------------------------------------------------------
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace ref05{
    template<typename REFSOL> class GivenData;
}

//------------------------------------------------------------------------------
/** Define a boundary value problem for Poisson's equation.
 *  The boundary value problem reads
 *  \f[
 *       -\nabla^2 u = f \quad x \in \Omega
 *  \f]
 *  with Dirichlet boundary condition
 *  \f[
 *          u = \bar{u}  \quad x \in \Gamma_D
 *  \f]
 *  and Neumann boundary condition
 *  \f[
 *         \nabla u \cdot \vec{n} = \bar{t} \quad x \in \Gamma_N
 *  \f]
 *  Here, it is assumed that \f$ x_i = 0 \f$ is the Dirichlet boundary, that is
 *  independent of the spatial dimension and the shape of the domain
 *  \f$ \Omega \f$. Any boundary degrees of freedom fulfilling \f$ x_i = 0 \f$
 *  will be constrained with Dirichlet boundary condition. The remaining
 *  boundary is a Neumann boundary on which the co-normal derivative is
 *  prescribed.
 *
 *  The given data \f$ f \f$, \f$ \bar{u} \f$ and \f$ \bar{u} \f$ are taken from
 *  a given reference solution which provides
 *
 *  -  function evaluation \f$ u(x) \f$ used for the Dirichlet BC
 *     \f$ \bar{u} \f$
 *  -  gradient evaluation \f$ \nabla u (x) \f$ used for the Neumann BC
 *     \f$ \bar{t} \f$
 *  - evaluate of the Laplacian \f$ \nabla^2 u(x) \f$ used for the force term
 *    \f$ f \f$
 *
 *
 *  \tparam REFSOL Type of reference solution
 */
template<typename REFSOL>
class ref05::GivenData
{
public:
    //! Template parameter: reference solution
    typedef REFSOL ReferenceSolution;
    
    //! @name Linear Algebra types
    //@{
    typedef typename ReferenceSolution::VecDoF    VecDoF;
    typedef typename ReferenceSolution::VecDim    VecDim;
    typedef typename ReferenceSolution::GradType  GradType;
    //@}

    //! Initialise reference to reference solution
    GivenData( const ReferenceSolution& rs )
        : referenceSolution_( rs ) { }

    //--------------------------------------------------------------------------
    /** Dirichlet constraint \f$ \bar{u} = u  \f$
     *  \param[in]     x      Coordinate
     *  \param[in,out] doFPtr Pointer to degree of freedom on boundary
     */
    template<typename DOF>
    void dirichletBC( const VecDim& x, DOF* doFPtr ) const
    {
        const double value = ( referenceSolution_.evaluate( x ) )[0];

        // make x_i = 0 a Dirichlet boundary
        bool onDirichletBdr = false;
        for ( unsigned d = 0; d < ReferenceSolution::dim; d++ ) {
            if ( std::abs( x[d] - 0. ) < coordTol ) onDirichletBdr = true;
        }

        if ( onDirichletBdr ) {
            if ( doFPtr -> isActive(0) )
                doFPtr -> constrainValue( 0, value );
        }
        return;
    }

    //--------------------------------------------------------------------------
    /** Force term \f$ f = - \nabla^2 u  \f$
     *  \param[in] x      Coordinate
     *  \return           The value of the force term
     */
    VecDoF forceFun( const VecDim& x ) const
    {
        const VecDoF result = - referenceSolution_.laplacian( x );
        return result;
    }

#if 0
    //--------------------------------------------------------------------------
    // For demonstration only, alternative way to prescribe body forces
    /** Force term \f$ f = - \nabla^2 u  \f$
     *  \param[in] gep    Pointer to geometry element
     *  \param[in] xi     Local evaluation coordinate
     *  \return           The value of the force term
     */
    template<typename GEOMELEMENT>
    VecDoF forceFun2( const GEOMELEMENT* gep, const VecDim& xi ) const
    {
        const VecDoF result =
            -referenceSolution_.laplacian( base::Geometry<GEOMELEMENT>()( gep, xi ) );
        return result;
    }
    
#endif
    //--------------------------------------------------------------------------
    /** Neumann boundary condition \f$ \bar{t} = \nabla u \cdot \vec{u} \f$
     *  \param[in] x      Coordinate
     *  \param[in] normal Normal vector
     *  \return           The value of the co-normal derivative 
     */
    VecDoF neumannBC( const VecDim& x, const VecDim& normal ) const
    {
        const GradType gradient = referenceSolution_.evaluateGradient( x );
        const VecDoF result = gradient.transpose() * normal;
        return result;
    }
    
private:
    //! Reference to reference solution :)
    const ReferenceSolution& referenceSolution_;
};


#endif
