//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   heat/PoissonDriver.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef heat_poissondriver_hpp
#define heat_poissondriver_hpp

//------------------------------------------------------------------------------
#include <fstream>
#include <boost/function.hpp>

#include <base/verify.hpp>
#include <base/Quadrature.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/auxi/functions.hpp>

#include <mat/thermal/IsotropicConstant.hpp>
#include <heat/Static.hpp>

//------------------------------------------------------------------------------
namespace heat{

    //! \ingroup driver
    template<typename BVP,
             typename MATERIAL = mat::thermal::IsotropicConstant>
    class PoissonDriver;
}

//------------------------------------------------------------------------------
/** Driver for the solution of a Poisson problem.
 *  The problem reads
 *  \f[
 *       - \nabla \cdot (D(u) \nabla u) = f \quad x \in \Omega
 *  \f]
 *  with Dirichlet boundary conditions
 *  \f[
 *        u = \bar{u} \quad x \in \Gamma_D
 *  \f]
 *  and Neumann boundary conditions
 *  \f[
 *       t(u)  = \bar{t} x \in \Gamma_N
 *  \f]
 *
 * 
 *  \tparam BVP      Generic boundary value problem
 *  \tparam MATERIAL Material type
 */
template<typename BVP, typename MATERIAL>
class heat::PoissonDriver
{
public:
    //! @name Template parameter
    //@{
    typedef BVP      BoundaryValueProblem;
    typedef MATERIAL Material;
    //@}

    //! Make sure that there are only scalar DoFs
    STATIC_ASSERT_MSG( (BoundaryValueProblem::doFSize==1),
                       "Only scalar DoFs are supported" );

    //! @name Function types
    //@{
    typedef typename BoundaryValueProblem::DirichletFun    DirichletFun;
    typedef typename BoundaryValueProblem::BodyForceFun    BodyForceFun;
    typedef typename BoundaryValueProblem::SurfaceForceFun SurfaceForceFun;
    //@}

    //! A function that always gives zero is used for default values
    typedef base::auxi::ReturnZeroVector<1> ZeroFun;

    //! The integral kernel for the bilinear form is this
    typedef heat::Static<Material,typename BoundaryValueProblem::UU::Tuple> Kernel;

    //! @name Constructors
    //@{

    //! Call with reference to a material
    PoissonDriver( BoundaryValueProblem& boundaryValueProblem,
                   const Material&       material )
        : boundaryValueProblem_( boundaryValueProblem ),
          material_(             material ),
          kernel_(               material_ )
    {
        this -> initialise_();
    }

    //! Call with material parameter
    PoissonDriver( BoundaryValueProblem& boundaryValueProblem,
                   const double materialParam = 1)
        : boundaryValueProblem_( boundaryValueProblem ),
          material_(             materialParam ),
          kernel_(               material_ )
    {
        this -> initialise_();
    }
    //@}

private:
    void initialise_()
    {
        bodyForceFun_    = ZeroFun();
        surfaceForceFun_ = ZeroFun();
        withNeumann_     = false;
        withBodyForce_   = false;
    }

public:

    //! Apply Dirichlet constraints
    void dirichlet( DirichletFun dirichletFun )
    {
        boundaryValueProblem_.applyDirichletConstraints( dirichletFun );
    }

    //! Register a surface force function
    void neumann( SurfaceForceFun surfaceForceFun )
    {
        surfaceForceFun_ = surfaceForceFun;
        boundaryValueProblem_.generateBoundaryMesh();
        withNeumann_ = true;
    }

    //! Register a body force function
    void force( BodyForceFun bodyForceFun )
    {
        bodyForceFun_  = bodyForceFun;
        withBodyForce_ = true;
    }

    //--------------------------------------------------------------------------
    /** Assemble the weak form of the Poisson problem.
     *  The components are
     *  - the bilinear form of the Laplace operator (or similar)
     *  - the body force domain term
     *  - the Neumann bc surface term
     *
     *  The last two terms are optional and only active if neumann() and/or
     *  force() have been called prior to this function.
     *
     *  \tparam SOLVER System solver type
     *  \tparam QUAD   Domain quadrature type
     *  \tparam SQUAD  Surface quadrature type
     *  \param[in,out] solver            Target of assembly routines
     *  \param[in]     quadrature        Numerical integration on the mesh
     *  \param[in]     surfaceQuadrature Numerical integration on the surface
     */
    template<typename SOLVER, typename QUAD, typename SQUAD>
    void assemble( SOLVER& solver,
                   const QUAD& quadrature,
                   const SQUAD& surfaceQuadrature )
    {
        // assemble the bilinear form of the laplace operator
        boundaryValueProblem_.assembleBilinearForm( kernel_, quadrature, solver );

        // assemble the forces from the body force term
        if ( withBodyForce_ )
            boundaryValueProblem_.applyBodyForce( quadrature, solver,
                                                  bodyForceFun_ );

        // assemble the forces from a boundary force term
        if ( withNeumann_ ) {
            boundaryValueProblem_.applyBoundaryForces( surfaceQuadrature,
                                                       solver, surfaceForceFun_ );
        }
    }

    //--------------------------------------------------------------------------
    /** Black-box solve of a Poisson problem.
     *  Perform the following steps
     *  -  numbering of degrees of freedom
     *  -  create quadratures by introspection of the polynomial degree
     *  -  create a solver object
     *  -  call assembly
     *  -  solve
     *  -  set field from solution
     *
     *  \param[in] symmetric Assume a symmetric system matrix
     */
    void solve( const bool symmetric = true )
    {
        // number the degrees of freedom
        const std::size_t numDoFs = boundaryValueProblem_.numberDoFs();
        
        // Quadrature
        static const unsigned kernelDegree = 2 * BoundaryValueProblem::fieldDeg;
        base::Quadrature<       kernelDegree,BoundaryValueProblem::shape> quadrature;
        base::SurfaceQuadrature<kernelDegree,BoundaryValueProblem::shape> surfaceQuadrature;

        // Solver
        base::solver::Eigen3 solver( numDoFs );
        boundaryValueProblem_.registerInSolver( solver );

        // assemble system
        this -> assemble( solver, quadrature, surfaceQuadrature );

        // Solve
        solver.finishAssembly();

        if ( symmetric ) solver.choleskySolve();
        else             solver.luSolve();

        // Pass back to field
        boundaryValueProblem_.setDoFsFromSolver( solver );
    }

    Kernel& getKernel() { return kernel_; }
    
protected:
    //! Handler of the scalar equation
    BoundaryValueProblem& boundaryValueProblem_;

private:
    //! Body force function
    BodyForceFun    bodyForceFun_;
    //! Surface force function
    SurfaceForceFun surfaceForceFun_;

    //! Object of the material behaviour
    Material material_;
    //! Bilinear form of the 'Laplacian'
    Kernel   kernel_;

    //--------------------------------------------------------------------------
    //! @name Computational flags
    //@{
    
    //! Check if a Neumann BC is to be applied
    bool            withNeumann_;

    //! Check if a body force is to be applied
    bool            withBodyForce_;
    //@}
};

#endif
