//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   solid/CompressibleDriver.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef solid_compressibledriver_hpp
#define solid_compressibledriver_hpp

//------------------------------------------------------------------------------
#include <fstream>
#include <boost/function.hpp>

#include <base/verify.hpp>
#include <base/Quadrature.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/auxi/functions.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Distribute.hpp>

#include <solid/HyperElastic.hpp>

//------------------------------------------------------------------------------
namespace solid{

    //! \ingroup driver
    template<typename BVP, typename MATERIAL>
    class CompressibleDriver;
}

//------------------------------------------------------------------------------
/** Driver for the solution of a compressible hyperelasticity.
 * 
 *  \tparam BVP      Generic boundary value problem
 *  \tparam MATERIAL Material type
 */
template<typename BVP, typename MATERIAL>
class solid::CompressibleDriver
{
public:
    //! @name Template parameter
    //@{
    typedef BVP      BoundaryValueProblem;
    typedef MATERIAL Material;
    //@}

    //! Make sure that there are only dim-vector DoFs
    STATIC_ASSERT_MSG( (BoundaryValueProblem::doFSize==
                        BoundaryValueProblem::dim),
                       "Only scalar DoFs are supported" );

    //! @name Function types
    //@{
    typedef typename BoundaryValueProblem::DirichletFun    DirichletFun;
    typedef typename BoundaryValueProblem::BodyForceFun    BodyForceFun;
    typedef typename BoundaryValueProblem::SurfaceForceFun SurfaceForceFun;
    //@}

    //! A function that always gives zero is used for default values
    typedef base::auxi::ReturnZeroVector<BoundaryValueProblem::dim> ZeroFun;

    //! The integral kernel for the bilinear form is this
    typedef solid::HyperElastic<Material,typename BoundaryValueProblem::UU::Tuple> Kernel;

    //! Call with reference to a material
    CompressibleDriver( BoundaryValueProblem& boundaryValueProblem,
                        const Material&       material )
        : boundaryValueProblem_( boundaryValueProblem ),
          material_(             material ),
          kernel_(               material_ )
    {
        this -> initialise_();
    }


private:
    void initialise_()
    {
        numDoFs_         = 0;
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
        boundaryValueProblem_.assembleBilinearForm( kernel_, quadrature, solver,
                                                    true );

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
    /**
     *  \param[in] symmetric Assume a symmetric system matrix
     */
    unsigned  solve( const unsigned stepNum,
                     const unsigned maxIter,
                     const double tolerance )
    {
        // number the degrees of freedom
        if ( numDoFs_ == 0 )
            numDoFs_ = boundaryValueProblem_.numberDoFs();

        // rescale constraints in every load step: (newValue / oldValue)
        const double  constraintFactor =
            (stepNum == 0 ?
             static_cast<double>( stepNum+1 ) :
             static_cast<double>( stepNum+1 )/ static_cast<double>(stepNum) );

        // scale constraints
        base::dof::scaleConstraints( boundaryValueProblem_.getField(),
                                     constraintFactor );
        
        // Quadrature
        static const unsigned kernelDegree = 2 * BoundaryValueProblem::fieldDeg;
        base::Quadrature<       kernelDegree,BoundaryValueProblem::shape> quadrature;
        base::SurfaceQuadrature<kernelDegree,BoundaryValueProblem::shape> surfaceQuadrature;

        unsigned iter = 0;
        while ( iter < maxIter ) {
        
            // Solver
            base::solver::Eigen3 solver( numDoFs_ );
            boundaryValueProblem_.registerInSolver( solver );

            // assemble system
            this -> assemble( solver, quadrature, surfaceQuadrature );

            // Solve
            solver.finishAssembly();

            // norm of residual 
            const double normR = solver.norm();

            // convergence via residual norm
            if ( normR < tolerance ) { 
                break;
            }
            
            solver.choleskySolve();

            // Pass back to field
            boundaryValueProblem_.addToDoFsFromSolver( solver );

            // norm of displacement increment
            const double normDU = solver.norm();
            iter++;
            
            // convergence via increment
            if ( normDU < tolerance ) break;
        }

        return iter;
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

protected:
    std::size_t numDoFs_;

private:
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
