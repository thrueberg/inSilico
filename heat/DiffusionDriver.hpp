//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   heat/DiffusionDriver.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef heat_diffusiondriver_hpp
#define heat_diffusiondriver_hpp

//------------------------------------------------------------------------------
#include <fstream>
#include <boost/function.hpp>

#include <base/verify.hpp>
#include <base/Quadrature.hpp>
#include <base/solver/Eigen3.hpp>

#include <heat/PoissonDriver.hpp>

//------------------------------------------------------------------------------
namespace heat{

    //! \ingroup driver
    template<typename BVP,
             typename MATERIAL = mat::thermal::IsotropicConstant>
    class DiffusionDriver;
}


//------------------------------------------------------------------------------
/** Solution of diffusion problem.
 *  This object inherts from heat::PoissonDriver and adds time-dependent
 *  functionality:
 *  - setting initial conditions
 *  - performing one time integration step
 *
 *  \tparam BVP       Underlying boundary value problem
 *  \tparam MATERIAL  Type of material
 */
template<typename BVP, typename MATERIAL>
class heat::DiffusionDriver
    : public heat::PoissonDriver<BVP>
{
public:
    //! @name Template parameter
    //@{
    typedef BVP      BoundaryValueProblem;
    typedef MATERIAL Material;
    //@}

    //! Basis class: Poisson problem driver
    typedef heat::PoissonDriver<BoundaryValueProblem> Poisson;

    //! @name Constructors call PoissonDriver
    //@{
    DiffusionDriver( BoundaryValueProblem& boundaryValueProblem,
                     const Material& material )
        : Poisson(   boundaryValueProblem, material ),
          numDoFs_(  0 )
    { }

    DiffusionDriver( BoundaryValueProblem& boundaryValueProblem,
                     const double materialParam = 1 )
        : Poisson(   boundaryValueProblem, materialParam ),
          numDoFs_(  0 )
    { }
    //@}

    //! Set initial conditions given by a function
    void initialCondition( typename BoundaryValueProblem::DirichletFun initialCond )
    {
        Poisson::boundaryValueProblem_.setInitialValue( initialCond );
    }
    
    /** Advance one step in time.
     *  The sub-steps are
     *  - numbering of degrees of freedom (if still un-numbered)
     *  - create quadratures
     *  - create solver
     *  - call heat::PoissonDriver<BoundaryValueProblem>::assemble()
     *  - add the time-dependent terms
     *  - solve and set field to new values
     *  
     */
    void advanceInTime( const unsigned numStep,
                        const double   stepSize )
    {
        // number DoFs
        if ( numDoFs_ == 0 )
            numDoFs_ = Poisson::boundaryValueProblem_.numberDoFs();
        
        // Quadrature
        static const unsigned kernelDegree = 2 * BoundaryValueProblem::fieldDeg;
        base::Quadrature<       kernelDegree,BoundaryValueProblem::shape> quadrature;
        base::SurfaceQuadrature<kernelDegree,BoundaryValueProblem::shape> surfaceQuadrature;

        // Solver
        base::solver::Eigen3 solver( numDoFs_ );
        Poisson::boundaryValueProblem_.registerInSolver( solver );

        // assemble weak form
        Poisson::assemble( solver, quadrature, surfaceQuadrature );

        // time integration
        Poisson::boundaryValueProblem_.applyTimeIntegrator( numStep, stepSize,
                                                            Poisson::getKernel(),
                                                            quadrature, solver );

        // Solve
        solver.finishAssembly();
        solver.choleskySolve();

        // Pass back to field
        Poisson::boundaryValueProblem_.setDoFsFromSolver( solver );
    }

private:
    std::size_t numDoFs_; //!< number of DoFs in the system
};

#endif
