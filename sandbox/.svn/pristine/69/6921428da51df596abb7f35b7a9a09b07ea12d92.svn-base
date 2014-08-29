//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ElasticityDriver.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef tfm_elasticitydriver_hpp
#define tfm_elasticitydriver_hpp

//------------------------------------------------------------------------------
#include <mat/hypel/StVenant.hpp>
#include <solid/CompressibleDriver.hpp>
#include "NeumannForce2.hpp"

//------------------------------------------------------------------------------
namespace tfm{

    template<typename BVP>
    class ElasticityDriver;
}

//------------------------------------------------------------------------------
template<typename BVP>
class tfm::ElasticityDriver
    : public solid::CompressibleDriver<BVP,mat::hypel::StVenant>
{
public:
    typedef BVP      BoundaryValueProblem;
    typedef mat::hypel::StVenant  Material;

    typedef solid::CompressibleDriver<BVP,Material> Basis;

    // other type of neumann BC
    typedef typename base::Vector<BoundaryValueProblem::dim-1>::Type SurfVecDim;
    typedef typename BoundaryValueProblem::BoundaryMesh::Element     SurfElement;
    typedef boost::function<typename BoundaryValueProblem::VecDoF(
        const SurfElement*, const SurfVecDim& )> SurfaceForceFun2;

    // other type of body forces
    typedef boost::function<typename BoundaryValueProblem::VecDoF(
        const typename BoundaryValueProblem::Mesh::Element*,
        const typename BoundaryValueProblem::VecDim& ) > BodyForceFun2;

    //! Call with reference to a material
    ElasticityDriver( BoundaryValueProblem& boundaryValueProblem,
                      Material& material )
        : Basis( boundaryValueProblem, material ),
          withNeumann2_( false ), withBodyForce2_( false )
    {
    }
    
public:

    //! Register a surface force function
    void neumann( SurfaceForceFun2 surfaceForceFun )
    {
        surfaceForceFun_ = surfaceForceFun;
        Basis::boundaryValueProblem_.generateBoundaryMesh();
        withNeumann2_ = true;
    }

    void force( BodyForceFun2 bodyForceFun )
    {
        bodyForceFun_ = bodyForceFun;
        withBodyForce2_ = true;
    }

    //--------------------------------------------------------------------------
    template<typename SMESH>
    void  solve( SMESH& surfaceMesh )
    {
        // linear elasticity -> start from zero
        base::dof::clearDoFs( Basis::boundaryValueProblem_.getField() );
        
        // number the degrees of freedom
        if ( Basis::numDoFs_ == 0 )
            Basis::numDoFs_ = Basis::boundaryValueProblem_.numberDoFs();

        // Quadrature
        static const unsigned kernelDegree = 2 * BoundaryValueProblem::fieldDeg;
        base::Quadrature<       kernelDegree,BoundaryValueProblem::shape> quadrature;
        base::SurfaceQuadrature<kernelDegree,BoundaryValueProblem::shape> surfaceQuadrature;

        // Solver
        base::solver::Eigen3 solver( Basis::numDoFs_ );

        Basis::assemble( solver, quadrature, surfaceQuadrature );

        // alternative surface force computation
        if ( withNeumann2_ ) {
            typename BoundaryValueProblem::SurfaceFieldBinder
                sfb( surfaceMesh,
                     Basis::boundaryValueProblem_.getField() );
        
            base::asmb::neumannForceComputation2<typename BoundaryValueProblem::SUU>(
                surfaceQuadrature, solver, sfb, surfaceForceFun_ );
        }

        // alternative body force computation
        if ( withBodyForce2_ ) {
            typename BoundaryValueProblem::FieldBinder
                fb( Basis::boundaryValueProblem_.getMesh(),
                    Basis::boundaryValueProblem_.getField() );
            
            base::asmb::bodyForceComputation2<typename BoundaryValueProblem::UU>(
                quadrature, solver, fb, bodyForceFun_ );
        }

        solver.finishAssembly();
        
        solver.choleskySolve();

        Basis::boundaryValueProblem_.addToDoFsFromSolver( solver );
    }
    
private:
    //! Surface force function
    SurfaceForceFun2 surfaceForceFun_;
    BodyForceFun2    bodyForceFun_;

    bool withNeumann2_;
    bool withBodyForce2_;
};

#endif
