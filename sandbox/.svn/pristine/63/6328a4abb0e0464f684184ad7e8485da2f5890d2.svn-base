//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   tfm/ImmersedElasticityDriver.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef tfm_immersedelasticitydriver_hpp
#define tfm_immersedelasticitydriver_hpp

//------------------------------------------------------------------------------
#include <fstream>
#include <boost/function.hpp>

#include <base/verify.hpp>
#include <base/Quadrature.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/auxi/functions.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Distribute.hpp>

#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>
#include <base/nitsche/Parameters.hpp>


#include <solid/HyperElastic.hpp>

//------------------------------------------------------------------------------
namespace tfm{

    //! \ingroup thomas
    template<typename BVP, typename MATERIAL>
    class ImmersedElasticityDriver;
}

//------------------------------------------------------------------------------
template<typename BVP, typename MATERIAL>
class tfm::ImmersedElasticityDriver
{
public:
    typedef BVP      BoundaryValueProblem;
    typedef MATERIAL Material;

    static const unsigned    dim   = BoundaryValueProblem::dim;
    static const base::Shape shape = BoundaryValueProblem::shape;

    typedef typename BoundaryValueProblem::BoundaryMesh SurfaceMesh;
    
    // Functions
    typedef typename BoundaryValueProblem::DirichletFun  DirichletFun;

    typedef typename base::Vector<BoundaryValueProblem::dim-1>::Type SurfVecDim;
    typedef typename SurfaceMesh::Element                            SurfElement;
    typedef boost::function<typename BoundaryValueProblem::VecDoF(
        const SurfElement*, const SurfVecDim& )> SurfaceForceFun;


    typedef boost::function<typename BoundaryValueProblem::VecDoF(
        const typename BoundaryValueProblem::Mesh::Element*,
        const typename BoundaryValueProblem::VecDim& ) > BodyForceFun;

    typedef base::cut::LevelSet<dim> LevelSet;
    typedef base::cut::Cell<shape>   Cell;


    //! A function that always gives zero is used for default values
    typedef base::auxi::ReturnZeroVector<BoundaryValueProblem::dim> ZeroFun;

    typedef solid::HyperElastic<Material,typename BoundaryValueProblem::UU::Tuple> Kernel;

    static const unsigned kernelDegEstimate = 5;
    typedef base::cut::Quadrature<kernelDegEstimate,  shape> CutQuadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;



    ImmersedElasticityDriver( BoundaryValueProblem& bvp1,
                              BoundaryValueProblem& bvp2,
                              Material&             mat1,
                              Material&             mat2 )
        : bvp1_( bvp1 ), bvp2_( bvp2 ),
          material1_( mat1 ), material2_( mat2 ),
          elastic1_(  material1_ ),
          elastic2_(  material2_ ),
          withImmersedTraction_( false ), withBodyForce_( false )
    {
        numDoFs1_ = 0;
        numDoFs2_ = 0;

        immersedTraction_ = ZeroFun();
        bodyForceFun_     = ZeroFun();
        
    }

    template<typename SURFMESH>
    void immerse( const SURFMESH& surfaceMesh )
    {
        base::cut::bruteForce( bvp1_.getMesh(), surfaceMesh, true, levelSet_ );
        base::cut::generateCutCells( bvp1_.getMesh(), levelSet_, cells_ );

        // quadrature for in- and outside computations
        CutQuadrature cutQuadrature1( cells_, true  );
        CutQuadrature cutQuadrature2( cells_, false );

        // find geometry association for the dofs
        std::vector<std::pair<std::size_t,BoundaryValueProblem::VecDim> > doFLocation;
        base::dof::associateLocation( bvp1_.getField(), doFLocation );

        // compute supports
        std::vector<double> supports1, supports2;
        base::cut::supportComputation( mesh, bvp1_.getField(), cutQuadrature1, supports1 );
        base::cut::supportComputation( mesh, bvp2_.getField(), cutQuadrature2, supports2 );

        base::cut::stabiliseBasis( mesh, bvp1_.getField(), supports1, doFLocation );
        base::cut::stabiliseBasis( mesh, bvp2_.getField(), supports2, doFLocation );

        numDoFs1_ = bvp1_.numberDoFs( 0 );
        numDoFs2_ = bvp2_.numberDoFs( numDoFs1 );
    }


    void dirichlet( DirichletFun dirichletFun )
    {
        bvp1_.applyDirichletConstraints( dirichletFun );
        bvp2_.applyDirichletConstraints( dirichletFun );
    }

    void immersedTraction( SurfaceForceFun surfaceForceFun )
    {
        //!
    }

    void adjointBodyForce( BodyForceFun bodyForceFun )
    {
        bodyForceFun_ = bodyForceFun;
        withBodyForce_ = true;
    }


    void solve(  )
    {
        // force linear elasticity by starting from zero
        base::clearDoFs( bvp1_.getField() );
        base::clearDoFs( bvp2_.getField() );

        // quadrature for in- and outside computations
        CutQuadrature cutQuadrature1( cells_, true  );
        CutQuadrature cutQuadrature2( cells_, false );

        base::solver::Eigen3 solver( numDoFs1_ + numDoFs2_ );

        // assemble bulk
        bvp1_.assembleBilinearForm( elastic1_, cutQuadrature1, solver );
        bvp2_.assembleBilinearForm( elastic2_, cutQuadrature2, solver );

        if ( withBodyForce_ ) {

            // Inside the cell, there is no measurement
            // typename BoundaryValueProblem::FieldBinder
            //     fb1( bvp1_.getMesh(), bvp1_.getField() );
            // base::asmb::bodyForceComputation2<typename BoundaryValueProblem::UU>(
            //     cutQuadrature1, solver, fb1, adjointBodyForce_ );

            typename BoundaryValueProblem::FieldBinder
                fb2( bvp2_.getMesh(), bvp2_.getField() );
            base::asmb::bodyForceComputation2<typename BoundaryValueProblem::UU>(
                cutQuadrature2, solver, fb2, adjointBodyForce_ );
        }

        typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,
                                               BoundaryValueProblem::Field,
                                               BoundaryValueProblem::Field> SFB;
        typedef typename SurfaceFieldBinder::template TupleBinder<1,1>::Type SU1U1;
        typedef typename SurfaceFieldBinder::template TupleBinder<1,2>::Type SU1U2;
        typedef typename SurfaceFieldBinder::template TupleBinder<2,1>::Type SU2U1;
        typedef typename SurfaceFieldBinder::template TupleBinder<2,2>::Type SU2U2;

        SFB sfb( , bvp1_.getField(), bvp2_.getField() );

        {
            base::nitsche::PrescribedParameters ip(  1.0, 0.5 );

            // U-U only
            base::nitsche::penaltyLHS<SU1U1>( surfaceQuadrature, solver, sfb_, ip,  penaltyFac );
            base::nitsche::penaltyLHS<SU2U2>( surfaceQuadrature, solver, sfb_, ip,  penaltyFac );
            base::nitsche::penaltyLHS<SU1U2>( surfaceQuadrature, solver, sfb_, ip, -penaltyFac );
            base::nitsche::penaltyLHS<SU2U1>( surfaceQuadrature, solver, sfb_, ip, -penaltyFac );
        }


        solver.finishAssembly();
        solver.choleskySolve();
        bvp1_.setDoFsFromSolver( solver );
        bvp2_.setDoFsFromSolver( solver );

    }



private:
    BoundaryValueProblem& bvp1_;
    BoundaryValueProblem& bvp2_;

    Material material1_;
    Material material2_;

    SurfaceForceFun immersedTraction_;
    BodyForceFun    adjointBodyForce_;


    std::size_t numDoFs1_, numDoFs2_;

    Kernel elastic1_, elastic2_;
    
private:
    bool withImmersedTraction_;
    bool withBodyForce_;

private:
    std::vector<LevelSet> levelSet_;
    std::vector<Cell>     cells_;
};

#endif
