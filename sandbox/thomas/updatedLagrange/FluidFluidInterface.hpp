#ifndef fluidfluidinterface_h
#define fluidfluidinterface_h

#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>
#include <base/nitsche/Parameters.hpp>

#include "SurfaceTension.hpp"
//------------------------------------------------------------------------------
template<typename SURFMESH, typename FLUID>
class FluidFluidInterface
{
public:
    //! Template parameter
    typedef SURFMESH SurfaceMesh;
    typedef FLUID    Fluid;

    //! Deduced types
    typedef typename Fluid::Velocity     Velocity;
    typedef typename Fluid::Pressure     Pressure;
    
    //! Surface field binder
    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,
                                           Velocity, Pressure,
                                           Velocity, Pressure> SurfaceFieldBinder;

    // U-U
    typedef typename SurfaceFieldBinder::template TupleBinder<1,1>::Type SU1U1;
    typedef typename SurfaceFieldBinder::template TupleBinder<1,3>::Type SU1U2;
    typedef typename SurfaceFieldBinder::template TupleBinder<3,1>::Type SU2U1;
    typedef typename SurfaceFieldBinder::template TupleBinder<3,3>::Type SU2U2;

    // U-P
    typedef typename SurfaceFieldBinder::template TupleBinder<1,2>::Type SU1P1;
    typedef typename SurfaceFieldBinder::template TupleBinder<1,4>::Type SU1P2;
    typedef typename SurfaceFieldBinder::template TupleBinder<3,2>::Type SU2P1;
    typedef typename SurfaceFieldBinder::template TupleBinder<3,4>::Type SU2P2;

    // P-U
    typedef typename SurfaceFieldBinder::template TupleBinder<2,1>::Type SP1U1;
    typedef typename SurfaceFieldBinder::template TupleBinder<2,3>::Type SP1U2;
    typedef typename SurfaceFieldBinder::template TupleBinder<4,1>::Type SP2U1;
    typedef typename SurfaceFieldBinder::template TupleBinder<4,3>::Type SP2U2;


    static const base::Shape shape = Velocity::Element::shape;
    typedef base::cut::Cell<shape> Cell;
    

    FluidFluidInterface( SurfaceMesh&  surfaceMesh,
                         Velocity&     velocity1,
                         Pressure&     pressure1, 
                         Velocity&     velocity2,
                         Pressure&     pressure2,
                         std::vector<Cell>& cells,
                         const double viscosity1,
                         const double viscosity2,
                         const double sigma )
        : surfaceMesh_(  surfaceMesh ),
          velocity1_(     velocity1 ),
          pressure1_(     pressure1 ),
        velocity2_(     velocity2 ),
        pressure2_(     pressure2 ),
        sfb_( surfaceMesh_, velocity1_, pressure1_, velocity2_, pressure2_ ),
        cells_( cells ),
        viscosity1_( viscosity1 ),
        viscosity2_( viscosity2 ),
        sigma_(      sigma ),
        stressDiv1_( viscosity1_ ),
        stressDiv2_( viscosity2_ ),
        gradP_()
    { }

    //--------------------------------------------------------------------------
    template<typename SOLVER>
    void registerInSolver( SOLVER& solver )
    {
        // U-U
        solver.template registerFields<SU1U1>( sfb_ );
        solver.template registerFields<SU1U2>( sfb_ );
        solver.template registerFields<SU2U1>( sfb_ );
        solver.template registerFields<SU2U2>( sfb_ );

        // U-P
        solver.template registerFields<SU1P1>( sfb_ );
        solver.template registerFields<SU1P2>( sfb_ );
        solver.template registerFields<SU2P1>( sfb_ );
        solver.template registerFields<SU2P2>( sfb_ );

        // P-U
        solver.template registerFields<SP1U1>( sfb_ );
        solver.template registerFields<SP1U2>( sfb_ );
        solver.template registerFields<SP2U1>( sfb_ );
        solver.template registerFields<SP2U2>( sfb_ );
    }

    //--------------------------------------------------------------------------
    template<typename SURFQUAD, typename SOLVER>
    void assemblePenaltyTerms( const SURFQUAD& surfaceQuadrature,
                               SOLVER& solver,
                               const double penaltyFac )
    {
        //base::nitsche::ImmersedInterface<Cell> ip( viscosity1_, viscosity2_, cells_ );
        base::nitsche::PrescribedParameters ip(  1.0, 0.5 );

        // U-U only
        base::nitsche::penaltyLHS<SU1U1>( surfaceQuadrature, solver, sfb_, ip,  penaltyFac );
        base::nitsche::penaltyLHS<SU2U2>( surfaceQuadrature, solver, sfb_, ip,  penaltyFac );
        base::nitsche::penaltyLHS<SU1U2>( surfaceQuadrature, solver, sfb_, ip, -penaltyFac );
        base::nitsche::penaltyLHS<SU2U1>( surfaceQuadrature, solver, sfb_, ip, -penaltyFac );

    }

    //--------------------------------------------------------------------------
    template<typename SURFQUAD, typename SOLVER>
    void assembleEnergyTerms( const SURFQUAD& surfaceQuadrature,
                              SOLVER& solver )
    {
        base::nitsche::PrescribedParameters ip(  1., 0.5 );
        //base::nitsche::ImmersedInterface<Cell> ip( viscosity1_, viscosity2_, cells_ );

        //----------------------------------------------------------------------
        // U1 -> U1
        base::nitsche::primalEnergyLHS<SU1U1>( stressDiv1_, surfaceQuadrature, solver,
                                               sfb_, ip, true, true );
        base::nitsche::dualEnergyLHS<SU1U1>(   stressDiv1_, surfaceQuadrature, solver,
                                               sfb_, ip, true, true );
        // U1 -> U2
        base::nitsche::primalEnergyLHS<SU2U1>( stressDiv1_, surfaceQuadrature, solver,
                                               sfb_, ip, true, false );
        base::nitsche::dualEnergyLHS<SU2U1>(   stressDiv1_, surfaceQuadrature, solver,
                                               sfb_, ip, true, false );

        // U2 -> U1
        base::nitsche::primalEnergyLHS<SU1U2>( stressDiv2_, surfaceQuadrature, solver,
                                               sfb_, ip, false, true );
        base::nitsche::dualEnergyLHS<SU1U2>(   stressDiv2_, surfaceQuadrature, solver,
                                               sfb_, ip, false, true );
        // U2 -> U2
        base::nitsche::primalEnergyLHS<SU2U2>( stressDiv2_, surfaceQuadrature, solver,
                                               sfb_, ip, false, false );
        base::nitsche::dualEnergyLHS<SU2U2>(   stressDiv2_, surfaceQuadrature, solver,
                                               sfb_, ip, false, false );

        //----------------------------------------------------------------------
        // P1 -> U1
        base::nitsche::primalEnergyLHS<SU1P1>( gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, true, true );
        base::nitsche::dualEnergyLHS<SU1P1>(   gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, true, true );
        // P1 -> U2
        base::nitsche::primalEnergyLHS<SU2P1>( gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, true, false );
        base::nitsche::dualEnergyLHS<SU2P1>(   gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, true, false );

        // P2 -> U1
        base::nitsche::primalEnergyLHS<SU1P2>( gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, false, true );
        base::nitsche::dualEnergyLHS<SU1P2>(   gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, false, true );
        // P2 -> U2
        base::nitsche::primalEnergyLHS<SU2P2>( gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, false, false );
        base::nitsche::dualEnergyLHS<SU2P2>(   gradP_, surfaceQuadrature, solver,
                                               sfb_, ip, false, false );

        
    }

    //--------------------------------------------------------------------------
    template<typename SURFQUAD, typename SOLVER>
    void surfaceTension( const SURFQUAD& surfaceQuadrature,
                         SOLVER& solver )

    {
        //base::nitsche::ImmersedInterface<Cell> ip( viscosity1_, viscosity2_, cells_ );
        base::nitsche::PrescribedParameters ip(  1., 0.5 );
        
        ::surfaceTension<SU1U2>( surfaceQuadrature, solver, sfb_, ip, sigma_, true  );
        ::surfaceTension<SU2U2>( surfaceQuadrature, solver, sfb_, ip, sigma_, false );

    }

    //--------------------------------------------------------------------------
    SurfaceFieldBinder& getBinder() { return sfb_; }

private:
    SurfaceMesh&  surfaceMesh_;
    Velocity&     velocity1_;
    Pressure&     pressure1_;
    Velocity&     velocity2_;
    Pressure&     pressure2_;

    SurfaceFieldBinder         sfb_;

    std::vector<Cell>& cells_;
    
    const double viscosity1_;
    const double viscosity2_;
    const double sigma_;
    
    const typename Fluid::StressDiv  stressDiv1_;
    const typename Fluid::StressDiv  stressDiv2_;
    const typename Fluid::GradP      gradP_;
};

#endif
