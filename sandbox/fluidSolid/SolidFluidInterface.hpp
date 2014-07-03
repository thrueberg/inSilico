#ifndef solidfluidinterface_h
#define solidfluidinterface_h

#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>
#include <base/nitsche/Parameters.hpp>

//------------------------------------------------------------------------------
template<typename SURFMESH, typename SOLID, typename FLUID>
class SolidFluidInterface
{
public:
    //! Template parameter
    typedef SURFMESH SurfaceMesh;
    typedef SOLID    Solid;
    typedef FLUID    Fluid;

    //! Deduced types
    typedef typename Solid::Displacement Displacement;
    typedef typename Fluid::Velocity     Velocity;
    typedef typename Fluid::Pressure     Pressure;
    
    //! Surface field binder
    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh, Displacement,
                                           Velocity, Pressure> SurfaceFieldBinder;
    
    typedef typename SurfaceFieldBinder::template TupleBinder<1,1>::Type SDD;
    typedef typename SurfaceFieldBinder::template TupleBinder<1,2>::Type SDU;
    typedef typename SurfaceFieldBinder::template TupleBinder<2,1>::Type SUD;
    typedef typename SurfaceFieldBinder::template TupleBinder<2,2>::Type SUU;

    typedef typename SurfaceFieldBinder::template TupleBinder<1,3>::Type SDP;
    typedef typename SurfaceFieldBinder::template TupleBinder<3,1>::Type SPD;
    typedef typename SurfaceFieldBinder::template TupleBinder<2,3>::Type SUP;

    //! Type-coupling for solid-fluid
    typedef typename base::asmb::DomainFieldElementPointerTuple<typename SDU::Tuple>::Type DU;
    typedef typename base::asmb::DomainFieldElementPointerTuple<typename SDP::Tuple>::Type DP;
    typedef fluid::StressDivergence<DU> StressDiv2;
    typedef fluid::PressureGradient<DP> GradP2;


    SolidFluidInterface( SurfaceMesh&  surfaceMesh,
                         Displacement& displacement,
                         Velocity&     velocity,
                         Pressure&     pressure, 
                         const double viscosity )
        : surfaceMesh_(  surfaceMesh ),
          displacement_( displacement ),
          velocity_(     velocity ),
          pressure_(     pressure ),
          sfb_( surfaceMesh_, displacement_, velocity_, pressure_ ),
          viscosity_(  viscosity ),
        stressDiv_(  viscosity_ ),
        gradP_(),
        stressDiv2_( viscosity_ ),
        gradP2_() { }

    //--------------------------------------------------------------------------
    template<typename SOLVER>
    void registerInSolver( SOLVER& solver )
    {
        solver.template registerFields<SUD>( sfb_ );
        solver.template registerFields<SDU>( sfb_ );
        solver.template registerFields<SDP>( sfb_ );
        solver.template registerFields<SPD>( sfb_ );
    }

    //--------------------------------------------------------------------------
    template<typename SURFQUAD, typename SOLVER>
    void assemblePenaltyTerms( const SURFQUAD& surfaceQuadrature,
                               SOLVER& solver,
                               const double penaltyFac,
                               const double dt )
    {
        base::nitsche::OuterBoundary p( penaltyFac );

        // Solid - fluid
        base::nitsche::penaltyLHS<SDD>( surfaceQuadrature, solver,
                                        sfb_, p,  1.0 / dt );
        base::nitsche::penaltyLHS<SDU>( surfaceQuadrature, solver,
                                        sfb_, p, -1.0 );
        base::nitsche::penaltyRHSInterface<SDD>( surfaceQuadrature, solver,
                                                 sfb_, p, 1.0/dt );

        // Fluid - solid
        base::nitsche::penaltyLHS<SUU>( surfaceQuadrature, solver,
                                        sfb_, p,  1.0 );
        base::nitsche::penaltyLHS<SUD>( surfaceQuadrature, solver,
                                        sfb_, p, -1.0 /dt );
        base::nitsche::penaltyRHSInterface<SUD>( surfaceQuadrature, solver,
                                                 sfb_, p, -1.0/dt );
    }

    //--------------------------------------------------------------------------
    template<typename SURFQUAD, typename SOLVER>
    void assembleEnergyTerms( const SURFQUAD& surfaceQuadrature,
                              SOLVER& solver,
                              const double dt )
    {
        // +sigma^f(u,p) n dd
        base::nitsche::PrescribedParameters bla1( 1., 1.0 );
        base::nitsche::primalEnergyLHS<SDU>( stressDiv2_,
                                             surfaceQuadrature, solver,
                                             sfb_, bla1 );
        base::nitsche::primalEnergyLHS<SDP>( gradP2_,
                                             surfaceQuadrature, solver,
                                             sfb_, bla1 );

        // -sigma^f(u,p) n du
        base::nitsche::PrescribedParameters bla2( 1., -1.0 );
        base::nitsche::primalEnergyLHS<SUU>( stressDiv_,
                                             surfaceQuadrature, solver,
                                             sfb_, bla2 );
        
        base::nitsche::primalEnergyLHS<SUP>( gradP_,
                                             surfaceQuadrature, solver,
                                             sfb_, bla2 );

        // -sigma^f(du,pm dp) n (Delta d/Delta t)
        base::nitsche::PrescribedParameters bla3(  1., -1.0/dt );
        base::nitsche::dualEnergyLHS<SDU>( stressDiv2_,
                                           surfaceQuadrature, solver,
                                           sfb_, bla3 );
        
        base::nitsche::dualEnergyLHS<SDP>( gradP2_,
                                           surfaceQuadrature, solver,
                                           sfb_, bla3 );
                
        // +sigma^f(du,pm dp) n u
        base::nitsche::PrescribedParameters bla4(  1., 1.0 );
        base::nitsche::dualEnergyLHS<SUU>( stressDiv_,
                                           surfaceQuadrature, solver,
                                           sfb_, bla4 );
        
        base::nitsche::dualEnergyLHS<SUP>( gradP_,
                                           surfaceQuadrature, solver,
                                           sfb_, bla4 );
                
        // Residual: -sigma^f(du,pm dp) n (Delta d/Delta t)
        base::nitsche::PrescribedParameters bla5( 1., 1.0 / dt );
        base::nitsche::energyResidualTransposed<SDU>( stressDiv2_,
                                                      surfaceQuadrature,
                                                      solver, sfb_, bla5 );
        
        base::nitsche::energyResidualTransposed<SDP>( gradP2_,
                                                      surfaceQuadrature,
                                                      solver, sfb_, bla5 );
    }

private:
    SurfaceMesh&  surfaceMesh_;
    Displacement& displacement_;
    Velocity&     velocity_;
    Pressure&     pressure_;

    SurfaceFieldBinder         sfb_;
    
    const double viscosity_;
    
    const typename Fluid::StressDiv  stressDiv_;
    const typename Fluid::GradP      gradP_;
    StressDiv2 stressDiv2_;
    GradP2     gradP2_;
};

#endif
