#ifndef fluid_h
#define fluid_h

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>

#include <base/Field.hpp>
#include <base/fe/Basis.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>

#include <fluid/Stokes.hpp>
#include <fluid/GalerkinLeastSquares.hpp>

#include <base/time/BDF.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

#include <fluid/Convection.hpp>
#include "Convection.hpp"

//------------------------------------------------------------------------------
template<typename MESH, bool STOKESSTABIL>
class Fluid
{
public:
    //! Template parameter
    typedef MESH Mesh;
    static const bool     stokesStabil = STOKESSTABIL;

    //! Shape function degrees
    static const unsigned degreeU      = (stokesStabil? 1 : 2 );
    static const unsigned degreeP      = 1;

    //! Multi-step method
    typedef base::time::BDF<1> MSM;

    //! Number of history terms
    static const unsigned numHist = MSM::numSteps;

    //! Deduced parameter
    static const unsigned    dim   = Mesh::Node::dim;
    static const base::Shape shape = Mesh::Element::shape;

    //! Field types
    typedef base::fe::Basis<shape,degreeU>         FEBasisU;
    typedef base::Field<FEBasisU,dim,numHist>      Velocity;
    typedef typename Velocity::DegreeOfFreedom     DoFU;
    
    typedef base::fe::Basis<shape,degreeP>         FEBasisP;
    typedef base::Field<FEBasisP,1>                Pressure;
    typedef typename Pressure::DegreeOfFreedom     DoFP;

    //! Field binder
    typedef base::asmb::FieldBinder<Mesh,Velocity,Pressure> FieldBinder;
    typedef typename FieldBinder::template TupleBinder<1,1>::Type   UU;
    typedef typename FieldBinder::template TupleBinder<1,2>::Type   UP;
    typedef typename FieldBinder::template TupleBinder<2,1>::Type   PU;
    typedef typename FieldBinder::template TupleBinder<2,2>::Type   PP;
    typedef typename FieldBinder::template TupleBinder<1,1,1>::Type UUU;

    //! Kernels
    typedef fluid::StressDivergence<  typename UU::Tuple>  StressDiv;
    typedef fluid::PressureGradient<  typename UP::Tuple>  GradP;
    typedef fluid::VelocityDivergence<typename PU::Tuple>  DivU;
    typedef ::Convection<typename UU::Tuple>               Convection1;
    typedef fluid::Convection<typename UUU::Tuple>         Convection2;

    // Stabilisation kernels
    typedef fluid::gls::StressDivergence<   typename UU::Tuple>  StressDivStabil;
    typedef fluid::gls::PressureGradient2<  typename UP::Tuple>  GradPStabil2;
    typedef fluid::gls::VelocityDivergence2<typename PU::Tuple>  DivUStabil2;
    typedef fluid::gls::PressureLaplace<    typename PP::Tuple>  PressureLaplace;


    //--------------------------------------------------------------------------
    Fluid( Mesh& mesh,
           const double viscosity, 
           const double alpha,
           const bool   dynamic = false,
           const double density = 1.0,
           const double stepSize = 1.0 )
        : mesh_(     mesh ),
          viscosity_( viscosity ), 
          alpha_( alpha ), dw_( true ),
          dynamic_( dynamic ),
          density_( density ),
          stepSize_( stepSize ),
          stressDiv_(  viscosity_ ),
          gradP_(),
          divU_(),
          convection1_( density_ ),
          convection2_( density_ ),
          stressDivStabil_( alpha_, viscosity_, dw_ ),
          gradPStabil2_(    alpha_, viscosity_, dw_ ),
          divUStabil2_(     alpha_, viscosity_ ),
          pressureLaplace_( alpha_ )
    {
        base::dof::generate<FEBasisU>( mesh_, velocity_ );
        base::dof::generate<FEBasisP>( mesh_, pressure_ );
    }

    //--------------------------------------------------------------------------
    template<typename SOLVER>
    void registerInSolver( SOLVER& solver )
    {
        FieldBinder fb( mesh_, velocity_, pressure_ );
        solver.template registerFields<UU>( fb );
        solver.template registerFields<UP>( fb );
        solver.template registerFields<PU>( fb );
        if ( stokesStabil ) solver.template registerFields<PP>( fb );
    }


    //--------------------------------------------------------------------------
    template<typename QUAD, typename SOLVER>
    void assembleBulk( const QUAD& quadrature,
                       SOLVER& solver, const unsigned iter = 0 )
    {
        FieldBinder fb( mesh_, velocity_, pressure_ );
        
        //
        const bool bla = false;
        
        // Stiffness matrices of the fluid
        base::asmb::stiffnessMatrixComputation<UU>( quadrature, solver,
                                                    fb, stressDiv_, bla );
        base::asmb::stiffnessMatrixComputation<UP>( quadrature, solver,
                                                    fb, gradP_, bla );
        base::asmb::stiffnessMatrixComputation<PU>( quadrature, solver,
                                                    fb, divU_, bla );

        // Stabilisation terms of the fluid
        if ( stokesStabil ){
            base::asmb::stiffnessMatrixComputation<UU>( quadrature, solver,
                                                        fb, stressDivStabil_ );
            base::asmb::stiffnessMatrixComputation<UP>( quadrature, solver,
                                                        fb, gradPStabil2_ );
            base::asmb::stiffnessMatrixComputation<PU>( quadrature, solver,
                                                        fb, divUStabil2_ );
            base::asmb::stiffnessMatrixComputation<PP>( quadrature, solver,
                                                        fb, pressureLaplace_ );
        }

        if ( dynamic_ ) {

            // compute inertia terms, d/dt, due to time integration
            base::time::computeInertiaTerms<UU,MSM>( quadrature, solver,
                                                     fb, stepSize_, 1,
                                                     density_, false );

            // semi-explicit convection
            if ( iter == 0 ) {
                base::asmb::stiffnessMatrixComputation<UU>( quadrature, solver,
                                                            fb, convection1_ );
            }
            else {
                base::asmb::stiffnessMatrixComputation<UUU>( quadrature, solver,
                                                             fb, convection2_ );

            }
        }

    }

    //--------------------------------------------------------------------------
    template<typename QUAD, typename SOLVER, typename FFUN>
    void bodyForce( const QUAD& quadrature,
                    SOLVER& solver, FFUN forceFun )
    {
        FieldBinder fb( mesh_, velocity_, pressure_ );
        
        base::asmb::bodyForceComputation<UU>( quadrature, solver, fb, forceFun );

        if ( stokesStabil ) {
            fluid::gls::bodyForceComputation<UU,PU>( quadrature, solver,
                                                     fb, forceFun,
                                                     alpha_, viscosity_, dw_ );
        }

    }

    //--------------------------------------------------------------------------
    Mesh&     getMesh()     { return mesh_; }
    Velocity& getVelocity() { return velocity_; }
    Pressure& getPressure() { return pressure_; }

private:
    // Fields
    Mesh&     mesh_;
    Velocity  velocity_;
    Pressure  pressure_;

    // Parameters
    const double viscosity_;
    
    const double alpha_;
    const bool   dw_;

    const bool   dynamic_;
    const double density_;
    const double stepSize_;
    
    // Fluid kernel
    const StressDiv   stressDiv_;
    const GradP       gradP_;
    const DivU        divU_;
    const Convection1 convection1_;
    const Convection2 convection2_;

    // Fluid stabilisation terms
    const StressDivStabil  stressDivStabil_;
    const GradPStabil2     gradPStabil2_;
    const DivUStabil2      divUStabil2_;
    const PressureLaplace  pressureLaplace_;
};


#endif
