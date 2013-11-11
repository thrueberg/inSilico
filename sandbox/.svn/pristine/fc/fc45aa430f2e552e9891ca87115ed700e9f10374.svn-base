#include <base/asmb/BodyForce.hpp>

//------------------------------------------------------------------------------
template<typename MESH, unsigned DEGREE, unsigned NHIST>
class Solute
{
public:
    //! @name Template parameter
    //@{
    typedef MESH Mesh;
    static const unsigned degree = DEGREE;
    static const unsigned nHist  = NHIST;
    //@}

    //!@name Solute Field
    //@{
    static const unsigned doFSize = 1;
    typedef base::fe::Basis<Mesh::Element::shape,degree> FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>           Field;
    typedef typename Field::DegreeOfFreedom              DoF;    
    //@}
    
    //! Constructor with mesh reference and material parameters
    Solute( Mesh& mesh, const double D, const double phi, const double k,
            const bool withConvection = true )
        : mesh_( mesh ), D_( D ), phi_( phi ), k_( k ),
          withConvection_( withConvection )
    {
        // empty
    }

    //! Generate degrees of freedom
    void generateDoFs( )
    {
        base::dof::generate<FEBasis>( mesh_, field_ );
    }

    //! Set initial concentration
    template<typename SETFUN>
    void setInitialConcentration( SETFUN setFun )
    {
        base::dof::setField( mesh_, field_, setFun );
    }
    
    //! Constrain the boundary
    template<typename MESHBOUNDARY, typename DIRI>
    void constrainBoundary( const MESHBOUNDARY& meshBoundary, DIRI diri )
    {
        base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                               meshBoundary.end(),
                                               mesh_, field_, diri );
    }

    //! Number dofs
    std::size_t numberDoFs()
    {
        numDoFs_ =
            base::dof::numberDoFsConsecutively( field_.doFsBegin(),
                                                field_.doFsEnd() );
        return numDoFs_;
    }

    //! Advance in time
    template<typename MSM,      typename DISPLACEMENT,
             typename PRESSURE, typename QUADRATURE,
             typename BODYFORCE>
    void advanceInTime( const QUADRATURE& quadrature,
                        DISPLACEMENT& displacement, PRESSURE& pressure,
                        BODYFORCE bodyForce, 
                        const double deltaT, const unsigned step )
    {
        // bind fields to solute
        typedef base::asmb::FieldBinder<Mesh,Field,DISPLACEMENT,PRESSURE> FieldBinder;
        typedef typename FieldBinder::template TupleBinder<1,1,2,3>::Type  TupleBinder;
        FieldBinder fieldBinder( mesh_, field_, displacement, pressure );

        // matrix kernels for the solute
        typedef heat::Laplace<typename TupleBinder::Tuple> Laplace;
        Laplace laplace( phi_ * D_ );

        // create a solver
        base::solver::Eigen3 solver( numDoFs_ );

        // compute inertia terms, d/dt, due to time integration
        base::time::computeInertiaTerms<TupleBinder,MSM>( quadrature, solver,
                                                          fieldBinder, deltaT, step,
                                                          phi_ );

        // convection term
        if ( withConvection_ ) {
            typedef Convection<typename TupleBinder::Tuple> Convection;
            Convection convection( phi_, k_, deltaT, step );
            base::asmb::stiffnessMatrixComputation<TupleBinder>( quadrature, solver,
                                                                 fieldBinder, convection );
        }
            
        // compute diffusion term
        base::asmb::stiffnessMatrixComputation<TupleBinder>( quadrature, solver,
                                                             fieldBinder, laplace );

        // compute body force term
        base::asmb::bodyForceComputation2<TupleBinder>( quadrature, solver, fieldBinder,
                                                        bodyForce );

        solver.finishAssembly();
        
        solver.superLUSolve();
        
        base::dof::setDoFsFromSolver( solver, field_ );
        
        std::for_each( field_.doFsBegin(), field_.doFsEnd(),
                       boost::bind( &DoF::pushHistory, _1 ) );
    }

    //! Accessors
    Field& soluteField() { return field_; }
    
private:
    Mesh&        mesh_;           //!< Access to mesh
    Field        field_;          //!< Field of solute

    const double D_;      //!< Diffusivity
    const double phi_;    //!< Fluid volume fraction
    const double k_;      //!< Hydraulic permeabilityx

    const bool withConvection_;
    
    //! Values of dof numbers
    std::size_t numDoFs_;
};

