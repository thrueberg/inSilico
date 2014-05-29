template<typename SURFFIELDTUPLE, unsigned HIST=0>
class Penalty;


//------------------------------------------------------------------------------

template<unsigned HIST,
         typename SURFACETUPLEBINDER,
         typename SURFACEQUADRATURE,
         typename SOLVER,
         typename BOUNDFIELD,
         typename PARAMETER>
void computePenaltyResidual( const SURFACEQUADRATURE& surfaceQuadrature,
                             SOLVER& solver, BOUNDFIELD& boundField,
                             const PARAMETER& parameter,
                             const double multiplier )
{
    // object to compute the LHS penalty term
    typedef Penalty<typename SURFACETUPLEBINDER::Tuple, HIST> Penalty;
    Penalty penalty( multiplier );

    // integrator and assembler object
    typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                        typename SURFACETUPLEBINDER::Tuple>
        SurfaceForceInt;
                
    typename SurfaceForceInt::ForceKernel surfaceForceKernel =
        boost::bind( &Penalty::residual, &penalty, _1, _2, _3, _4 );
    SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                     surfaceQuadrature, solver );
                
    // apply
    typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
    typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
    for ( ; iter != end; ++iter ) {

        const double factor = multiplier * parameter.penaltyWeight( iter );
                    
        penalty.setFactor( factor );
                    
        surfaceForceInt( SURFACETUPLEBINDER::makeTuple( *iter ) );
    }
            
    return;

}

//------------------------------------------------------------------------------

template<typename SURFFIELDTUPLE, unsigned HIST>
class Penalty
{
public:
    //! Template parameter
    typedef SURFFIELDTUPLE SurfFieldTuple;

    //! @name Extract element types
    //@{
    typedef typename SurfFieldTuple::GeomElement  SurfaceElement;
    typedef typename SurfFieldTuple::TestElement  TestElement;
    typedef typename SurfFieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim GlobalVecDim;

    //! Type of domain element
    typedef typename SurfaceElement::DomainElement    DomainElement;

    //! Size of a DoF
    static const unsigned doFSize = TestElement::DegreeOfFreedom::size;

    //! Type of DoF value vector
    typedef typename base::Vector<doFSize,base::number>::Type VecDof;

    //--------------------------------------------------------------------------
    Penalty( const double factor ) : factor_( factor ) { }

    void setFactor( const double factor ) { factor_ = factor; }

    //--------------------------------------------------------------------------
    void residual( const SurfFieldTuple& surfFieldTuple,
                   const LocalVecDim&    eta,
                   const double          weight,
                   base::VectorD&        result ) const
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();
        const TestElement*    testEp  = surfFieldTuple.testElementPtr();
        const TrialElement*   trialEp = surfFieldTuple.trialElementPtr();
        
        // Get pointer to domain element
        const DomainElement* domainEp = surfEp -> getDomainElementPointer();

        // compute mesh size
        const double h = base::mesh::Size<DomainElement>::apply( domainEp );

        // Get surface metric
        GlobalVecDim dummy;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, dummy );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );
                
        // Evaluate the shape function
        typename TestElement::FEFun::FunArray testFunValues;
        (testEp -> fEFun()).evaluate( domainEp, xi, testFunValues );

        // deduce the size of every contribution
        const unsigned numRowBlocks = static_cast<unsigned>( testFunValues.size() );

        // Evaluate boundary coondition
        const VecDof  u = base::post::evaluateFieldHistory<HIST>( domainEp, trialEp, xi );

        // scalar multiplier of the whole entry
        const double aux = weight * detG * (factor_ / h);

        // Loop over shape functions
        for ( unsigned i = 0; i < numRowBlocks; i++ ) {

            for ( unsigned d = 0; d < doFSize; d++ ) {

                result[ i*doFSize + d ] += testFunValues[i] * u[d] * aux;
            }
        }

        return;
    }


private:
    double factor_;
};

//------------------------------------------------------------------------------
template<typename KERNEL, typename SURFFIELDTUPLE>
class Energy;

template<typename SURFACETUPLEBINDER,
         typename KERNEL,     typename SURFACEQUADRATURE,
         typename SOLVER,     typename BOUNDFIELD,
         typename PARAMETER>
void energyRHS2( const KERNEL&            kernel,
                 const SURFACEQUADRATURE& surfaceQuadrature,
                 SOLVER&                  solver, 
                 const BOUNDFIELD&        boundField,
                 const PARAMETER&         parameter,
                 const bool               inOut = true,
                 const bool               plusMinus = true )
{
    // object to compute the LHS penalty term
    typedef ::Energy<KERNEL,
                     typename SURFACETUPLEBINDER::Tuple> Energy;
    Energy energy( kernel );

    // integrator and assembler object
    typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                        typename SURFACETUPLEBINDER::TransposedTuple>
        SurfaceForceInt;
    typename SurfaceForceInt::ForceKernel surfaceForceKernel =
        boost::bind( &Energy::rhs, &energy, _1, _2, _3, _4 );
    
    SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                     surfaceQuadrature, solver );

    // apply
    typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
    typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
    for ( ; iter != end; ++iter ) {
                    
        const double sign  = (plusMinus ? 1.0 : -1.0 );
        const double kappa = sign * parameter.energyWeight( iter, inOut );
        energy.setKappa( kappa );
        
        surfaceForceInt( SURFACETUPLEBINDER::makeTransposedTuple( *iter ) );
    }

    return;
}

template<typename KERNEL, typename SURFFIELDTUPLE>
class Energy
{
public:
    //! @name Template parameter
    //@{
    typedef KERNEL         Kernel;
    typedef SURFFIELDTUPLE SurfFieldTuple;
    //@}

    //! @name Extract element types
    //@{
    typedef typename SurfFieldTuple::GeomElement  SurfaceElement;
    typedef typename SurfFieldTuple::TestElement  TestElement;
    typedef typename SurfFieldTuple::TrialElement TrialElement;
    //@}

    //! Type of field element tuple generated from the surface tuple
    typedef typename
    base::asmb::DomainFieldElementPointerTuple<SurfFieldTuple>::Type DomainFieldTuple;


    typedef typename SurfFieldTuple::TransposedTuple TransposedSurfFieldTuple;

    typedef typename
    base::asmb::DomainFieldElementPointerTuple<TransposedSurfFieldTuple>::Type
    TransposedDomainFieldTuple;

    //! Local coordinate
    typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim GlobalVecDim;

    //! Type of domain element
    typedef typename SurfaceElement::DomainElement    DomainElement;

    //! Size of a DoF
    static const unsigned doFSize = TestElement::DegreeOfFreedom::size;

    //! Type of DoF value vector
    typedef typename base::Vector<doFSize,base::number>::Type VecDof;

    //! Type of BC function
    typedef boost::function<VecDof( const GlobalVecDim& )> BCFun;

    //! Constructor with kernel object and multiplier
    Energy( const Kernel& kernel, const double kappa = 1.0 )
        : kernel_( kernel ), kappa_( kappa ) { }

    void setKappa( const double kappa ) { kappa_ = kappa; }

    //--------------------------------------------------------------------------
public:
    //--------------------------------------------------------------------------
    void rhs( const TransposedSurfFieldTuple& surfFieldTuple, 
              const LocalVecDim&              eta,
              const double                    weight,
              base::VectorD&                  result )
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();

        // Get surface metric
        GlobalVecDim normal;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, normal );
        
        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );

        // get domain field element pointer tuple (still transposed)
        const TransposedDomainFieldTuple domainFieldTupleT =
            base::asmb::DomainFieldElementPointerTuple<TransposedSurfFieldTuple>::
            convert( surfFieldTuple );

        // un-transpose the tuple
        const DomainFieldTuple domainFieldTuple =
            domainFieldTupleT.transpose();

        // compute co-normal from kernel using domain tuple
        base::MatrixD coNormal;
        kernel_.coNormalDerivative( domainFieldTuple, xi,
                                    normal, coNormal );

        // Evaluate solution
        const VecDof  u = base::post::evaluateField( surfEp -> getDomainElementPointer(),
                                                     surfFieldTuple.trialElementPtr(),
                                                     xi );

        // scalar multiplier
        const double scalar = -1. * detG * weight * kappa_;

        const unsigned otherSize = static_cast<unsigned>(coNormal.cols() );
        
        for ( unsigned i = 0; i < otherSize; i++ ) {
            double sum = 0.;
            for ( unsigned d = 0; d < doFSize; d++ ) {
                sum += scalar * coNormal( d, i ) * u[d];
            }
            result[i] += sum;
        }

        return;
    }

private:
    const Kernel&  kernel_; //!< Delivers conormal and dual conormal derivative
    double         kappa_;  //!< Factor for interface binding
};
