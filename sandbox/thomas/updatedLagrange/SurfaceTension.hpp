#ifndef surfacetension_h
#define surfacetension_h

template<typename SURFFIELDTUPLE> class SurfaceTension;


//----------------------------------------------------------------------
template<typename SURFACETUPLEBINDER,
         typename SURFACEQUADRATURE,
         typename SOLVER,     typename BOUNDFIELD,
         typename PARAMETER>
void surfaceTension( const SURFACEQUADRATURE& surfaceQuadrature,
                     SOLVER&                  solver, 
                     const BOUNDFIELD&        boundField,
                     const PARAMETER&         parameter,
                     const double             sigma, 
                     const bool               inOut = true,
                     const bool               plusMinus = true )
{
    // object to compute the LHS penalty term
    typedef SurfaceTension<typename SURFACETUPLEBINDER::Tuple> ST;
    ST surfaceTension( sigma );

    // integrator and assembler object
    typedef base::asmb::ForceIntegrator<SURFACEQUADRATURE,SOLVER,
                                        typename SURFACETUPLEBINDER::Tuple>
        SurfaceForceInt;
    typename SurfaceForceInt::ForceKernel surfaceForceKernel =
        boost::bind( &ST::residual, &surfaceTension, _1, _2, _3, _4 );
    
    SurfaceForceInt surfaceForceInt( surfaceForceKernel,
                                     surfaceQuadrature, solver );
            
    // apply
    typename BOUNDFIELD::FieldIterator iter = boundField.elementsBegin();
    typename BOUNDFIELD::FieldIterator end  = boundField.elementsEnd();
    for ( ; iter != end; ++iter ) {
        
        const double sign  = (plusMinus ? 1.0 : -1.0 );
        const double kappa = sign * parameter.energyWeight( iter, not inOut ); //!<!!!
        surfaceTension.setKappa( kappa );
        surfaceForceInt( SURFACETUPLEBINDER::makeTuple( *iter ) );
    }
    
    return;
}

//------------------------------------------------------------------------------
template<typename SURFFIELDTUPLE>
class SurfaceTension
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

    //! Type of field element tuple generated from the surface tuple
    typedef typename
    base::asmb::DomainFieldElementPointerTuple<SurfFieldTuple>::Type DomainFieldTuple;

    //! Local coordinate
    typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim GlobalVecDim;

    //! Type of domain element
    typedef typename SurfaceElement::DomainElement    DomainElement;

    //! Dimension
    static const unsigned globalDim = base::GeomTraits<SurfaceElement>::globalDim;


    //--------------------------------------------------------------------------
    SurfaceTension( const double sigma, const double kappa = 1.0 )
        : sigma_( sigma ), kappa_( kappa ) { }

    void setKappa( const double kappa ) { kappa_ = kappa; }

    //--------------------------------------------------------------------------
    void residual( const SurfFieldTuple& surfFieldTuple,
                   const LocalVecDim&    eta,
                   const double          weight,
                   base::VectorD&        result )
    {
        // extract test and trial elements from tuple
        const SurfaceElement* surfEp  = surfFieldTuple.geomElementPtr();

        const TestElement* testEp     = surfFieldTuple.testElementPtr();

        // Get surface metric
        GlobalVecDim normal;
        const double detG =
            base::SurfaceNormal<SurfaceElement>()( surfEp, eta, normal );
        
        // Tangential projector
        typename base::Matrix<globalDim,globalDim>::Type P;
        for ( unsigned i = 0; i < globalDim; i++ )
            for ( unsigned j = 0; j < globalDim; j++ )
                P(i,j) = (i==j ? 1. : 0.) - normal[i] * normal[j];


        // Get local domain coordinate
        typename DomainElement::GeomFun::VecDim xi =
            surfEp -> localDomainCoordinate( eta );

        // Standard gradient of test functions
        std::vector<GlobalVecDim> testGradX;
        (testEp -> fEFun()).evaluateGradient( surfEp -> getDomainElementPointer(),
                                              xi, testGradX );

        // Tangential gradient of identity
        typename base::Matrix<globalDim,globalDim>::Type gradId =
            base::constantMatrix<globalDim,globalDim>( 0. );

        // nodal coordinates of the surface element
        boost::array<GlobalVecDim,SurfaceElement::numNodes> nodalX =
            base::NodalCoordinates<SurfaceElement>()( surfEp );
        // Tangential gradients of geometry
        std::vector<GlobalVecDim> geomGradX;
        (surfEp -> geomFun()).evaluateGradient( surfEp, eta, geomGradX );

        for ( unsigned L = 0; L < SurfaceElement::numNodes; L++ )
            gradId += nodalX[L] * (geomGradX[L].transpose());

        //
        const double scalar = -detG * weight * kappa_ * sigma_;

        for ( unsigned K = 0; K < testGradX.size(); K++ ) {

            // tangential derivative of test-gradient
            GlobalVecDim testGradTan;
            testGradTan.noalias() = P * testGradX[K];

            for ( unsigned i = 0; i < globalDim; i++ ) {

                double entry = 0.;
                for ( unsigned j = 0; j < globalDim; j++ )
                    entry += testGradTan[j] * gradId( i, j );

                entry *= scalar;

                result[ K * globalDim + i ] += entry;
            }
        }
        
        return;
    }

    

private:
    const double sigma_; //!< Physical value of the surface tension
    double       kappa_; //!< Method parameter for interface conditions
};


#endif
