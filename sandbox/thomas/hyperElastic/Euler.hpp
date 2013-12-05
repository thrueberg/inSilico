
class AnalyticEuler1D
{
public:
    AnalyticEuler1D() { }

    void setFactor( const double ubar )
    {
        alpha_ = ubar;
        beta_  = 3. * alpha_ * alpha_ + 3. * alpha_ + 1.;
    }

    
public:
    double X( const double x ) const
    {
        const double arg = beta_*x + alpha_*alpha_*alpha_;
        return std::pow( arg, 1./3. ) - alpha_;
    }

    double u( const double x ) const
    {
        return x - X(x);
    }


    double f( const double x ) const
    {
        const double arg = beta_*x + alpha_*alpha_*alpha_;
        return beta_/3. * std::pow( arg, -2./3. );
    }

    double df( const double x ) const
    {
        const double arg = beta_*x + alpha_*alpha_*alpha_;
        return -2.*beta_*beta_/9. * std::pow( arg, -5./3. );
    }

private:
    double alpha_;
    double beta_;
};

//------------------------------------------------------------------------------
template<unsigned DIM>
class AnalyticEulerTensor
{
public:

    static const unsigned dim = DIM;
    typedef typename base::Vector<dim>::Type VecDim;
    
    AnalyticEulerTensor( const double lambda, const double mu,
                         boost::array<double,DIM> ubar )
        : lambda_( lambda ), mu_( mu )
    {
        for ( unsigned d = 0; d < DIM; d++ )
            oneDimensional_[d].setFactor( ubar[d] );
    }

public:
    mat::Tensor f( const VecDim& X ) const
    {
        mat::Tensor f = mat::Tensor::Identity();
        for ( unsigned d = 0; d < DIM; d++ )
            f(d,d) = oneDimensional_[d].f( X[d] );

        return f;
    }

    
    mat::Tensor sigma( const VecDim& x ) const
    {
        const mat::Tensor ff = f( x );
        mat::Tensor F;
        const double j = mat::inverseAndDet( ff, F );
        mat::Tensor b;
        mat::leftCauchyGreen( F, b );

        return (mu_*j) * b -
            (lambda_ * std::log(j) + mu_) * j * mat::Tensor::Identity();
    }

    VecDim force( const VecDim& x ) const
    {
        VecDim ff, df;
        double j = 1.;
        for ( unsigned d = 0; d < dim; d++ ) {
            ff[d] = oneDimensional_[d].f(  x[d] );
            df[d] = oneDimensional_[d].df( x[d] );
            j *= ff[d];
        }

        VecDim result;
        for ( unsigned d = 0; d < dim; d++ )
            result[d] = j * df[d] / ff[d] *
                ( lambda_ * (1. + std::log(j)) + mu_ * (1. + 1./ff[d]/ff[d]) );
        
        return result;
    }

    VecDim solution( const typename base::Vector<DIM>::Type x ) const
    {
        VecDim result;
        for ( unsigned d = 0; d < dim; d++ )
            result[d] = oneDimensional_[d].u( x[d] );
        
        return result;
    }

private:
    const double lambda_;
    const double mu_;

    boost::array<AnalyticEuler1D,dim> oneDimensional_;
};


//------------------------------------------------------------------------------
template<typename GEOMELEMENT, typename FIELDELEMENT>
void deformationGradient( const GEOMELEMENT*  geomEp,
                          const FIELDELEMENT* fieldEp,
                          const typename FIELDELEMENT::FEFun::VecDim& xi,
                          typename mat::Tensor& f )
{
    // Compute the deformation gradient
    f = mat::Tensor::Identity();

    // Evaluate the displacement gradient
    const typename base::Matrix<GEOMELEMENT::Node::dim,
                                FIELDELEMENT::DegreeOfFreedom::size,
                                double>::Type
        gradU = base::post::evaluateFieldGradient( geomEp, fieldEp, xi );
        
    // Add displacement gradient to tensor (note the transposition)
    f.block( 0, 0, gradU.cols(), gradU.rows() ) -= gradU.transpose();

    ASSERT_MSG( mat::determinant( f ) > 0., "J=det(f) > 0 is mandatory" );
}

//------------------------------------------------------------------------------
template<typename MATERIAL, typename FIELDTUPLE>
class HyperElasticEuler
{
public:
    typedef MATERIAL     Material;
    typedef FIELDTUPLE   FieldTuple;

    typedef typename FieldTuple::GeomElement      GeomElement;
    typedef typename FieldTuple::TestElement      TestElement;
    typedef typename FieldTuple::TrialElement     TrialElement;

    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Number of DoFs per vectorial entry
    static const unsigned nDoFs = TestElement::DegreeOfFreedom::size;


    HyperElasticEuler( Material& material )
        : material_( material )
    {
        
    }

    void tangentStiffness( const FieldTuple&     fieldTuple,
                           const LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of test and trial functions
        std::vector<GlobalVecDim> testGradX, trialGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );
        trialGradX = testGradX;

        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testGradX.size() );
        const unsigned numColBlocks = static_cast<unsigned>( trialGradX.size() );

        // Get deformation gradient of the current state
        typename mat::Tensor f;
        deformationGradient( geomEp, trialEp, xi, f );

        // Material evaluations
        mat::Tensor F; mat::inverse( f, F );
        mat::ElastTensor cElMixed;
        material_.mixedElasticityTensor( F, cElMixed );
        
        // loop over the test and trial functions
        for ( unsigned M = 0; M < numRowBlocks; M++ ) { // test functions
            for ( unsigned N = 0; N < numColBlocks; N++ ) { // trial functions

                // loop over the vector components of test and trials
                for ( unsigned i = 0; i < nDoFs; i++ ) { // test fun comp
                    for ( unsigned K = 0; K < nDoFs; K++ ) { // trial fun comp

                        double sum = 0.;
                        for ( unsigned j = 0; j < nDoFs; j++ ) {

                            // first Voigt index
                            const unsigned ij = mat::Voigt::apply( i, j );
                            
                            for ( unsigned L = 0; L < nDoFs; L++ ) {

                                // second Voigt index
                                const unsigned KL = mat::Voigt::apply( K, L );

                                // compute 'material gradient' \phi^N_{,m} F_{mL}
                                double materialGradient = 0.;
                                for ( unsigned m = 0; m < nDoFs; m++ ) {
                                    materialGradient += trialGradX[N][m] * F(m,L);
                                }

                                // \phi^M_{,j} c_{ijKL} (\phi^N_{,m} F_{mL}
                                sum +=
                                    testGradX[M][j] * cElMixed(ij, KL) * materialGradient;
                                
                            } // L
                        } // j

                        sum *= detJ * weight;//  * / J; //??
                        matrix( M * nDoFs + i, N * nDoFs + K ) += sum;
                        
                    } // K
                }// i

                
            } // N
        } // M

        return;

    }

    //--------------------------------------------------------------------------
    void residualForce( const FieldTuple&   fieldTuple,
                        const LocalVecDim&  xi,
                        const double        weight,
                        base::VectorD&      vector ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Evaluate gradient of test and trial functions
        std::vector<GlobalVecDim> testGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );

        const std::size_t numRowBlocks = testGradX.size();

        // Get deformation gradient of the current state
        typename mat::Tensor f;
        deformationGradient( geomEp, trialEp, xi, f );

        // Material evaluations
        typename mat::Tensor sigma;      // 2nd PK
        typename mat::Tensor F; mat::inverse( f, F );
        material_.cauchyStress( F, sigma );

        // loop over the test functions
        for ( std::size_t M = 0; M < numRowBlocks; M++ ) { // test functions

            // loop over the vector components of trials
            for ( unsigned i = 0; i < nDoFs; i++ ) { // test fun comp

                double sum = 0.;
                for ( unsigned j = 0; j < nDoFs; j++ )
                    sum += sigma( i, j ) * testGradX[M][j];

                sum *= detJ * weight;
                
                vector[ M * nDoFs + i ] += sum;
            }
        }

        return;
    }


private:
    const Material& material_;
};

