#ifndef hyperelastic_hpp
#define hyperelastic_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
#include <base/auxi/EqualPointers.hpp>
// mat includes
#include <mat/TensorAlgebra.hpp>
// solid includes
#include <solid/Deformation.hpp>

//------------------------------------------------------------------------------
template<typename MATERIAL,
         typename FIELDTUPLE>
class HyperElastic;

//--------------------------------------------------------------------------
template<typename GEOMELEMENT, typename FIELDELEMENT>
void inverseDeformationGradient( const GEOMELEMENT*  geomEp,
                                 const FIELDELEMENT* fieldEp,
                                 const typename FIELDELEMENT::FEFun::VecDim& xi,
                                 typename mat::Tensor& Finv )
{
    // Compute the deformation gradient
    Finv = mat::Tensor::Identity();

    // Evaluate the displacement gradient
    const typename base::Matrix<GEOMELEMENT::Node::dim,
                                FIELDELEMENT::DegreeOfFreedom::size,
                                double>::Type
        GradU = base::post::evaluateFieldGradientHistory<1>( geomEp,
                                                             fieldEp, xi );
    
    // Add displacement gradient to tensor (note the transposition)
    Finv.block( 0, 0, GradU.cols(), GradU.rows() ) -= GradU.transpose();

    
    ASSERT_MSG( mat::determinant( Finv ) > 0., "J=det(F) > 0 is mandatory" );
}


//------------------------------------------------------------------------------
template<typename MATERIAL, typename FIELDTUPLE>
class HyperElastic
{
public:
    //! @name Template parameters
    //@{
    typedef MATERIAL     Material;
    typedef FIELDTUPLE   FieldTuple;
    //@}

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement      GeomElement;
    typedef typename FieldTuple::TestElement      TestElement;
    typedef typename FieldTuple::TrialElement     TrialElement;
    //@}

    //! Number of DoFs per vectorial entry
    static const unsigned nDoFs = TestElement::DegreeOfFreedom::size;

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //! Type of global vector
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim GlobalVecDim;

    //! Constructor with form and test functions 
    HyperElastic( const Material&   material )
       : material_(   material )
    { }

    //--------------------------------------------------------------------------
    /**  Contribution to the element stiffness matrix in a quadrature rule.
     *   The element stiffness matrix for the 
     *   \f[
     *       K[M d+i, N d+k] = \int \phi^M_{,J} C^{eff}_{iJkL} \phi^N_{,L} dX
     *   \f]
     *   This object adds the weighted integrand evaluated at a local coordinate
     *   \f$\xi\f$ to a provided storage.
     *   \param[in]  fieldTuple Tuple of field element pointers
     *   \param[in]  xi       Local coordinate: quadrature point
     *   \param[in]  weight   Weight corresponding to the quadrature point
     *   \param[out] matrix   Result storage
     */
    void tangentStiffness( const FieldTuple&     fieldTuple,
                           const LocalVecDim&    xi,
                           const double          weight,
                           base::MatrixD&        matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // bubnov flag
        const bool bubnov =
            base::auxi::EqualPointers<TestElement,TrialElement>::apply( testEp,
                                                                        trialEp );
        
        // Evaluate gradient of test and trial functions
        std::vector<GlobalVecDim> testGradX, trialGradX;
        const double detJ =
            (testEp -> fEFun()).evaluateGradient( geomEp, xi, testGradX );
        
        if ( bubnov ) trialGradX = testGradX;
        else
            (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );

        // Sizes and sanity checks
        const unsigned numRowBlocks = static_cast<unsigned>( testGradX.size() );
        const unsigned numColBlocks = static_cast<unsigned>( trialGradX.size() );
        assert( static_cast<unsigned>( matrix.rows() ) == numRowBlocks * nDoFs );
        assert( static_cast<unsigned>( matrix.cols() ) == numColBlocks * nDoFs );

        // Get deformation gradient of the current state
        mat::Tensor H;  // (t,t+dt)
        solid::deformationGradient( geomEp, trialEp, xi, H );
        mat::Tensor g; // (t,0)
        inverseDeformationGradient( geomEp, trialEp, xi, g );
        mat::Tensor G;    // (0,t)
        mat::inverse( g, G );
        mat::Tensor F; // (0,t+dt)
        F.noalias() = H * G;

        const double blub = mat::determinant( F );
        if ( blub < 0.01 ) {
            Eigen::IOFormat bla( Eigen::StreamPrecision, 0, ", ", ", " );
            std::cout << geomEp -> getID() << ": "
                      << blub << "| " << xi.transpose() << " -- "
                      << (base::Geometry<GeomElement>()( geomEp, xi ) ).transpose() << "\n"
                      << "F: " << F.format( bla ) << "\n"
                      << "H: " << H.format( bla ) << "\n"
                      << "G: " << G.format( bla ) << "\n"
                      << "\n\n";

            return;
        }
        

        mat::Tensor sigma;
        mat::ElastTensor Csp;
        material_.cauchyStress(            F, sigma );
        material_.spatialElasticityTensor( F, Csp );
        
        // loop over the test and trial functions
        for ( unsigned M = 0; M < numRowBlocks; M++ ) { // test functions
            for ( unsigned N = 0; N < numColBlocks; N++ ) { // trial functions

                // loop over the vector components of test and trials
                for ( unsigned i = 0; i < nDoFs; i++ ) { // test fun comp
                    for ( unsigned k = 0; k < nDoFs; k++ ) { // trial fun comp

                        // Inner sum: \phi^M_{,j} C^{eff}_{ijkl} \phi^N_{,l}
                        double sum = 0.;
                        for ( unsigned alpha = 0; alpha < nDoFs; alpha++ ) {
                            for ( unsigned beta = 0; beta < nDoFs; beta++ ) {

                                sum +=
                                    testGradX[M][alpha] *
                                    (this -> effectiveElasticity_( H, sigma, Csp,
                                                                   i, alpha, k, beta ) ) *
                                    trialGradX[N][beta];
                                
                            }
                        }

                        sum *= detJ * weight;
                        matrix( M * nDoFs + i, N * nDoFs + k ) += sum;
                        
                    } // j
                }// i
            } // N
        } // M

        return;
    }

    //--------------------------------------------------------------------------
    /** Internal force computation for the latest displacement field.
     */
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
        assert( static_cast<unsigned>( vector.size() ) == numRowBlocks * nDoFs );

        // incremental deformation gradients
        mat::Tensor H;  // (t,t+dt)
        solid::deformationGradient( geomEp, trialEp, xi, H);
        mat::Tensor g; // (t,0)
        inverseDeformationGradient( geomEp, trialEp, xi, g );
        mat::Tensor G;    // (0,t)
        mat::inverse( g, G );

        // total deformation gradient
        mat::Tensor F; // (0,t+dt)
        F.noalias() = H * G;


        // Material evaluations
        mat::Tensor sigma;
        material_.cauchyStress( F, sigma );
        const double detH = mat::determinant( H );
        mat::Tensor Hinv;
        mat::inverse( H, Hinv );
        // shift to latest configuration
        sigma = detH * sigma * Hinv.transpose();
        
        // loop over the test functions
        for ( std::size_t M = 0; M < numRowBlocks; M++ ) { // test functions

            // loop over the vector components of trials
            for ( unsigned i = 0; i < nDoFs; i++ ) { // test fun comp

                double sum = 0.;
                for ( unsigned alpha = 0; alpha < nDoFs; alpha++ )
                    sum += sigma(i,alpha) * testGradX[M][alpha];

                sum *= detJ * weight;
                
                vector[ M * nDoFs + i ] += sum;
            }
        }

        return;
    }

private:
    //--------------------------------------------------------------------------
    //! Effective elasticity coefficient in updated reference configuration
    double effectiveElasticity_( const typename mat::Tensor&      H,
                                 const typename mat::Tensor&      sigma,
                                 const typename mat::ElastTensor& Csp,
                                 const unsigned i, const unsigned alpha,
                                 const unsigned k, const unsigned beta ) const
    {
        mat::Tensor Hinv;
        const double detH = mat::inverseAndDet( H, Hinv );

        double c_iakb = 0.;

        for ( unsigned j = 0; j < nDoFs; j++ ) {
            const unsigned ij = mat::Voigt::apply( i, j );
            
            for ( unsigned el = 0; el < nDoFs; el++ ) {
                const unsigned kl = mat::Voigt::apply( k, el );

                const double cEff_ijkl =
                    (i==k ? sigma(j,el) : 0.) + Csp( ij, kl );

                c_iakb += detH * Hinv(alpha,j) * cEff_ijkl * Hinv(beta,el);

            }
        }

        return c_iakb;
    }


public:
    //--------------------------------------------------------------------------
    void coNormalDerivative( const FieldTuple&  fieldTuple,
                             const LocalVecDim& xi,
                             const GlobalVecDim& normal,
                             base::MatrixD& result ) const
    {
        // Extract element pointer from tuple
        const GeomElement*   geomEp  = fieldTuple.geomElementPtr();
        const TrialElement*  trialEp = fieldTuple.trialElementPtr();

        // evaluate gradient of trial functions
        std::vector<GlobalVecDim> trialGradX;
        (trialEp -> fEFun()).evaluateGradient( geomEp, xi, trialGradX );

        const unsigned numColBlocks = static_cast<unsigned>( trialGradX.size() );

        result = base::MatrixD::Zero( +nDoFs, numColBlocks * nDoFs );

        // incremental deformation gradients
        mat::Tensor H;  // (t,t+dt)
        solid::deformationGradient( geomEp, trialEp, xi, H);
        mat::Tensor g; // (t,0)
        inverseDeformationGradient( geomEp, trialEp, xi, g );
        mat::Tensor G;    // (0,t)
        mat::inverse( g, G );

        // total deformation gradient
        mat::Tensor F; // (0,t+dt)
        F.noalias() = H * G;

        // Material evaluations
        mat::Tensor sigma;
        mat::ElastTensor Csp;
        material_.cauchyStress(            F, sigma );
        material_.spatialElasticityTensor( F, Csp );

        // loop over field function gradients (test or trial)
        for ( unsigned M = 0; M < numColBlocks; M++ ) {
            // loop over rows
            for ( unsigned i = 0; i < nDoFs; i++ )  {
                // loop over column degrees of freedom
                for ( unsigned k = 0; k < nDoFs; k++ ) {
                    
                    double entry = 0.;
                    // contract effective elastity with gradient and normal
                    for ( unsigned alpha = 0; alpha < nDoFs; alpha++ ) {
                        for ( unsigned beta = 0; beta < nDoFs; beta++ ) {

                            entry +=
                                (this -> effectiveElasticity_(H, sigma, Csp, i, alpha, k, beta) ) *
                                trialGradX[M][beta] * normal[alpha];
                        }
                    }
                    result( i, M * nDoFs + k ) = entry;
                }
            }
        }
        
        return;
    }

public:
//--------------------------------------------------------------------------
    void boundaryResidual( const FieldTuple&   fieldTuple,
                           const LocalVecDim&  xi,
                           const GlobalVecDim& normal,
                           base::VectorD&      result ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // incremental deformation gradients
        mat::Tensor H;  // (t,t+dt)
        solid::deformationGradient( geomEp, trialEp, xi, H);
        mat::Tensor g; // (t,0)
        inverseDeformationGradient( geomEp, trialEp, xi, g );
        mat::Tensor G;    // (0,t)
        mat::inverse( g, G );

        // total deformation gradient
        mat::Tensor F; // (0,t+dt)
        F.noalias() = H * G;


        // Material evaluations
        mat::Tensor sigma;
        material_.cauchyStress( F, sigma );
        const double detH = mat::determinant( H );
        mat::Tensor Hinv;
        mat::inverse( H, Hinv );
        // shift to latest configuration
        sigma = detH * sigma * Hinv.transpose();

        // Evaluate test functions
        typename TestElement::FEFun::FunArray testFun;
        (testEp -> fEFun()).evaluate( geomEp, xi, testFun );

        const std::size_t numRowBlocks = testFun.size();

        // initialise result container
        result = base::VectorD::Zero( numRowBlocks * nDoFs );

        for ( std::size_t M = 0; M < testFun.size(); M++ ) {
            for ( unsigned i = 0; i < nDoFs; i++) {

                double ti = 0.;
                for ( unsigned alpha = 0; alpha < nDoFs; alpha++ )
                    ti += sigma(i, alpha) * normal[alpha];

                result[ M * nDoFs + i ] = ti * testFun[M];
            }
        }

        return;
    }

private:
    const Material&    material_;   //!< Material behaviour
};

//--------------------------------------------------------------------------
template<typename FIELDTUPLE, typename MATERIAL>
class Traction
    : public boost::function<
    typename base::Vector<FIELDTUPLE::GeomElement::Node::dim>::Type(
        const FIELDTUPLE&,
        const typename base::GeomTraits<typename
        FIELDTUPLE::GeomElement>::LocalVecDim&,
        const typename base::GeomTraits<typename
        FIELDTUPLE::GeomElement>::GlobalVecDim&) >
{
public:
    //! Template parameter: a field tuple
    typedef FIELDTUPLE FieldTuple;

    //! @name Derive element pointer types
    //@{
    typedef typename FieldTuple::GeomElement              GeomElement;
    typedef typename FieldTuple::template Binder<1>::Type DisplacementElementPtr;
    //@}
    
    //! Spatial dimension
    static const unsigned dim = GeomElement::Node::dim;
    
    //! @name Result and coordinate types
    //@{ 
    typedef typename base::Matrix<dim,dim>::Type                MatrixType;
    typedef typename base::Vector<dim>::Type                    ResultType;
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim GlobalVecDim;
    //@}
    
    //! Constructor with viscosity
    Traction( const MATERIAL& material ) : material_( material ) { }

    //----------------------------------------------------------------------
    ResultType operator()( const FieldTuple& fieldTuple,
                           const LocalVecDim& xi,
                           const GlobalVecDim& normal ) const
    {
        // extract field pointers
        const GeomElement*       geomEp     = fieldTuple.geomElementPtr();
        const DisplacementElementPtr dispEp = fieldTuple.template get<1>();

        // incremental deformation gradients
        mat::Tensor H;  // (t,t+dt)
        solid::deformationGradient( geomEp, dispEp, xi, H);
        mat::Tensor g; // (t,0)
        inverseDeformationGradient( geomEp, dispEp, xi, g );
        mat::Tensor G;    // (0,t)
        mat::inverse( g, G );

        // total deformation gradient
        mat::Tensor F; // (0,t+dt)
        F.noalias() = H * G;

        // Material evaluations
        mat::Tensor sigma;
        material_.cauchyStress( F, sigma );

        // compute resulting Cauchy stress
        MatrixType s;
        for ( unsigned d1 = 0; d1 < dim; d1++ ) {
            for ( unsigned d2 = 0; d2 < dim; d2++ ) {

                s(d1,d2) = sigma( d1, d2 );
            }
        }

        ResultType traction;
        traction.noalias() = s * normal;
        return traction;
    }
        
private:
    const MATERIAL material_;
};

#endif
