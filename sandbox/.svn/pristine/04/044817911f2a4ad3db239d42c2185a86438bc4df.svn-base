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
        typename mat::Tensor F;
        solid::deformationGradient( geomEp, trialEp, xi, F );

        // Material evaluations
        typename mat::Tensor S;      // 2nd PK
        typename mat::ElastTensor C; // elasticity tensor
        material_.secondPiolaKirchhoff( F, S );
        material_.materialElasticityTensor( F, C );

        // loop over the test and trial functions
        for ( unsigned M = 0; M < numRowBlocks; M++ ) { // test functions
            for ( unsigned N = 0; N < numColBlocks; N++ ) { // trial functions

                // loop over the vector components of test and trials
                for ( unsigned i = 0; i < nDoFs; i++ ) { // test fun comp
                    for ( unsigned k = 0; k < nDoFs; k++ ) { // trial fun comp

                        // Inner sum: \phi^M_{,J} C^{eff}_{iJkL} \phi^N_{,L}
                        double sum = 0.;
                        for ( unsigned J = 0; J < nDoFs; J++ ) {
                            for ( unsigned L = 0; L < nDoFs; L++ ) {

                                sum +=
                                    testGradX[M][J] *
                                    (this -> effectiveElasticity_( F, S, C, i, J, k, L ) ) *
                                    trialGradX[N][L];
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

        // 
        typename mat::Tensor Finc;  // (t,t+dt)
        solid::deformationGradient( geomEp, trialEp, xi, Finc);
        typename mat::Tensor FprevInv; // (t,0)
        inverseDeformationGradient( geomEp, trialEp, xi, FprevInv );
        typename mat::Tensor Fprev;    // (0,t)
        mat::inverse( FprevInv, Fprev );
        //solid::deformationGradientHistory<1>( geomEp, trialEp, xi, Fprev);

        typename mat::Tensor Ftot; // (0,t+dt)
        Ftot.noalias() = Finc * Fprev;

        const double J = mat::determinant( Fprev );
        
        // Material evaluations
        typename mat::Tensor Stot;      // 2nd PK
        material_.secondPiolaKirchhoff( Ftot, Stot );

        // First Piola-Kirchhoff stress tensor P = FS
        typename mat::Tensor Ptot;
        Ptot.noalias() = Ftot * Stot;

        typename mat::Tensor Pcurr;
        Pcurr.noalias() = 1./J * Ptot * Fprev.transpose();

        if ( geomEp -> getID() == 0 ) 
            std::cout << geomEp -> getID() << ":  "
                      << Ftot(0,0) << "  " << Stot(0,0) << "  "
                      << Ptot(0,0) << std::endl;

        // loop over the test functions
        for ( std::size_t M = 0; M < numRowBlocks; M++ ) { // test functions

            // loop over the vector components of trials
            for ( unsigned i = 0; i < nDoFs; i++ ) { // test fun comp

                double sum = 0.;
                for ( unsigned J = 0; J < nDoFs; J++ )
                    sum += Pcurr( i, J ) * testGradX[M][J];

                sum *= detJ * weight;
                
                vector[ M * nDoFs + i ] += sum;
            }
        }

        return;
    }

private:

    //--------------------------------------------------------------------------
    /** Computation of effective elasticity tensor entries.
     *  In the final expression of the linearised virtual strain energy, the
     *  so-called effective elasticity tensor pops up which is composed of
     *  the geometrical (or initial) stress contribution and the material
     *  contribution. In index notation, we have
     *  \f[
     *      C^{eff}_{iJkL} = \delta_{ik} S_{JL} + F_{iA} C_{AJBL} F_{kB}
     *  \f]
     *  using the second Piola-Kirchhoff stress tensor \f$ S \f$, the
     *  deformation gradient tensor \f$ F \f$, and the material elasticity
     *  tensor \f$ C \f$.
     *  Note that \f$ C \f$ is stored in the Voigt notation and therefore the
     *  access index pairs, e.g. \f$ 1 \leq A,J \leq 3 \f$, is converted to
     *  a Voigt index \f$ 1 \leq v \leq 6 \f$.
     *
     *  \param[in]  F        Deformation gradient
     *  \param[in]  S        2nd Piola-Kirchhoff stress tensor
     *  \param[in]  C        Material elasticity tensor
     *  \param[in]  i,J,k,L  Current indices of interest
     */
    double effectiveElasticity_( const typename mat::Tensor&      F,
                                 const typename mat::Tensor&      S,
                                 const typename mat::ElastTensor& C,
                                 const unsigned i, const unsigned J,
                                 const unsigned k, const unsigned L ) const
    {
        // geometrical stress contribution
        double result = ( i == k ? S( J, L ) : 0. );

        // material contribution
        
        for ( unsigned A = 0; A < nDoFs; A++ ) { 
            // first Voigt index
            const unsigned voigt1 = mat::Voigt::apply( A, J );
            
            for ( unsigned B = 0; B < nDoFs; B ++ ) {
                // second Voigt index
                const unsigned voigt2 = mat::Voigt::apply( B, L );
                
                // material elasticity 
                const double cAJBL = C( voigt1, voigt2 );

                // convert with deformation gradient
                result += F( i, A ) * cAJBL * F( k, B );
            }
        }

        return result;
    }

public:
    //--------------------------------------------------------------------------
    /** The linearised co-normal derivative of field functions.
     *  In reference configuration, the co-normal derivative becomes
     *  \f[
     *        B_N(u) = P(u) \cdot N
     *  \f]
     *  with the first Piola-Kirchhoff stress tensor \f$ P \f$ and the normal
     *  vector in reference configuration \f$ N \f$. Applying a Newton method,
     *  the directional derivative becomes
     *  \f[
     *       (C^{eff}(u^n) : \nabla_X \Delta u) \cdot N
     *  \f]
     *  i.e. the effective elasticity tensor evaluated for the current
     *  deformation state, contracted with the material gradient of the
     *  displacement increment (the direction of the derivative) and multiplied
     *  by the normal vector.
     *  The discrete counter-part gives the matrix
     *  \f[
     *      B[i,Md+k] = C^{eff}_{iJkL}(u^n) \phi^M_{,L} N_J
     *  \f]
     *  Since the dual co-normal derivative is also frequently needed, the
     *  implementation is deferred to normalDerivativeHelper_() .
     *  \param[in]   fieldTuple Tuple of field element pointers
     *  \param[in]   xi         Local evaluation coordinate
     *  \param[in]   normal     Normal vector in reference coordinates.
     *  \param[out]  result     Result container
     */
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


        // Get deformation gradient of the current state
        typename mat::Tensor F;
        solid::deformationGradient( geomEp, trialEp, xi, F );

        // Material evaluations
        typename mat::Tensor S;      // 2nd PK
        typename mat::ElastTensor C; // elasticity tensor
        material_.secondPiolaKirchhoff( F, S );
        material_.materialElasticityTensor( F, C );

        // loop over field function gradients (test or trial)
        for ( unsigned M = 0; M < numColBlocks; M++ ) {
            // loop over rows
            for ( unsigned i = 0; i < nDoFs; i++ )  {
                // loop over column degrees of freedom
                for ( unsigned k = 0; k < nDoFs; k++ ) {
                    
                    double entry = 0.;
                    // contract effective elastity with gradient and normal
                    for ( unsigned J = 0; J < nDoFs; J++ ) {
                        for ( unsigned L = 0; L < nDoFs; L++ ) {

                            entry +=
                                (this -> effectiveElasticity_(F, S, C, i, J, k, L) ) *
                                trialGradX[M][L] * normal[J];
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
    /**
     */
    void coNormalDerivativeResidual( const FieldTuple&   fieldTuple,
                                     const LocalVecDim&  xi,
                                     const GlobalVecDim& normal,
                                     base::VectorD&      result ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        // Get deformation gradient of the current state
        typename mat::Tensor F;
        solid::deformationGradientHistory<0>( geomEp, trialEp, xi, F );
        
        // Material evaluations
        typename mat::Tensor S;      // 2nd PK
        material_.secondPiolaKirchhoff( F, S );

        // First Piola-Kirchhoff stress tensor P = FS
        typename mat::Tensor P;
        P.noalias() = F * S;

        // Evaluate pressure trial functions
        typename TestElement::FEFun::FunArray testFun;
        (testEp -> fEFun()).evaluate( geomEp, xi, testFun );

        const std::size_t numRowBlocks = testFun.size();

        result = base::VectorD::Zero( numRowBlocks * nDoFs );

        for ( std::size_t M = 0; M < testFun.size(); M++ ) {
            for ( unsigned i = 0; i < nDoFs; i++) {

                double entry = 0.;
                for ( unsigned J = 0; J < nDoFs; J++ )
                    entry += P(i,J) * normal[J];

                result[ M * nDoFs + i ] = entry;
            }
        }
        
    }

    
private:
    const Material&    material_;   //!< Material behaviour
};

#endif
