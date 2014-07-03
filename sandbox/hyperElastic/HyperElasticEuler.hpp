//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   HyperElasticEuler.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef solid_hyperelasticeuler_hpp
#define solid_hyperelasticeuler_hpp

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
//#include <solid/Deformation.hpp>

#include <mat/hypel/SpatialWrapper.hpp>

//------------------------------------------------------------------------------
template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
void deformationGradientHistory( const GEOMELEMENT*  geomEp,
                                 const FIELDELEMENT* fieldEp,
                                 const typename FIELDELEMENT::FEFun::VecDim& xi,
                                 typename mat::Tensor& F )
{
    // Compute the deformation gradient
    mat::Tensor Finv = mat::Tensor::Identity();

    // Evaluate the displacement gradient
    const typename base::Matrix<GEOMELEMENT::Node::dim,
                                FIELDELEMENT::DegreeOfFreedom::size,
                                double>::Type
        gradU = base::post::evaluateFieldGradientHistory<0>( geomEp,
                                                             fieldEp, xi );
    
    // Evaluate the displacement gradient
    const typename base::Matrix<GEOMELEMENT::Node::dim,
                                FIELDELEMENT::DegreeOfFreedom::size,
                                double>::Type
        gradU2 = base::post::evaluateFieldGradientHistory<1>( geomEp,
                                                              fieldEp, xi );

        
    // Add displacement gradient to tensor
    Finv.block( 0, 0, gradU.rows(), gradU.cols() ) -= (gradU + gradU2).transpose();
    // first part, i.e. Grad(Delta u) shall be zero after geometry update

    mat::inverse( Finv, F );

    VERIFY_MSG( mat::determinant( F ) > 0., "J=det(F) > 0 is mandatory" );
}

//------------------------------------------------------------------------------
template<typename GEOMELEMENT, typename TRIALELEMENT, typename MATERIAL>
mat::Tensor cauchy( const GEOMELEMENT*  geomEp,
                    const TRIALELEMENT* trialEp,
                    const MATERIAL& material,
                    const typename TRIALELEMENT::FEFun::VecDim& xi )
{
    // get deformation gradient
    mat::Tensor F;
    ::deformationGradientHistory<0>( geomEp, trialEp, xi, F );
    
    // Material evaluations
    typename mat::Tensor sigma; 
    material.cauchy( F, sigma );
    
    return sigma;
}

//------------------------------------------------------------------------------
//! Overload call of solid::cauchy for element centroid
template<typename GEOMELEMENT, typename TRIALELEMENT, typename MATERIAL>
mat::Tensor cauchy( const GEOMELEMENT*  geomEp,
                    const TRIALELEMENT* trialEp,
                    const MATERIAL& material )
{
    return cauchy<GEOMELEMENT,TRIALELEMENT,MATERIAL>(
        geomEp, trialEp, material,
        base::ShapeCentroid<GEOMELEMENT::shape>::apply() );
}



//------------------------------------------------------------------------------

template<typename MATERIAL, typename FIELDTUPLE>
class HyperElasticEuler
    : boost::noncopyable
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

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    typedef mat::hypel::SpatialWrapper<Material> SpatialWrapper;
    

    //! Constructor with form and test functions 
    HyperElasticEuler( const Material&   material )
        : material_(   material ), bla_( material )
    { }

    //--------------------------------------------------------------------------
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
        //solid::deformationGradientHistory<1>( geomEp, trialEp, xi, F );
        ::deformationGradientHistory<0>( geomEp, trialEp, xi, F );

        // Material evaluations
        typename mat::Tensor sigma;    // Cauchy
        typename mat::ElastTensor cEl; // elasticity tensor
        //material_.cauchy( F, sigma );
        //material_.spatialElasticityTensor( F, cEl );
        material_.secondPiolaKirchhoff( F, sigma );
        material_.materialElasticityTensor( F, cEl );
        

        // loop over the test and trial functions
        for ( unsigned M = 0; M < numRowBlocks; M++ ) { // test functions
            for ( unsigned N = 0; N < numColBlocks; N++ ) { // trial functions

                // loop over the vector components of test and trials
                for ( unsigned i = 0; i < nDoFs; i++ ) { // test fun comp
                    for ( unsigned k = 0; k < nDoFs; k++ ) { // trial fun comp


                        double sum = 0.;
                        for ( unsigned j = 0; j < nDoFs; j++ ) {
                            for ( unsigned el = 0; el < nDoFs; el++ ) {

                                sum +=
                                    testGradX[M][j] *
                                    (this -> effectiveElasticity2_( F, sigma, cEl,
                                                                    i, j, k, el ) ) *
                                    trialGradX[N][el];
                            }
                        }

                        sum *= detJ * weight;
                        matrix( M * nDoFs + i, N * nDoFs + k ) += sum;
                        
                    } // j
                }// i
            } // N
        } // M

        //std::cout << sigma << std::endl << mat::determinant(F) << std::endl << std::endl;
        //std::cout << cEl << std::endl;
        
        return;
    }

    //--------------------------------------------------------------------------
    /** Internal force computation for the latest displacement field.
     *  Delegates call to residualForceHistory() with HIST = 0.
     */
    void residualForce( const FieldTuple&   fieldTuple,
                        const LocalVecDim&  xi,
                        const double        weight,
                        base::VectorD&      vector ) const
    {
        this -> residualForceHistory<0>( fieldTuple,
                                         xi, weight, vector );
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST>
    void residualForceHistory( const FieldTuple&   fieldTuple,
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

        // Get deformation gradient of the current state
        typename mat::Tensor F;
        //solid::deformationGradientHistory<HIST>( geomEp, trialEp, xi, F );
        solid::deformationGradientHistory<1>( geomEp, trialEp, xi, F );
        //::deformationGradientHistory<0>( geomEp, trialEp, xi, F );
        
        // Material evaluations
        typename mat::Tensor sigma; // cauchy
        material_.cauchy( F, sigma );

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
    double effectiveElasticity2_( const typename mat::Tensor&      F,
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

        //std::cout << "(" << i << J << k << L << "): " << result  << ", " << S(J,L) << std::endl;
        //std::cout << "F = " << F << "\n\n" << "S = " << S << "\n\n"
        //          << "C = " << C << "\n\n" << result << std::endl;

        return result;
    }


    //--------------------------------------------------------------------------
    double effectiveElasticity_( const typename mat::Tensor&      F,
                                 const typename mat::Tensor&      sigma,
                                 const typename mat::ElastTensor& cEl,
                                 const unsigned i, const unsigned j,
                                 const unsigned k, const unsigned el ) const
    {
        const unsigned voigt1 = mat::Voigt::apply( i, j );
        const unsigned voigt2 = mat::Voigt::apply( k, el );
        
        // geometrical stress contribution
        const double result =
            ( i == k ? sigma( j, el ) : 0. ) + cEl( voigt1, voigt2 );

        //std::cout << "(" << i << j << k << el << "): " << result
        //          << "  --- " << F(0,0) << " -> " << sigma(0,0)
        //          << std::endl;
        
        return result;
    }

    
private:
    const SpatialWrapper material_;

    const Material bla_;
};

#endif
