#ifndef stabil_hpp
#define stabil_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/geometry.hpp>
#include <base/linearAlgebra.hpp>
// base/auxi includes
#include <base/auxi/EqualPointers.hpp>
// base/kernel includes
#include <base/kernel/KernelFun.hpp>

#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace detail_{

    template<typename FIELDTUPLE, typename QUAD>
    struct CollectQuadraturePoints
        : public boost::function< void( const FIELDTUPLE&,
                                        const typename QUAD::VecDim&,
                                        const double,
                                        std::vector< std::pair<double,
                                        typename QUAD::VecDim> >& ) >
    {
        typedef typename QUAD::VecDim LocalVecDim;
        typedef std::pair<double,LocalVecDim> WeightedPoint;

        void operator()( const FIELDTUPLE&  fieldTuple, 
                         const LocalVecDim& xi,
                         const double       weight,
                         std::vector<WeightedPoint>& weightedPoints ) const
        {
            weightedPoints.push_back( std::make_pair( weight, xi ) );
        }

    };
    

}


//------------------------------------------------------------------------------
template<typename FIELDTUPLE,
         typename QUADRATURE,
         unsigned DOFSIZE =
         FIELDTUPLE::TrialElement::DegreeOfFreedom::size>
class Stabil1
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;
    typedef QUADRATURE Quadrature;
    static const unsigned doFSize   = DOFSIZE;
    STATIC_ASSERT_MSG( doFSize == 1, "Only implement for pressure" );

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    typedef void result_type;

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    //--------------------------------------------------------------------------
    void operator()( const FieldTuple&  fieldTuple,
                     const Quadrature&  quadrature, 
                     base::MatrixD&     matrix ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();

        //----------------------------------------------------------------------
        // compute averages
        double sizeElem = 0.;

        boost::array<double,TrialElement::FEFun::numFun> trialAvg;
        trialAvg.assign( 0. );
        boost::array<double, TestElement::FEFun::numFun>  testAvg;
        testAvg.assign( 0. );


        // collect quadrature points
        typedef detail_::CollectQuadraturePoints<FieldTuple,Quadrature> Collector;
        Collector collector;
        typedef std::vector<typename Collector::WeightedPoint> WeightedPoints;
        WeightedPoints weightedPoints;
        quadrature.apply( collector, fieldTuple, weightedPoints );
        

        typename WeightedPoints::iterator qIter = weightedPoints.begin();
        typename WeightedPoints::iterator qEnd  = weightedPoints.end();
        for ( ; qIter != qEnd; ++qIter ) {

            const double      weight = qIter -> first;
            const LocalVecDim xi     = qIter -> second;
            
            // element jacobian
            const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );

            // element size
            sizeElem += weight * detJ;

            // Evaluate trial functions
            typename TrialElement::FEFun::FunArray trialFun;
            (trialEp ->  fEFun()).evaluate( geomEp, xi, trialFun );

            // Evaluate test functions
            typename TestElement::FEFun::FunArray  testFun;
            (testEp -> fEFun()).evaluate( geomEp, xi, testFun );

            for ( std::size_t t = 0; t < trialFun.size(); t++ )
                trialAvg[t] += trialFun[t] * weight * detJ;
            
            for ( std::size_t t = 0; t < trialFun.size(); t++ )
                testAvg[t]  +=  testFun[t] * weight * detJ;

        }


        // compute the matrix entries
        qIter = weightedPoints.begin();
        qEnd  = weightedPoints.end();
        for ( ; qIter != qEnd; ++qIter ) {
        
            const double      weight = qIter -> first;
            const LocalVecDim xi     = qIter -> second;
            
            // element jacobian
            const double detJ = base::Jacobian<GeomElement>()( geomEp, xi );
        
            // Evaluate trial functions
            typename TrialElement::FEFun::FunArray trialFun;
            (trialEp ->  fEFun()).evaluate( geomEp, xi, trialFun );

            // Evaluate test functions
            typename TestElement::FEFun::FunArray  testFun;
            (testEp -> fEFun()).evaluate( geomEp, xi, testFun );

            // Sizes 
            const unsigned numRowBlocks = static_cast<unsigned>(  testFun.size() );
            const unsigned numColBlocks = static_cast<unsigned>( trialFun.size() );

            //matrix += scalar * (testGradX.transpose() * formGradX);
            for ( unsigned M = 0; M < numRowBlocks; M++ ) {
                for ( unsigned N = 0; N < numColBlocks; N++ ) {

                    const double entry =
                        (testFun[ M] - testAvg[ N]/sizeElem) *
                        (trialFun[N] - trialAvg[M]/sizeElem) * weight * detJ;

                    matrix( M, N ) -= entry;
                }
            }
        
        }
    }

    //! For compatibility with the base::asmb::computeStiffnessMatrix routine
    void tangentStiffness( const FieldTuple&  fieldTuple,
                           const Quadrature&  quadrature, 
                           base::MatrixD&     matrix ) const
    {
        return this -> operator()( fieldTuple, quadrature, matrix );
    }


};

#endif
