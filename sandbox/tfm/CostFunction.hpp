//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   CostFunction.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef tfm_costfunction_hpp
#define tfm_costfunction_hpp

#include <base/shape.hpp>
#include <base/linearAlgebra.hpp>
#include <base/Quadrature.hpp>
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace tfm{

    template<typename MESH, typename FIELD, typename TRACTION>
    class CostFunction;
    
}

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD, typename TRACTION>
class tfm::CostFunction
{
public:
    typedef typename base::Vector<FIELD::DegreeOfFreedom::size>::Type VecDoF;
    typedef typename base::Vector<MESH::Node::dim>::Type              VecDim;
    typedef boost::function<bool( const VecDim& )> FilterFun;

    
    CostFunction( MESH& mesh,
                  FIELD& currentU,
                  FIELD& targetU,
                  TRACTION& tractionHandler,              
                  const double delta,
                  FilterFun filterFun = boost::bind( tfm::AlwaysPositive<MESH>(), _1 )  )
        : mesh_(            mesh ),
          currentU_(        currentU ),
          targetU_(         targetU ),
          tractionHandler_( tractionHandler ),
          delta_(           delta ),
          filterFun_( filterFun )        
    { }

    double evaluate() const
    {
        // Quadratures
        static const unsigned kernelDegree = 2 * FIELD::Element::FEFun::degree;
        static const base::Shape shape     = FIELD::Element::shape;
        typedef base::Quadrature<       kernelDegree,shape> Quadrature;
        typedef base::SurfaceQuadrature<kernelDegree,shape> SurfaceQuadrature;

        //----------------------------------------------------------------------
        // compute
        double result = 0.;

        // Domain part
        const std::size_t numElements = std::distance( mesh_.elementsBegin(),
                                                       mesh_.elementsEnd() );
        Quadrature quadrature;

        for ( std::size_t e = 0; e < numElements; e++ ) {

            typename MESH::Element*  geomEp = mesh_.elementPtr( e );
            typename FIELD::Element* cuPtr  = currentU_.elementPtr( e );
            typename FIELD::Element* tuPtr  = targetU_.elementPtr(  e );

            typename Quadrature::Iter qIter = quadrature.begin();
            typename Quadrature::Iter qEnd  = quadrature.end();
            for ( ; qIter != qEnd; ++qIter ) {

                typename MESH::Node::VecDim x =
                    base::Geometry<typename MESH::Element>()( geomEp, qIter -> second );


                const VecDoF diffU =
                    ( filterFun_( x ) ? 
                      base::post::evaluateField( geomEp, cuPtr, qIter -> second ) -
                      base::post::evaluateField( geomEp, tuPtr, qIter -> second )
                      : base::constantVector<FIELD::DegreeOfFreedom::size>( 0. ) );

                result += base::dotProduct( diffU, diffU ) * (qIter -> first);
            }
            
        }

        // Regularisation part
        const std::size_t numElements2 = std::distance(
            tractionHandler_.getMesh().elementsBegin(),
            tractionHandler_.getMesh().elementsEnd() );
        SurfaceQuadrature surfaceQuadrature;
        
        for ( std::size_t e = 0; e < numElements2; e++ ) {

            typename TRACTION::SurfaceMesh::Element*  geomEp =
                tractionHandler_.getMesh().elementPtr( e );

            typename SurfaceQuadrature::Iter qIter = surfaceQuadrature.begin();
            typename SurfaceQuadrature::Iter qEnd  = surfaceQuadrature.end();
            for ( ; qIter != qEnd; ++qIter ) {

                const VecDoF t = tractionHandler_.apply( geomEp, qIter -> second );

                result += delta_ * 
                    base::dotProduct( t, t ) * (qIter -> first);
            }
            
        }
        

        return 0.5 * result;
    }
    
private:
    MESH&     mesh_;
    FIELD&    currentU_;
    FIELD&    targetU_;
    TRACTION& tractionHandler_;

    const double delta_;

    FilterFun filterFun_;
 
};

#endif
