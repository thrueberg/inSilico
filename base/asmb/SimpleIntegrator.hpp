//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SimpleIntegrator.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_simpleintegrator_hpp
#define base_asmb_simpleintegrator_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
// boost includes
#include <boost/bind.hpp>
#include <boost/function.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{
        
        template<typename QUAD, typename RESULT, typename FIELDTUPLE>
        class SimpleIntegrator;

        //--------------------------------------------------------------------------
        // Convenience function to perform the integration
        template<typename FIELDTUPLEBINDER,
                 typename QUADRATURE, typename RESULT,
                 typename FIELDBINDER, typename KERNEL>
        void simplyIntegrate( const QUADRATURE&  quadrature,
                              RESULT&            result, 
                              const FIELDBINDER& fieldBinder,
                              const KERNEL&      kernelObj )
        {
            // define integrator
            typedef SimpleIntegrator<QUADRATURE,RESULT,
                                     typename FIELDTUPLEBINDER::Tuple>
                SimpleIntegrator;

            // bind kernel function (obj has operator() defined)
            typename SimpleIntegrator::Kernel
                kernel = boost::bind( kernelObj, _1, _2, _3, _4 );

            // create integrator
            SimpleIntegrator simpleInt( kernel, quadrature, result );

            // apply to all elements
            typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
            typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
            for ( ; iter != end; ++iter ) {
                simpleInt( FIELDTUPLEBINDER::makeTuple( *iter ) );
            }
        
            return;
        }

        
    } // namespace asmb    
} // namespace base

//------------------------------------------------------------------------------
/** Integrate without assembly.
 *  Perform integration over the given field, but only pass down a result
 *  container. Possible application is the computation of the mesh volume in
 *  which a floating point  number is passed on and the integration result is
 *  summed up.
 *  \tparam QUAD        Type of quadrature to use
 *  \tparam RESULT      Result storage
 *  \tparam FIELDTUPLE  Tuple of field elements
 */
template<typename QUAD, typename RESULT, typename FIELDTUPLE>
class base::asmb::SimpleIntegrator
    : public boost::function<void( const FIELDTUPLE& )>
{
public:
    //! @name Template parameter
    //@{
    typedef QUAD        Quadrature;
    typedef RESULT      Result;
    typedef FIELDTUPLE  FieldTuple;
    //@}

    //! Type of integration kernel
    typedef boost::function< void( const FieldTuple&,
                                   const typename Quadrature::VecDim&,
                                   const double,
                                   Result& )>  Kernel;

    //--------------------------------------------------------------------------
    //! Constructor setting all references
    SimpleIntegrator( Kernel&              kernel,
                      const Quadrature&    quadrature,
                      Result&              result )
        : kernel_(      kernel ),
          quadrature_(  quadrature ),
          result_(      result )
    { }

    //--------------------------------------------------------------------------
    //! Simply let quadrature integrate over the kernel function
    void operator()( const FieldTuple& fieldTuple )
    {
        // apply quadrature
        quadrature_.apply( kernel_, fieldTuple, result_ );
        
        return;
    }

private:
    Kernel&              kernel_;     //!< Function producing the result
    const Quadrature&    quadrature_; //!< Quadrature
    Result&              result_;     //!< Result container

};


#endif

