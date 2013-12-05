//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ComputeSupport.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_computesupport_hpp
#define base_cut_computesupport_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <algorithm>
// base  includes
#include <base/shape.hpp>
#include <base/linearAlgebra.hpp>

#include <base/kernel/Measure.hpp>
#include <base/kernel/Mass.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename FIELDTUPLE, typename CUTQUAD> class ComputeSupport;

        
        template<typename FIELDTUPLEBINDER,
                 typename FIELDBINDER, typename CUTQUAD>
        void supportComputation( const FIELDBINDER& fieldBinder,
                                 const CUTQUAD& cutQuadrature, 
                                 std::vector<double>& supports,
                                 const bool inside = true )
        {
            ComputeSupport<typename FIELDTUPLEBINDER::Tuple,CUTQUAD>
                supportComputer( cutQuadrature, supports, inside );

            typename FIELDBINDER::FieldIterator fIter = fieldBinder.elementsBegin();
            typename FIELDBINDER::FieldIterator fEnd  = fieldBinder.elementsEnd();
            for ( ; fIter != fEnd; ++fIter )
                supportComputer( FIELDTUPLEBINDER::makeTuple( *fIter ) );
            
        }

        
        
    }
}

//------------------------------------------------------------------------------
template<typename FIELDTUPLE, typename CUTQUAD>
class base::cut::ComputeSupport
{
public:
    //! Template parameters
    //@{
    typedef FIELDTUPLE FieldTuple;
    typedef CUTQUAD    CutQuadrature;
    //@}

    typedef typename FieldTuple::TestElement FieldElement;
    
    //!

    ComputeSupport( const CutQuadrature& cutQuadrature, 
                    std::vector<double>& doFSupports,
                    const bool inside  = true )
        : cutQuadrature_( cutQuadrature ), doFSupports_( doFSupports ),
          inside_( true )
    {
        std::fill( doFSupports_.begin(), doFSupports_.end(), 0. );
    }

    void operator()( const FieldTuple& fieldTuple )
    {
        FieldElement* fElement = fieldTuple.testElementPtr();
        
        // add to supports
        typename FieldElement::DoFPtrConstIter doFIter = fElement -> doFsBegin();
        typename FieldElement::DoFPtrConstIter doFEnd  = fElement -> doFsEnd();
        
        
#ifdef MASS
        const std::size_t nDoFs = std::distance( doFIter, doFEnd );
        base::MatrixD massMat = base::MatrixD::Zero( nDoFs, nDoFs );
        base::kernel::Mass<FieldTuple,1> mass( 1.0 );
        cutQuadrature_.apply( mass, fieldTuple, massMat );
        
        for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ )
            doFSupports_[ (*doFIter) -> getID() ] += std::sqrt( massMat( d, d ) );
#else
        base::kernel::Measure<FieldTuple> measure;
        double area = 0.;
        cutQuadrature_.apply( measure, fieldTuple, area );

        for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ )
            doFSupports_[ (*doFIter) -> getID() ] += area;
#endif

        return;
    }

private:
    const CutQuadrature&     cutQuadrature_;
    std::vector<double>&     doFSupports_;
    const bool               inside_;
};


#endif
