//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   cut/ScaledField.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_scaledfield_hpp
#define base_cut_scaledfield_hpp

//------------------------------------------------------------------------------
#include <base/shape.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Element.hpp>
#include <base/dof/Field.hpp>
#include <base/cut/ScaledBasis.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        
        namespace detail_{
            
            template<typename FEBASIS, unsigned DOFSIZE, unsigned NHIST>
            struct FieldTraits
            {
                typedef base::dof::DegreeOfFreedom<DOFSIZE,NHIST>    DegreeOfFreedom;
                typedef typename FEBASIS::FEFun                      OriginalFieldFun;
                typedef base::cut::ScaledBasis<OriginalFieldFun>     FieldFun;
                typedef base::dof::Element<DegreeOfFreedom,FieldFun> Element;
                typedef base::dof::Field<Element>                    Type;
            };
        }

        //----------------------------------------------------------------------
        template<typename FEBASIS, unsigned DOFSIZE, unsigned NHIST = 0>
        struct ScaledField
            : public detail_::FieldTraits<FEBASIS,DOFSIZE,NHIST>::Type
        {
        public:
            typedef typename detail_::FieldTraits<FEBASIS,DOFSIZE,NHIST>::Type Basis;

            //
            void scaleAndTagBasis( const std::vector<double>& supportAreas,
                                   const double lowerThreshold )
            {
                typename Basis::ElementPtrIter elemIter = Basis::elementsBegin();
                typename Basis::ElementPtrIter elemEnd  = Basis::elementsEnd();
                for ( ; elemIter != elemEnd; ++elemIter ) {
                    typename Basis::Element::DoFPtrIter doFIter = (*elemIter) -> doFsBegin();
                    typename Basis::Element::DoFPtrIter doFEnd  = (*elemIter) -> doFsEnd();
                    for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ ) {
                        //
                        const double supportArea =
                            supportAreas[ (*doFIter) -> getID() ];
                        
                        if ( supportArea < lowerThreshold ) {
                            (*doFIter) -> deactivateAll();
                        }
                        else {
                            const double factor = 1. / supportArea;
                            ( (*elemIter) -> fEFun() ).setScalar( d, factor );
                        }
                        
                    }
                }
            }

            //
            void tagBasis( const std::vector<double>& supportAreas,
                           const double lowerThreshold )
            {
                typename Basis::ElementPtrIter elemIter = Basis::elementsBegin();
                typename Basis::ElementPtrIter elemEnd  = Basis::elementsEnd();
                for ( ; elemIter != elemEnd; ++elemIter ) {
                    typename Basis::Element::DoFPtrIter doFIter = (*elemIter) -> doFsBegin();
                    typename Basis::Element::DoFPtrIter doFEnd  = (*elemIter) -> doFsEnd();
                    for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ ) {
                        //
                        const double supportArea =
                            supportAreas[ (*doFIter) -> getID() ];
                        
                        if ( supportArea < lowerThreshold ) {
                            (*doFIter) -> deactivateAll();
                        }
                        else {
                            const double factor = 1.;
                            ( (*elemIter) -> fEFun() ).setScalar( d, factor );
                        }

                    }
                }
            }

            
            
        };
        
    }
}

#endif
