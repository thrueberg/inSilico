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

            //! Helper to define the basis field
            template<typename FEBASIS, unsigned DOFSIZE, unsigned NHIST>
            struct ScaledFieldTraits
            {
                typedef base::dof::DegreeOfFreedom<DOFSIZE,NHIST>    DegreeOfFreedom;
                typedef typename FEBASIS::FEFun                      OriginalFieldFun;
                typedef base::cut::ScaledBasis<OriginalFieldFun>     FieldFun;
                typedef base::dof::Element<DegreeOfFreedom,FieldFun> Element;
                typedef base::dof::Field<Element>                    Type;
            };
        }

        template<typename FEBASIS, unsigned DOFSIZE, unsigned NHIST = 0>
        struct ScaledField;
    }
}


//------------------------------------------------------------------------------
/** Overload field with scaled shape functions.
 *  A possible solution to the cut-cell conditioning problem is to re-scale the
 *  shape functions inverse proportionally to their active support size.
 *  Effectively, this boils down to a diagonal pre-conditioning. Moreover, DoFs
 *  with a support size below a given threshold are deactivated. 
 *  
 */
template<typename FEBASIS, unsigned DOFSIZE, unsigned NHIST>
struct base::cut::ScaledField
    : public base::cut::detail_::ScaledFieldTraits<FEBASIS,DOFSIZE,NHIST>::Type
{
public:
    //! Type of field to inherit from
    typedef typename
    base::cut::detail_::ScaledFieldTraits<FEBASIS,DOFSIZE,NHIST>::Type Basis;

    //! Activate all DoFs and set all scaling values to one.
    void activateAll()
    {
        typename Basis::DoFPtrIter doFIter = Basis::doFsBegin();
        typename Basis::DoFPtrIter doFEnd  = Basis::doFsEnd();
        for ( ; doFIter != doFEnd; ++doFIter ) {
            (*doFIter) -> activateAll();
            (*doFIter) -> clearConstraints();
        }

        typename Basis::ElementPtrIter elemIter = Basis::elementsBegin();
        typename Basis::ElementPtrIter elemEnd  = Basis::elementsEnd();
        for ( ; elemIter != elemEnd; ++elemIter ) {
            ( (*elemIter) -> fEFun() ).reset();
        }
    }

    /** Set the basis scaling factors to inverse support size, scale the values
     *  stored in the DoF objects, and tag DoFs with support size below a given
     *  threshold to be inactive.
     *  
     */
    void scaleAndTagBasis( const std::vector<double>& supportAreas,
                           const double lowerThreshold )
    {
        // markers in order to avoid duplicate scaling
        std::vector<bool> mark( std::distance( Basis::doFsBegin(),
                                               Basis::doFsEnd() ), false );

        // go through all elements of the field
        typename Basis::ElementPtrIter elemIter = Basis::elementsBegin();
        typename Basis::ElementPtrIter elemEnd  = Basis::elementsEnd();
        for ( ; elemIter != elemEnd; ++elemIter ) {

            // go through all DoFs of each element
            typename Basis::Element::DoFPtrIter doFIter = (*elemIter) -> doFsBegin();
            typename Basis::Element::DoFPtrIter doFEnd  = (*elemIter) -> doFsEnd();
            for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ ) {
                        
                // obtain DoF ID and the size of the support
                const std::size_t doFID = (*doFIter) -> getID();
                const double supportArea = supportAreas[ doFID ];

                // set scaling factor in shape function basis
                if ( supportArea >= lowerThreshold ) {
                    ( (*elemIter) -> fEFun() ).setScalar( d, 1./supportArea );
                }

                // if DoF has not been treated so far
                if ( not mark[ doFID ] ) {

                    if ( supportArea < lowerThreshold ) {
                        // deactivate if support is below threshold
                        (*doFIter) -> deactivateAll();
                    }
                    else {
                        // scale the values stored in the DoF
                        (*doFIter) -> scaleAllValues( supportArea );
                    }

                    // avoid duplication
                    mark[doFID] = true;
                }

                        
            } // all DoFs per element
        } // all elements

        return;
    }

    /** Reset the scaling to one.
     *  For post-processing and in dynamic computations it is convenient
     *  to reverse the scaling of the basis functions and DoF values.
     */
    void unscaleValues( const std::vector<double>& supportAreas,
                        const double lowerThreshold )
    {
        // rescale all active dof values
        typename Basis::DoFPtrIter doFIter = Basis::doFsBegin();
        typename Basis::DoFPtrIter doFEnd  = Basis::doFsEnd();
        for ( ; doFIter != doFEnd; ++doFIter ) {

            const std::size_t doFID = (*doFIter) -> getID();
            const double supportArea = supportAreas[ doFID ];

            if ( supportArea > lowerThreshold ) {
                (*doFIter) -> scaleAllValues( 1./supportArea );
            }
        }

        // reset the scalar factor of the basis functions
        typename Basis::ElementPtrIter elemIter = Basis::elementsBegin();
        typename Basis::ElementPtrIter elemEnd  = Basis::elementsEnd();
        for ( ; elemIter != elemEnd; ++elemIter ) {
            ( (*elemIter) -> fEFun() ).reset();
        }
    }

            
    /**  Only deactivate the DoFs with a support that is below the given
     *   threshold.
     */
    void tagBasis( const std::vector<double>& supportAreas,
                   const double lowerThreshold )
    {
        // go through all field elements
        typename Basis::ElementPtrIter elemIter = Basis::elementsBegin();
        typename Basis::ElementPtrIter elemEnd  = Basis::elementsEnd();
        for ( ; elemIter != elemEnd; ++elemIter ) {
            // go through all DoFs of each element
            typename Basis::Element::DoFPtrIter doFIter = (*elemIter) -> doFsBegin();
            typename Basis::Element::DoFPtrIter doFEnd  = (*elemIter) -> doFsEnd();
            for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ ) {

                // get size of DoF's support
                const double supportArea =
                    supportAreas[ (*doFIter) -> getID() ];
                        
                // activate if support is too small
                if ( supportArea < lowerThreshold ) {
                    (*doFIter) -> deactivateAll();
                }

            }
        }
    }
};
        

#endif
