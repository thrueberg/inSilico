//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   cut/tagBasis.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_cut_tagbasis_hpp
#define base_cut_tagbasis_hpp

//------------------------------------------------------------------------------
#include <vector>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename FIELD>
        void tagBasis( FIELD& field, 
                       const std::vector<double>& supportAreas,
                       const double lowerThreshold );
        
    }
}

//------------------------------------------------------------------------------
/**  Deactivate the DoFs with a support that is below the given threshold.
 *   \tparam FIELD Type of field to treat
 *   \param[in,out] field          Field which is being tagged
 *   \param[in]     supportSizes   The size of every DoF's support
 *   \param[in]     lowerThreshold User-defined threshold for the tagging
 */
template<typename FIELD>
void base::cut::tagBasis( FIELD& field, 
                          const std::vector<double>& supportSizes,
                          const double lowerThreshold )
{
    // go through all field elements
    typename FIELD::ElementPtrIter elemIter = field.elementsBegin();
    typename FIELD::ElementPtrIter elemEnd  = field.elementsEnd();
    for ( ; elemIter != elemEnd; ++elemIter ) {
        
        // go through all DoFs of each element
        typename FIELD::Element::DoFPtrIter doFIter = (*elemIter) -> doFsBegin();
        typename FIELD::Element::DoFPtrIter doFEnd  = (*elemIter) -> doFsEnd();
        for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ ) {

            // get size of DoF's support
            const double supportSize = supportSizes[ (*doFIter) -> getID() ];
                        
            // activate if support is too small
            if ( supportSize < lowerThreshold ) (*doFIter) -> deactivateAll();
        }
    }
}

#endif
