//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   sampleStructured.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_samplestructured_hpp
#define base_samplestructured_hpp

//------------------------------------------------------------------------------
// std includes
#include <bitset>
// boost includes
#include <boost/function.hpp>
#include <boost/bind.hpp>
// base includes
#include <base/MultiIndex.hpp>
// base/mesh includes
#include <base/mesh/Structured.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename FIELD, typename SAMPLEOP, typename OUTITER>
        void sampleStructuredField( const FIELD&   field,
                                    SAMPLEOP&      sampling,
                                    OUTITER        outputIter,
                                    const typename
                                    base::MultiIndex<FIELD::Element::dim>::Type&
                                    gridSizes, 
                                    const unsigned resolution = 1,
                                    const bool     discont = false );

        //! Special case: Field=Grid and sampling of geometry
        template<typename GRID, typename OUTITER>
        void sampleGridGeometry( const GRID& grid,
                                 OUTITER        outputIter,
                                 const unsigned resolution = 1,
                                 const bool     discont = false )
        {
            typedef typename GRID::Element Element;
            
            // Function object for evaluating the element geometry
            typedef typename base::GeomTraits<Element>::LocalVecDim LVD;
            typedef boost::function< typename GRID::Node::VecDim(
                const Element*, const LVD&) > GeomEval;
        
            GeomEval geometry = boost::bind( base::Geometry<Element>(), _1, _2 );


            return sampleStructuredField( grid, geometry, outputIter,
                                          grid.gridSizes(), resolution, discont );
        }
        

    }
}

//------------------------------------------------------------------------------
/** A story to be told.
 */
template<typename FIELD, typename SAMPLEOP, typename OUTITER>
void base::mesh::sampleStructuredField( const FIELD&   field,
                                        SAMPLEOP&      sampling,
                                        OUTITER        outputIter,
                                        const typename
                                        base::MultiIndex<FIELD::Element::dim>::Type&
                                        gridSizes, 
                                        const unsigned resolution,
                                        const bool     discont  )
{
    // Deduce element type
    typedef typename FIELD::Element Element;

    // Local space dimension
    static const unsigned dim = Element::dim;

    // Multi-index structures
    typedef base::MultiIndex<dim>          MultiIndex;
    typedef typename MultiIndex::Type      MultiIndexType;

    // Vector of local coordinates
    typedef typename base::GeomTraits<Element>::LocalVecDim  LocalVecDim;

    // number of sampling points per direction
    const unsigned numSamplesPerDir = ( resolution == 0 ? 1 : resolution );

    // Samples multi-index
    MultiIndexType samples = MultiIndexType::Constant( numSamplesPerDir );

    // Subsample sizes
    const LocalVecDim subSizes =  1. / samples.template cast<double>();

    // increase the local sampling to capture right boundary
    samples += 1;

    // For mid-point sampling, shift the values
    const LocalVecDim shift = ( resolution == 0 ?
                                LocalVecDim::Constant( 0.5 ) :
                                LocalVecDim::Constant( 0.  ) );

    // Sub-sampled grid size
    const MultiIndexType gridSamples = samples * gridSizes;

    // Total number of sampling points
    const std::size_t numTotal = MultiIndex::length( gridSamples );

    // go through all elements and sampling points
    for ( std::size_t es = 0; es < numTotal; es ++ ) {

        // Wrap global counter back to global multi-index
        const MultiIndexType gM = MultiIndex::wrap( es, gridSamples );

        // Compute element counter by division
        const MultiIndexType eM = gM / samples;

        // Compute sampling counter by subtraction
        const MultiIndexType sM = gM - (samples * eM );

        //------------------------------------------------------------------
        // Check if point is to be sampled
        bool sampleValue = false;
            
        // In case of a discontinous field, always sample
        if ( discont or ( resolution == 0 ) ) sampleValue = true;

        // If this is not the case check if the point is
        // an interior sampling point or on the right edge
        // of a right-most element
        if ( not sampleValue ) {

            // Flags for last element and sample points per direction
            std::bitset<dim> lastElement;
            std::bitset<dim> lastSample;
            for ( unsigned d = 0; d < dim; d ++ ) {
                if ( eM[d] == gridSizes[d]-1 ) lastElement.set(d);
                if ( sM[d] == samples[d]  -1 ) lastSample.set(d);
            }

            // for interior sampling, always perform
            if ( lastSample.none() ) sampleValue = true;
            else {
                // In the very last (corner) element, evaluate all points
                if ( lastElement.count() == dim ) sampleValue = true;

                // If the flags are identical, evaluate
                if ( lastElement == lastSample ) sampleValue = true;
            }
        }

        // Sample
        if ( sampleValue ) {

            // Generate a local coordinate from multi-index
            LocalVecDim xi;
            for ( unsigned d = 0; d < dim; d++ )
                xi[d] =  static_cast<double>( sM[d] ) * subSizes[d];
            xi += shift;

            // linear index
            const std::size_t linearIndex = MultiIndex::unwrap( eM, gridSizes );
            
            // Get access to an element
            const Element* const ep = field.elementPtr( linearIndex );

            // sample operation
            *outputIter++ = sampling( ep, xi );
        }
           
    }// end loop over sampling points

    return;
}

#endif
