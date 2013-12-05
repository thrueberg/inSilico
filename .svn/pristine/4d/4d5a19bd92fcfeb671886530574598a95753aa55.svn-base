//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   gp/Writer.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_gp_writer_hpp
#define base_io_gp_writer_hpp

//------------------------------------------------------------------------------
// std   includes
#include <ostream>
#include <vector>
// base includes
#include <base/shape.hpp>
#include <base/types.hpp>
// base/io/gp includes
#include <base/io/gp/ElementFilter.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace gp{

            struct Writer;
        }
    }
}

//------------------------------------------------------------------------------
/** Writes the element coordinates to a stream digestable for GnuPlot.
 *
 */
struct base::io::gp::Writer
{

    template<typename EITER>
    static void apply( EITER begin, EITER end, std::ostream& out )
    {
        // Deduce type of element
        typedef typename
            base::TypeReduction<typename EITER::value_type>::Type Element;
    
        // @name Deduced attributes for writing 
        //@{
        static const base::Shape elementShape     = Element::shape;
        static const unsigned    nNodesPerElement = Element::numNodes;
        static const unsigned    coordDim         = Element::Node::dim;
        //@}

        // Filter for the element coordinates
        typedef base::io::gp::ElementFilter<elementShape,nNodesPerElement> EF;
        EF elemFilter;

        // Go through range of elements
        for ( EITER eIter = begin; eIter != end; ++eIter ) {

            // Go through the element filter
            typename EF::FilterIter fBegin = elemFilter.begin();
            typename EF::FilterIter fEnd   = elemFilter.end();
            for ( ; fBegin != fEnd; ++fBegin ) {

                const int index = *fBegin;

                // Positive entry: write node coordinates and ID
                if ( index >= 0 ) {

                    std::vector<double> x( coordDim );
                    (*eIter) -> nodePtr( index ) -> getX( x.begin() );
 

                    for ( unsigned d = 0; d < coordDim; d++ )
                        out << x[d] << " ";

                    out << (*eIter) -> nodePtr( index ) -> getID();
                }
                out << '\n';
            }

            out << '\n';
        }

        return;
    }

};
#endif
