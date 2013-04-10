//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   sgf/Writer.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_sgf_writer_hpp
#define base_io_sgf_writer_hpp

//------------------------------------------------------------------------------
// std   includes
#include <ostream>
// base includes
#include <base/geometry.hpp>
// base/mesh includes
#include <base/mesh/sampleStructured.hpp>
// base/io includes
#include <base/io/OStreamIterator.hpp>


//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace sgf{
            
            template<typename GRID> class Writer;

            namespace detail_{

                //! Helper to write the vectors as 3D
                template<typename VEC>
                void writeVec( const VEC & vec, std::ostream & out )
                {
                    for ( int d = 0; d < vec.size(); d++ )
                        out << vec[d] << " ";
                    for ( int d = vec.size(); d < 3; d++ )
                        out << "0 ";
                    out << '\n';
                }
            }

            
        }
    }
}

//------------------------------------------------------------------------------
/** Writes the grid data (geometry) in the SGF file format
 *
 *  For a description of this format, see base::io::sgf.
 *
 *  \tparam GRID Type of grid to be written
 */
template<typename GRID>
class base::io::sgf::Writer
{
public:

    /** Overloaded function call operator to write the grid geometry
     *  \param[in] grid       The grid to be written out
     *  \param[in] sgf        Output stream for the writing
     *  \param[in] resolution Optional resolution of the geometry sampling
     */
    void operator()( const GRID   & grid,
                     std::ostream & sgf,
                     const unsigned resolution = 1 ) const
    {
        const typename GRID::MultiIndexType gridSizes = grid.gridSizes();

        // write header
        for ( unsigned d = 0; d < GRID::dim; d++ )
            sgf << resolution * gridSizes[d] << " ";
        for ( unsigned d = GRID::dim; d < 3; d++ ) sgf << "0 ";
        sgf << '\n';

        // write nodal coordinates
        {
            typedef typename
                base::GeomTraits<typename GRID::Element>::GlobalVecDim GVD;

            // sample structured grid with sub-sampling
            base::mesh::sampleGridGeometry(
                grid,
                base::io::OStreamIterator< GVD >( sgf,
                                                  &detail_::writeVec<GVD> ),
                resolution );
        }

        return;
    }

};
#endif
