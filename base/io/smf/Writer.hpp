//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smf/Writer.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_smf_writer_hpp
#define base_io_smf_writer_hpp

//------------------------------------------------------------------------------
// std   includes
#include <ostream>
#include <vector>
#include <algorithm>
// boost includes
#include <boost/bind.hpp>
#include <boost/ref.hpp>
// base includes
#include <base/shape.hpp>
// base/io includes
#include <base/io/OStreamIterator.hpp>
// raw includes
#include <base/io/raw/ascii.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace smf{
            
            template<typename MESH> class Writer;

            //! Convenience function to write a mesh in SMF format
            template<typename MESH>
            void writeMesh( const MESH& mesh, std::ostream& smf )
            {
                Writer<MESH> writer;
                writer( mesh, smf );
            }
            
        }
    }
}

//------------------------------------------------------------------------------
/** Writes the mesh data (geometry) in the SMF file format
 *
 *  For a description of this format, see base::io::smf.
 *
 *  \tparam MESH Type of mesh to be written
 */
template<typename MESH>
class base::io::smf::Writer
{
public:
    //! Template parameter: type of mesh
    typedef MESH Mesh;

    //! @name Deduced attributes for writing 
    //@{
    static const base::Shape elementShape     = Mesh::Element::shape;
    static const unsigned    nNodesPerElement = Mesh::Element::numNodes;
    static const unsigned    coordDim         = Mesh::Node::dim;
    //@}

    /** Overloaded function call operator to write the mesh geometry
     *  \param[in] mesh  The mesh to be written out
     *  \param[in] smf   Output stream for the writing
     */
    void operator()( const Mesh & mesh, std::ostream & smf ) const
    {
        // write header
        this -> writeHeader_( smf );

        // write number of nodes and elements
        {
            const std::size_t nNodes    = std::distance( mesh.nodesBegin(),
                                                         mesh.nodesEnd() );
            const std::size_t nElements = std::distance( mesh.elementsBegin(),
                                                         mesh.elementsEnd() );
            smf << nNodes << "  " << nElements << "\n";
        }
        
        // write nodal coordinates
        this -> writeNodes_( smf, mesh );

        // write elements' connectivities
        this -> writeElements_( smf, mesh );
    }

private:
    //--------------------------------------------------------------------------
    /** Write the header with the element type description
     * \param[in] smf output stream
     */
    void writeHeader_( std::ostream & smf ) const
    {
        const char headerChar = '!';
        
        smf << headerChar << " elementShape "
            << base::ShapeName<elementShape>::apply() << "\n"
            << headerChar << " elementNumPoints "
            << nNodesPerElement << "\n";
    }

    //--------------------------------------------------------------------------
    /** Write the nodal coordinates (always with trailing zeros as if 3D).
     *  \param[in] smf  Output stream
     *  \param[in] mesh Reference to mesh
     */
    void writeNodes_( std::ostream & smf, const Mesh & mesh ) const
    {
        // Use custom ostream iterator for writing the coordinates
        typedef typename Mesh::Node Node;
        
        std::copy( mesh.nodesBegin(), mesh.nodesEnd(), 
                   base::io::OStreamIterator<Node*>(
                       smf,
                       &base::io::raw::template writeNodeCoordinates<Node> ) );
        
    }

    //--------------------------------------------------------------------------
    /** Write the elements' connectivities in terms of node indices
     *  \param[in] smf  Output stream
     *  \param[in] mesh Reference to mesh
     */
    void writeElements_( std::ostream & smf, const Mesh & mesh ) const
    {
        // Use custom ostream iterator for writing the connectivity
        typedef typename Mesh::Element Element;
        
        std::copy( mesh.elementsBegin(), mesh.elementsEnd(),
                   base::io::OStreamIterator<Element*>(
                       smf,
                       &base::io::raw::template writeElementConnectivity<Element,
                                                                         Element::numNodes> ) );
    }

};
#endif
