//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   LegacyWriter.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_vtk_legacywriter_hpp
#define base_io_vtk_legacywriter_hpp

//------------------------------------------------------------------------------
// base/mesh includes
#include <base/mesh/sampleStructured.hpp>
// base/io includes
#include <base/io/OStreamIterator.hpp>
// base/io/raw includes
#include <base/io/raw/ascii.hpp>
// base/io/vtk includes
#include <base/io/vtk/CellTraits.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace vtk{

            class LegacyWriter;

            namespace detail_{

                //--------------------------------------------------------------
                template<typename DATUM>
                struct IsScalar
                {
                    static const bool value = (DATUM::SizeAtCompileTime == 1);
                };

                template<>
                struct IsScalar<double>
                {
                    static const bool value = true;
                };

                //--------------------------------------------------------------
                template<typename DATUM>
                struct IsTensor
                {
                    static const bool value = (DATUM::SizeAtCompileTime == 9);
                };

                template<>
                struct IsTensor<double>
                {
                    static const bool value = false;
                };


            }
            
        }
    }
}

//------------------------------------------------------------------------------
/** Write output in an ASCII VTK Legacy file format.
 *  See www.vtk.org/VTK/img/file-formats.pdf for a description.
 *  
 */
class base::io::vtk::LegacyWriter
{
public:
    LegacyWriter( std::ostream & vtk,
                  const unsigned gridResolution = 1 )
        : vtk_( vtk ),
          gridResolution_( gridResolution ),
          pointDataOpen_( false ),
          cellDataOpen_(  false )
    {
        vtk_ << "# vtk DataFile Version 2.0\n"
             << "Generated by inSilico's base::io::vtk::LegacyWriter\n"
             << "ASCII\n";
        return;
    }

public:

    //! @name Geometry output
    //@{
    template<typename MESH>
    void writeUnstructuredGrid( const MESH& mesh );

    template<typename MESH>
    void writeStructuredGrid( const MESH& mesh );
    //@}

    //! @name Data field output
    //@{
    template<typename VALITER>
    void writePointData( VALITER first, VALITER last, const std::string & name )
    {
        this -> writeData_<VALITER>( first, last, name, false );
    }

    template<typename VALITER>
    void writeCellData( VALITER first, VALITER last, const std::string & name )
    {
        this -> writeData_<VALITER>( first, last, name, true );
    }
    //@}

private:
    template<typename VALITER>
    void writeData_( VALITER first, VALITER last, const std::string & name,
                     const bool isCellData );

private:
    std::ostream & vtk_; //!< Output stream
    const unsigned gridResolution_; //!< Possible resolution for grids
    bool pointDataOpen_; //!< State flag for open point data section
    bool cellDataOpen_;  //!< State flag for open cell data section
};

//------------------------------------------------------------------------------
/** Write a given mesh as unstructured grid.
 *  \param[in] mesh  Mesh to be written.
 */
template<typename MESH>
void base::io::vtk::LegacyWriter::writeUnstructuredGrid( const MESH& mesh )
{
    // node range iterators and distance
    typename MESH::NodePtrConstIter node    = mesh.nodesBegin();
    typename MESH::NodePtrConstIter nodeEnd = mesh.nodesEnd();
    const std::size_t numNodes = std::distance( node, nodeEnd );

    // header for type of grid and point coordinates
    vtk_ << "DATASET UNSTRUCTURED_GRID" << "\n"
         << "POINTS " << numNodes << " float" << "\n";

    // write the points by querying the nodes
    typedef typename MESH::Node Node;
    std::copy( node, nodeEnd,
               base::io::OStreamIterator<Node*>(
                   vtk_,
                   &base::io::raw::template writeNodeCoordinates<Node> ) );

    // element traits
    typedef base::io::vtk::CellType<MESH::Element::shape,
                                    MESH::Element::numNodes> CT;
    typedef base::io::vtk::CellNumOutputNodes<MESH::Element::shape,
                                              MESH::Element::numNodes> CNON;

    // element range iterators and distance
    typename MESH::ElementPtrConstIter element    = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter elementEnd = mesh.elementsEnd();
    const std::size_t numElements = std::distance( element, elementEnd );

    // expected size of array
    const std::size_t sizeOfArray = numElements * (+CNON::value + 1);

    // write elements
    vtk_ << "CELLS " << numElements << " " << sizeOfArray << "\n";
    for ( ; element != elementEnd; ++element ) {
        vtk_ << +CNON::value << " ";
        base::io::raw::writeElementConnectivity<typename MESH::Element,
                                                CNON::value>( *element, vtk_ );
    }

    // write cell types
    vtk_ << "CELL_TYPES " << numElements << "\n";
    element = mesh.elementsBegin();
    for ( ; element != elementEnd; ++ element ) vtk_ << CT::value << "\n";

    return;
}

//------------------------------------------------------------------------------
/** Write the given mesh as a structured grid.
 *  A delegate call to base::mesh::sampleGridGeometry() yields the grid
 *  coordinates of the given mesh. These are then written as 3D coordinates
 *  to the stream.
 *  If a higher resolution is wanted, this has to be set in the constructor of
 *  this class.
 *  \param[in] grid  Grid to be sampled and written
 */
template<typename GRID>
void base::io::vtk::LegacyWriter::writeStructuredGrid( const GRID& grid )
{
    // get grid dimensions
    const typename GRID::MultiIndexType gridSizes = grid.gridSizes();

    // write header with grid dimensions
    // (here the values are the number of points per direction)
    vtk_ << "DATASET STRUCTURED_GRID" << "\n"
         << "DIMENSIONS ";
    for ( unsigned d = 0; d < GRID::dim; d++ )
        vtk_ << gridSizes[d] + 1 << " ";
    for ( unsigned d = GRID::dim; d < 3; d++ )
        vtk_ << "1 "; // 1 point in the orthogonal directions
    vtk_ << "\n";

    // get the coordinates via geometry sampling
    typedef typename GRID::Node::VecDim VecDim;
    std::vector<VecDim> coordinates;
    base::mesh::sampleGridGeometry( grid,
                                    std::back_inserter( coordinates ),
                                    gridResolution_ );
    

    // Write points
    vtk_ << "POINTS " << coordinates.size() << " float" << "\n";
    for ( typename std::vector<VecDim>::iterator iter = coordinates.begin();
          iter != coordinates.end(); ++iter ) {

        // write as 3D coordinates
        for ( unsigned d = 0; d < GRID::Node::dim; d++ )
            vtk_ << (*iter)[d] << " ";
        for ( unsigned d = GRID::Node::dim; d < 3; d++ )
            vtk_ << "0 ";
        vtk_ << "\n";
    }

    return;
}

//------------------------------------------------------------------------------
/** Write data field (point or cell datum) to file.
 *  The data field is represented by a range of iterators and will be given a
 *  name for the VTK viewer. An additional flag decides if the datum is
 *  written as a point or a cell datum.
 *  \param[in] first, last   Range of iterators
 *  \param[in] name          Human-readable descriptor of the datum
 *  \param[in] isCellData    Flag which evaluates to true for cell data
 */
template<typename VALITER>
void base::io::vtk::LegacyWriter::writeData_( VALITER first,
                                              VALITER last, 
                                              const std::string & name,
                                              const bool isCellData )
{
    // Deduce type of point datum
    typedef typename VALITER::value_type DoFValue;
    
    // Compute the number of dofs
    const std::size_t numDoFs = std::distance( first, last );

    // Check if datum is scalar field
    const bool isScalar = detail_::IsScalar<DoFValue>::value;
    const bool isTensor = detail_::IsTensor<DoFValue>::value;

    if ( isCellData ) {
    
        if ( not cellDataOpen_ ) {
            vtk_ << "CELL_DATA " << numDoFs << "\n";
            cellDataOpen_ = true;
        }
    }
    else{
        
        if ( not pointDataOpen_ ) {
            vtk_ << "POINT_DATA " << numDoFs << "\n";
            pointDataOpen_ = true;
        }
    }

    // Write data type identifier
    vtk_ << ( isScalar ? "SCALARS" : (isTensor ? "TENSORS" : "VECTORS" ) )
         << " " << name << " float " << "\n";

    // Strange VTK policy
    if ( isScalar ) vtk_ << "LOOKUP_TABLE default \n";

    // write the data
    std::copy( first, last, 
               base::io::OStreamIterator<DoFValue>(
                   vtk_, &base::io::raw::template writeDoFValues<DoFValue> ) );

}


#endif
