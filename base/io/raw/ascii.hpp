//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ascii.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_raw_ascii_hpp
#define base_io_raw_ascii_hpp

//------------------------------------------------------------------------------
// std   includes
#include <ostream>
#include <istream>
#include <vector>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace raw{

            //------------------------------------------------------------------
            namespace detail_{

                template<typename DOFVALUE>
                struct WriteDoFValue
                {
                    static void apply( const DOFVALUE& value,
                                       std::ostream& out )
                    {
                        for ( int d1 = 0; d1 < value.rows(); d1 ++ ) {
                            for ( int d2 = 0; d2 < value.cols(); d2 ++ ) {
                                out << value(d1,d2) << " ";
                            }
                        }
                        // hack to avoid 2D vectors
                        if ( value.size() == 2 ) out << "0";
                    }
                };

                template<>
                struct WriteDoFValue<double>
                {
                    static void apply( const double& value,
                                       std::ostream& out )
                    {
                        out << value;
                    }

                };

            }
            
            //------------------------------------------------------------------
            /** Write nodal coordinates.
             *  \tparam NODE Type of node with coordinates
             */
            template<typename NODE>
            void writeNodeCoordinates( const NODE* np, std::ostream& out )
            {
                std::vector<double> x( NODE::dim );
                np -> getX( x.begin() );
                
                for ( unsigned d = 0; d < NODE::dim; d ++ )
                    out << x[d] << " ";
                for ( unsigned d = NODE::dim; d < 3; d ++ )
                    out << "0 ";
                out << "\n";
            }

            //------------------------------------------------------------------
            /** Read nodal coordinates.
             *  \tparam NODE Type of node
             */
            template<typename NODE>
            void readNodeCoordinates( NODE* np, std::istream& inp )
            {
                std::vector<double> coord;
                for ( unsigned d = 0; d < NODE::dim; d ++ ) {
                    double x;
                    inp >> x;
                    coord.push_back( x );
                }
                
                np -> setX( coord.begin() );
            }

            //------------------------------------------------------------------
            /** Write element's connectivity to stream
             *  \tparam ELEMENT      Type of element
             *  \tparam NUMOUTNODES  Number of nodes to be written (serendipity)
             */
            template<typename ELEMENT, unsigned NUMOUTNODES>
            void writeElementConnectivity( const ELEMENT* ep, std::ostream& out)
            {
                typename ELEMENT::NodePtrConstIter first = ep -> nodesBegin();
                typename ELEMENT::NodePtrConstIter  last = ep -> nodesEnd();

                for ( unsigned n = 0;
                      ( first != last ) and ( n < NUMOUTNODES); ++first, n++ ) {
                    //! Collect node ID
                    const std::size_t nodeID = (*first) -> getID();
                    out << nodeID << " ";
                }

                out << "\n";
            }

            //------------------------------------------------------------------
            /** Read the connectivity of an element and pass node pointers.
             *  \tparam MESH Type of mesh
             */
            template<typename MESH>
            void readElementConnectivity( const MESH& mesh,
                                          typename MESH::Element* ep,
                                          std::istream& inp )
            {
                typename MESH::Element::NodePtrIter elemNIter = ep -> nodesBegin();
                typename MESH::Element::NodePtrIter elemNEnd  = ep -> nodesEnd();

                // go through geometry nodes of element
                for ( ; elemNIter != elemNEnd; ++elemNIter ) {

                    // Read ID of vertex
                    std::size_t vertexID;
                    inp >> vertexID;

                    // assign node pointer to element
                    *elemNIter = mesh.nodePtr( vertexID );
                }
                return;
            }


            //------------------------------------------------------------------
            /** Write dof values to stream
             *  \tparam DOFVAL Type of value to be written to stream
             */
            template<typename DOFVAL>
            void writeDoFValues( const DOFVAL& dv, std::ostream & out )
            {
                detail_::WriteDoFValue<DOFVAL>::apply( dv, out );
                
                out << "\n";
            }
            
                        
        }
    }
}

#endif
