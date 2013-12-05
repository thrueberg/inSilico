//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   fe/Policies.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_fe_policies_hpp
#define base_fe_policies_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/shape.hpp>
//------------------------------------------------------------------------------
namespace base{
    namespace fe{

        //----------------------------------------------------------------------
        /** Define the positions in a hiearchically ordered array that
         *  correspond to a given n-face.
         *  This helper object defines the begin and end positions of the
         *  range and the stride which is the number of entites per n-face.
         *  An entity will be either an FE degree of freedom or a geometry node.
         *
         *  \note \a Note: The notation is here to refer to such an entity as
         *                 DoF, but in praxis this could also be a geometry node.
         *                 As both are stored in the same hierarchical ordering,
         *                 this helper can be used for real FE-DoFs and to nodes.
         *
         *  \tparam FELEMENT  Type of finite element which defines these values
         *  \tparam NFACE     Type of n-face which is of interest
         */
        template<typename FELEMENT, base::NFace NFACE>
        struct DoFArrayPositions;

        //! \cond SKIPDOX
        
        //! Vertices come first
        template<typename FELEMENT>
        struct DoFArrayPositions<FELEMENT,base::VERTEX>
        {
            static const int begin  = 0;
            static const int stride = FELEMENT::numDoFsPerVertex;
            static const int   end  = begin + FELEMENT::numVertexDoFs;
        };

        //! Edges follow after vertices
        template<typename FELEMENT>
        struct DoFArrayPositions<FELEMENT,base::EDGE>
        {
            static const int begin  = DoFArrayPositions<FELEMENT,base::VERTEX>::end;
            static const int stride = FELEMENT::numDoFsPerEdge;
            static const int   end  = begin + FELEMENT::numEdgeDoFs;
        };

        //! Faces follow after edges
        template<typename FELEMENT>
        struct DoFArrayPositions<FELEMENT,base::FACE>
        {
            static const int begin  = DoFArrayPositions<FELEMENT,base::EDGE>::end;
            static const int stride = FELEMENT::numDoFsPerFace;
            static const int   end  = begin + FELEMENT::numFaceDoFs;
        };

        //! Cells follow after faces
        template<typename FELEMENT>
        struct DoFArrayPositions<FELEMENT,base::CELL>
        {
            static const int begin  = DoFArrayPositions<FELEMENT,base::FACE>::end;
            static const int stride = FELEMENT::numDoFsPerCell;
            static const int   end  = begin + FELEMENT::numCellDoFs;
        };
        //! \endcond
        //======================================================================
        
        //----------------------------------------------------------------------
        /** For a given face number, collect the indices which belong to that
         *  face \em including the boundaries of that face.
         *  As an example, think of the task to get all geometry nodes (or
         *  degrees of freedom) which lie on the 3rd face of a HEX. This face,
         *  which is a QUAD, has four vertices and four edges belonging to it.
         *  Using
         *  \code
         *        FaceExtraction<FELEMENT,base::FACE>::apply( 3, ids );
         *  \endcode
         *  will fill the vector<unsigned> ids with the indices which belong
         *  (in this order!) to the vertices, edges and face of the 3rd face.
         *  Here, FELEMENT hides a finite element which in this example
         *  is a HEX.
         *
         *  \tparam FELEMENT  Type of finite element, e.g. Lagrangian
         *  \tparam NFACE     Type of n-face to extract from element
         */
        template<typename FELEMENT, base::NFace NFACE>
        struct FaceExtraction;

        //! Insert specific vertex DoF to faces dofs
        template<typename FELEMENT>
        struct FaceExtraction<FELEMENT,base::VERTEX>
        {
            static void apply( const unsigned vertexNumber,
                               std::vector<unsigned>& faceDoFIDs )
            {
                faceDoFIDs.push_back( vertexNumber );
            }
        };

        //! Copy the vertex DoFs and the DoFs on edge for given edge number
        template<typename FELEMENT>
        struct FaceExtraction<FELEMENT,base::EDGE>
        {
            static void apply( const unsigned edgeNumber,
                               std::vector<unsigned>& faceDoFIDs )
            {
                // get vertex dofs of this edge
                typedef 
                    base::mesh::ElementFaces<FELEMENT::shape, base::EDGE>
                    ElementEdges;
                    
                for ( unsigned v = 0; v < 2; v ++ ) {
                    const unsigned vertexNum = ElementEdges::index( edgeNumber, v );
                    faceDoFIDs.push_back( vertexNum );
                }

                // get dofs on edge itself
                typedef DoFArrayPositions<FELEMENT,base::EDGE> DAP;

                for ( int d = 0; d < DAP::stride; d++ )
                    faceDoFIDs.push_back( DAP::begin + edgeNumber*DAP::stride + d );
            }
        };

        /** Copy the vertex DoF IDs of the given face,
         *  find all edges of the face and copy their DoF IDs,
         *  and copy all DoF IDs of the face interior
         */
        template<typename FELEMENT>
        struct FaceExtraction<FELEMENT,base::FACE>
        {
            static void apply( const unsigned faceNumber,
                               std::vector<unsigned>& faceDoFIDs )
            {
                // get vertex dofs of this edge
                typedef 
                    base::mesh::ElementFaces<FELEMENT::shape, base::FACE>
                    ElementFaces;

                static const unsigned numVertices =
                    ElementFaces::ElementFaceTraits::numVertices;

                for ( unsigned v = 0; v < numVertices; v ++ ) {
                    const unsigned vertexNumber = ElementFaces::index( faceNumber, v );
                    faceDoFIDs.push_back( vertexNumber );
                }
                    
                // get dofs on of edges
                typedef
                    base::mesh::FaceEdges<FELEMENT::shape>
                    FaceEdges;

                typedef DoFArrayPositions<FELEMENT,base::EDGE> EdgeDAP;

                for ( unsigned e = 0; e < FaceEdges::numFaceEdges; e ++ ){
                    const int edgeNumber = FaceEdges::index( e, faceNumber );
                    const int edgeOrient = FaceEdges::sign(  e, faceNumber );
                    for ( int d = 0; d < EdgeDAP::stride; d++ ) {
                        // take care of orientation of edge
                        const unsigned ctr =
                            (edgeOrient > 0 ? d : EdgeDAP::stride - d - 1);
                        
                        faceDoFIDs.push_back( EdgeDAP::begin +
                                              edgeNumber * EdgeDAP::stride
                                              + ctr );
                    }

                }

                // get dofs on face itself
                typedef DoFArrayPositions<FELEMENT,base::FACE> FaceDAP;

                for ( int d = 0; d < FaceDAP::stride; d++ )
                    faceDoFIDs.push_back( FaceDAP::begin +
                                          faceNumber*FaceDAP::stride + d );
            }
        };

        //! Copy all element DoFs to faceDoFs for face=CELL
        template<typename FELEMENT>
        struct FaceExtraction<FELEMENT,base::CELL>
        {
            static void apply( const unsigned dummy,
                               std::vector<unsigned>& faceDoFIDs )
            {
                typedef DoFArrayPositions<FELEMENT,base::CELL> CellDAP;
                for ( unsigned i = 0; i < CellDAP::end; i ++ )
                    faceDoFIDs.push_back( i );
            }
        };

        
    }
}
#endif
