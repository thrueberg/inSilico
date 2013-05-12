//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   marching.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_marching_hpp
#define base_cut_marching_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <utility>
#include <map>
#include <bitset>
// boost includes
#include <boost/array.hpp>
// base/cut includes
#include <base/cut/HyperCube.hpp>
#include <base/cut/DecomposeHyperCube.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        
        template<unsigned DIM>
        void marching( const boost::array<double,
                                          base::cut::HyperCube<DIM>::numVertices>&
                       signedDistances,
                       std::vector<typename base::Vector<DIM,double>::Type>& nodes,
                       std::vector<typename Simplex<DIM-1,unsigned>::Type>& surface,
                       std::vector<typename Simplex<DIM  ,unsigned>::Type>& volumeIn,
                       std::vector<typename Simplex<DIM  ,unsigned>::Type>& volumeOut );

        //----------------------------------------------------------------------
        //! Define a simplex to be a fixed-size array
        template<unsigned DIM, typename VAL=unsigned>
        struct Simplex
        {
            typedef boost::array<VAL,DIM+1> Type;
        };

        //----------------------------------------------------------------------
        /** Edge with comparison operator, allows unique identification of
         *  vertices and cut points.
         */
        class Edge
        {
        public:
            //! Construct from one vertex
            Edge( const unsigned v1 )
                : vertices_( std::make_pair( v1, base::invalidInt ) ) { }

            //! Construct from the line between two vertices
            Edge( const unsigned v1, const unsigned v2 )
                : vertices_( std::make_pair( (v1 < v2? v1 : v2),
                                             (v1 < v2? v2 : v1) ) ) { }

            //! Use build-in comparison for unique sorting
            bool operator<( const Edge& other ) const
            {
                return vertices_ < other.vertices_;
            }

        private:
            const std::pair<unsigned,unsigned> vertices_;
        };

        //----------------------------------------------------------------------
        //! @name Dimension-dependent algorithms
        //@{
        void bisect( const Simplex<1,unsigned>::Type& indexSimplex,
                     const Simplex<1,double>::Type&   distances,
                     std::vector<base::Vector<1,double>::Type>& nodes,
                     std::map<Edge,unsigned>& uniqueNodes,
                     std::vector<Simplex<0,unsigned>::Type>& surface,
                     std::vector<Simplex<1,unsigned>::Type>& volumeIn, 
                     std::vector<Simplex<1,unsigned>::Type>& volumeOut );

        void marchingTri( const Simplex<2,unsigned>::Type& indexSimplex,
                          const Simplex<2,double>::Type&   distances,
                          std::vector<base::Vector<2,double>::Type>& nodes,
                          std::map<Edge,unsigned>& uniqueNodes,
                          std::vector<Simplex<1,unsigned>::Type>& surface,
                          std::vector<Simplex<2,unsigned>::Type>& volumeIn, 
                          std::vector<Simplex<2,unsigned>::Type>& volumeOut );

        void marchingTet( const Simplex<3,unsigned>::Type& indexSimplex,
                          const Simplex<3,double>::Type&   distances,
                          std::vector<base::Vector<3,double>::Type>& nodes,
                          std::map<Edge,unsigned>& uniqueNodes,
                          std::vector<Simplex<2,unsigned>::Type>& surface,
                          std::vector<Simplex<3,unsigned>::Type>& volumeIn, 
                          std::vector<Simplex<3,unsigned>::Type>& volumeOut );
        //@}

        //----------------------------------------------------------------------
        namespace detail_{
            
            //------------------------------------------------------------------
            // Apply the right marching scheme in function of the dimension
            template<unsigned DIM> struct ApplyMarchingSimplex;

            //------------------------------------------------------------------
            // DIM = 1
            template<>
            struct ApplyMarchingSimplex<1>
            {
                static void apply( const Simplex<1,unsigned>::Type& indexSimplex,
                                   const Simplex<1,double>::Type&   distances,
                                   std::vector<base::Vector<1,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<Simplex<0,unsigned>::Type>& surface,
                                   std::vector<Simplex<1,unsigned>::Type>& volumeIn, 
                                   std::vector<Simplex<1,unsigned>::Type>& volumeOut )
                {
                    bisect( indexSimplex, distances, nodes, uniqueNodes,
                            surface, volumeIn, volumeOut );
                }
            };

            //------------------------------------------------------------------
            // DIM = 2
            template<>
            struct ApplyMarchingSimplex<2>
            {
                static void apply( const Simplex<2,unsigned>::Type& indexSimplex,
                                   const Simplex<2,double>::Type&   distances,
                                   std::vector<base::Vector<2,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<Simplex<1,unsigned>::Type>& surface,
                                   std::vector<Simplex<2,unsigned>::Type>& volumeIn, 
                                   std::vector<Simplex<2,unsigned>::Type>& volumeOut )
                {
                    marchingTri( indexSimplex, distances, nodes, uniqueNodes,
                                 surface, volumeIn, volumeOut );
                }
            };

            //------------------------------------------------------------------
            // DIM = 3
            template<>
            struct ApplyMarchingSimplex<3>
            {
                static void apply( const Simplex<3,unsigned>::Type& indexSimplex,
                                   const Simplex<3,double>::Type&   distances,
                                   std::vector<base::Vector<3,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<Simplex<2,unsigned>::Type>& surface,
                                   std::vector<Simplex<3,unsigned>::Type>& volumeIn, 
                                   std::vector<Simplex<3,unsigned>::Type>& volumeOut )
                {
                    marchingTet( indexSimplex, distances, nodes, uniqueNodes,
                                 surface, volumeIn, volumeOut );
                }
            };
        
        } // namespace detail_
    }
}

//------------------------------------------------------------------------------
/** Surface and volume reconstruction from signed distances at vertices.
 *  Main interface function to the dimension-dependent marching algorithmes.
 *  Input are the vertex values of signed distances of a unit hyper-cube, output
 *  are the nodes of the cube plus the cut points and the connectivities of the
 *  surface and volume simplices.
 *  This function decomposes the given hypercube into simplex elements and calls
 *  the dimension-dependent marching algorithm for every simplex. Using a map
 *  of lines between vertices, the duplication of nodes is avoided.
 *  For the 3D case, see the famous Marching Tetrahedrons algorithm:
 *  http://en.wikipedia.org/wiki/Marching_tetrahedrons
 *  The 2D and 1D cases are trivial reductions of that concept. Note that
 *  contrary to the original algorithm, here also the volume portions of the
 *  cut hypercube are identified and equipped with DIM-simplex elements.
 *  \tparam DIM The spatial dimension of the problem
 *  \param[in]  signedDistances   Vertex values of signed distances of hypercube
 *  \param[out] nodes             Coordinates of vertices and cut points
 *  \param[out] surface           Surface connectivity
 *  \param[out] volumeIn          Connectivity of the volume inside the domain
 *  \param[out] volumeOut         Connectivity of the volume outside the domain
 */
template<unsigned DIM>
void base::cut::marching( const boost::array<double,
                                             base::cut::HyperCube<DIM>::numVertices>&
                          signedDistances,
                          std::vector<typename base::Vector<DIM,double>::Type>& nodes,
                          std::vector<typename Simplex<DIM-1,unsigned>::Type>& surface,
                          std::vector<typename Simplex<DIM  ,unsigned>::Type>& volumeIn,
                          std::vector<typename Simplex<DIM  ,unsigned>::Type>& volumeOut )
{
    typedef base::cut::HyperCube<DIM>            HC;
    typedef base::cut::DecomposeHyperCube<DIM>   DHC;
            
    // remember the edges / vertices that are already in the
    std::map<Edge,unsigned> uniqueNodes;
            
    // fill nodes with the vertices of the HyperCube
    for ( unsigned v = 0; v < HC::numVertices; v++ ) {
        nodes.push_back( HC::vertex( v ) );
        Edge edge( v );
        uniqueNodes[ edge ] = v;
    }

    // go through all sub-simplices
    for ( unsigned s = 0; s < DHC::numSimplices; s++ ) {

        // generate index set of volume simplex and the distances
        typename Simplex<DIM,unsigned>::Type indexSimplex;
        typename Simplex<DIM,double  >::Type distances; 
        for ( unsigned v = 0; v < DIM+1; v++ ) {
            indexSimplex[v] = DHC::apply( s, v );
            distances[v] = signedDistances[ indexSimplex[v] ];
        }

        // apply marching to volume simplex
        detail_::ApplyMarchingSimplex<DIM>::apply( indexSimplex, distances, nodes,
                                                   uniqueNodes,
                                                   surface,
                                                   volumeIn, volumeOut );
    }

    return;
}

//------------------------------------------------------------------------------
// Helper functions used by the algorithms
namespace base{
    namespace cut{
        namespace detail_{

            //------------------------------------------------------------------
            /** Generate a cut point between two given vertices by linear
             *  interpolation of the signed distances. For uniqueness, the
             *  edge between the vertices is used as a key value of the cut
             *  point in a map.
             *  \tparam DIM Spatial dimension
             *  \param[in]     i1, i2  Indices of the end vertices of the edge
             *  \param[in]     d1, d2  Signed distances of these vertices
             *  \param[in,out] nodes   Vertex and cut point coordinates
             *  \param[in,out] uniqueNodes Map for uniqueness of the cut points
             */
            template<unsigned DIM>
            unsigned intersect( const unsigned i1, const unsigned i2,
                                const double   d1, const double   d2,
                                std::vector<typename base::Vector<DIM,double>::Type>& nodes,
                                std::map<Edge,unsigned>& uniqueNodes )
            {
                // New node lies on this edge
                Edge edge( i1, i2 );
                // Find edge in map of unique nodes
                std::map<Edge,unsigned>::iterator iter = uniqueNodes.find( edge );
                // If already exists pass back number of existing node
                if ( iter != uniqueNodes.end() ) {
                    return iter -> second;
                }
                // else create a new node
                else {
                    // with the number
                    const unsigned numNewNode = nodes.size();
                    // Store in map (cannot exist already!)
                    uniqueNodes[ edge ] = numNewNode;
                    // line coordinate of the intersection point by linear interpolation
                    const double xi = d1 / (d1 - d2);
                    // new node
                    const typename base::Vector<DIM,double>::Type newNode =
                        (1.-xi) * nodes[ i1 ] + xi * nodes[ i2 ];
                    // Store new node
                    nodes.push_back( newNode );
                    // pass back the number of the new node to the caller
                    return numNewNode;
                }
                return base::invalidInt;
            }

            //------------------------------------------------------------------
            //! Change orientation of a simplex
            template<unsigned DIM>
            void reverseSimplex( typename Simplex<DIM,unsigned>::Type& toBeReversed )
            {
                const unsigned tmp = toBeReversed[0];
                toBeReversed[0] = toBeReversed[1];
                toBeReversed[1] = tmp;
                return;
            }

        }
    }
}

//------------------------------------------------------------------------------
// Include implementation files
#include <base/cut/bisect.ipp>
#include <base/cut/marchingTri.ipp>
#include <base/cut/marchingTet.ipp>


#endif

