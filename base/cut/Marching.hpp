//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Marching.hpp
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

        template<unsigned DIM> struct MarchingCube;
        template<unsigned DIM> struct MarchingSimplex;

        //------------------------------------------------------------------------
        //! Marching scheme in function of the shape
        template<base::Shape SHAPE>
        struct Marching
            : base::IfElse<SHAPE==base::SimplexShape<ShapeDim<SHAPE>::value>::value,
                           MarchingSimplex<ShapeDim<SHAPE>::value>,
                           MarchingCube<   ShapeDim<SHAPE>::value> >::Type
        { };

        
        //----------------------------------------------------------------------
        //! Define a simplex to be a fixed-size array
        template<unsigned DIM, typename VAL=unsigned>
        struct Simplex
        {
            typedef boost::array<VAL,DIM+1> Type;
        };

        template<unsigned DIM>
        struct USimplex
        {
            typedef typename Simplex<DIM,unsigned>::Type Type;

            static Type create( const unsigned i1,
                                const unsigned i2 = base::invalidInt,
                                const unsigned i3 = base::invalidInt,
                                const unsigned i4 = base::invalidInt );
        };

        template<>
        USimplex<0>::Type USimplex<0>::create( const unsigned i1, const unsigned i2,
                                               const unsigned i3, const unsigned i4 )
        {
            Type result;
            result[0] = i1;
            return result;
        }

        template<>
        USimplex<1>::Type USimplex<1>::create( const unsigned i1, const unsigned i2,
                                               const unsigned i3, const unsigned i4 )
        {
            Type result;
            result[0] = i1; result[1] = i2;
            return result;
        }

        template<>
        USimplex<2>::Type USimplex<2>::create( const unsigned i1, const unsigned i2,
                                               const unsigned i3, const unsigned i4 )
        {
            Type result;
            result[0] = i1; result[1] = i2; result[2] = i3;
            return result;
        }

        template<>
        USimplex<3>::Type USimplex<3>::create( const unsigned i1, const unsigned i2,
                                               const unsigned i3, const unsigned i4 )
        {
            Type result;
            result[0] = i1; result[1] = i2; result[2] = i3; result[3] = i4;
            return result;
        }


        template<unsigned DIM>
        struct DSimplex
        {
            typedef typename Simplex<DIM,double>::Type Type;
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
            //! (v1,v2) denotes an oriented edge
            const std::pair<unsigned,unsigned> vertices_;
        };

        //----------------------------------------------------------------------
        //! @name Dimension-dependent algorithms
        //@{
        void bisect( const USimplex<1>::Type& indexSimplex,
                     const DSimplex<1>::Type&   distances,
                     std::vector<base::Vector<1,double>::Type>& nodes,
                     std::map<Edge,unsigned>& uniqueNodes,
                     std::vector<USimplex<0>::Type>& surface,
                     std::vector<USimplex<1>::Type>& volumeIn, 
                     std::vector<USimplex<1>::Type>& volumeOut );

        void marchingTri( const USimplex<2>::Type& indexSimplex,
                          const DSimplex<2>::Type&   distances,
                          std::vector<base::Vector<2,double>::Type>& nodes,
                          std::map<Edge,unsigned>& uniqueNodes,
                          std::vector<USimplex<1>::Type>& surface,
                          std::vector<USimplex<2>::Type>& volumeIn, 
                          std::vector<USimplex<2>::Type>& volumeOut );

        void marchingTet( const USimplex<3>::Type& indexSimplex,
                          const DSimplex<3>::Type&   distances,
                          std::vector<base::Vector<3,double>::Type>& nodes,
                          std::map<Edge,unsigned>& uniqueNodes,
                          std::vector<USimplex<2>::Type>& surface,
                          std::vector<USimplex<3>::Type>& volumeIn, 
                          std::vector<USimplex<3>::Type>& volumeOut );
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
                static void apply( const USimplex<1>::Type& indexSimplex,
                                   const DSimplex<1>::Type&   distances,
                                   std::vector<base::Vector<1,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<USimplex<0>::Type>& surface,
                                   std::vector<USimplex<1>::Type>& volumeIn, 
                                   std::vector<USimplex<1>::Type>& volumeOut )
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
                static void apply( const USimplex<2>::Type& indexSimplex,
                                   const DSimplex<2>::Type&   distances,
                                   std::vector<base::Vector<2,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<USimplex<1>::Type>& surface,
                                   std::vector<USimplex<2>::Type>& volumeIn, 
                                   std::vector<USimplex<2>::Type>& volumeOut )
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
                static void apply( const USimplex<3>::Type& indexSimplex,
                                   const DSimplex<3>::Type&   distances,
                                   std::vector<base::Vector<3,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<USimplex<2>::Type>& surface,
                                   std::vector<USimplex<3>::Type>& volumeIn, 
                                   std::vector<USimplex<3>::Type>& volumeOut )
                {
                    marchingTet( indexSimplex, distances, nodes, uniqueNodes,
                                 surface, volumeIn, volumeOut );
                }
            };

            //------------------------------------------------------------------
            
        
        } // namespace detail_

    }
}

//------------------------------------------------------------------------------
/** Reconstruct surface and partial volumina for a hyper-cube.
 *  Note that the input and output arrays are based on a lexicographic ordering.
 *  \tparam DIM Dimension of the shape manifold
 */
template<unsigned DIM>
struct base::cut::MarchingCube
{
    static const unsigned dim         = DIM;
    static const unsigned numVertices = base::cut::HyperCube<dim>::numVertices;

    typedef typename base::Vector<dim,double>::Type     VecDim;
    typedef typename base::cut::USimplex<dim  >::Type   VolSimplex;
    typedef typename base::cut::USimplex<dim-1>::Type   SurfSimplex;

    static void apply( const boost::array<double,numVertices>& signedDistances,
                       std::vector<VecDim>& nodes,
                       std::vector<SurfSimplex>& surface,
                       std::vector<VolSimplex>&  volumeIn,
                       std::vector<VolSimplex>&  volumeOut )
    {
        typedef base::cut::HyperCube<dim>            HC;
        typedef base::cut::DecomposeHyperCube<dim>   DHC;
            
        // remember the edges / vertices that are already in the
        std::map<base::cut::Edge,unsigned> uniqueNodes;
            
        // fill nodes with the vertices of the HyperCube
        for ( unsigned v = 0; v < numVertices; v++ ) {
            // store 
            nodes.push_back( HC::vertex(v) );
            base::cut::Edge edge( v );
            uniqueNodes[ edge ] = v;
        }

        // go through all sub-simplices
        for ( unsigned s = 0; s < DHC::numSimplices; s++ ) {

            // generate index set of volume simplex and the distances
            VolSimplex                   indexSimplex;
            typename DSimplex<dim>::Type distances;
            
            for ( unsigned v = 0; v < dim+1; v++ ) {
                indexSimplex[v] = DHC::apply( s, v );
                distances[v] = signedDistances[ indexSimplex[v] ];
            }

            // apply marching to volume simplex
            detail_::ApplyMarchingSimplex<dim>::apply( indexSimplex, distances, nodes,
                                                       uniqueNodes,
                                                       surface,
                                                       volumeIn, volumeOut );
        }

        return;
    }

};

//------------------------------------------------------------------------------
template<unsigned DIM>
struct base::cut::MarchingSimplex
{
    static const unsigned dim         = DIM;
    static const unsigned numVertices = dim+1;

    typedef typename base::Vector<dim,double>::Type     VecDim;
    typedef typename base::cut::USimplex<dim  >::Type   VolSimplex;
    typedef typename base::cut::USimplex<dim-1>::Type   SurfSimplex;

    typedef base::LagrangeShapeFun<1,base::SimplexShape<dim>::value>
    LinearLagrange;

    static void apply( const boost::array<double,numVertices>& signedDistances,
                       std::vector<VecDim>& nodes,
                       std::vector<SurfSimplex>& surface,
                       std::vector<VolSimplex>&  volumeIn,
                       std::vector<VolSimplex>&  volumeOut )
    {
        // remember the edges / vertices that are already in the
        std::map<base::cut::Edge,unsigned> uniqueNodes;

        // get the vertices of the simplex
        typename base::cut::Simplex<dim,VecDim>::Type vertices;
        LinearLagrange::supportPoints( vertices );
        
        // fill nodes with the vertices of the HyperCube
        for ( unsigned v = 0; v < numVertices; v++ ) {
            // store
            nodes.push_back( vertices[v] );
            base::cut::Edge edge( v );
            uniqueNodes[ edge ] = v;
        }

        // generate index set of volume simplex and the distances
        VolSimplex                   indexSimplex;
        typename DSimplex<dim>::Type distances;
            
        for ( unsigned v = 0; v < numVertices; v++ ) {
            indexSimplex[v] = v; 
            distances[v] = signedDistances[ v ];
        }

        // apply marching to volume simplex
        detail_::ApplyMarchingSimplex<dim>::apply( indexSimplex, distances, nodes,
                                                   uniqueNodes,
                                                   surface,
                                                   volumeIn, volumeOut );
        return;
    }

};


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
                                std::map<base::cut::Edge,unsigned>& uniqueNodes )
            {
                // New node lies on this edge
                base::cut::Edge edge( i1, i2 );
                // Find edge in map of unique nodes
                std::map<base::cut::Edge,unsigned>::iterator iter = uniqueNodes.find( edge );
                // If already exists pass back number of existing node
                if ( iter != uniqueNodes.end() ) {
                    return iter -> second;
                }
                // else create a new node
                else {
                    // with the number
                    const unsigned numNewNode = static_cast<unsigned>( nodes.size() );
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
            void reverseSimplex( typename USimplex<DIM>::Type& toBeReversed )
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

