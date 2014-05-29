//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   MarchingUtils.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_cut_marchingutils_hpp
#define base_cut_marchingutils_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <utility>
#include <map>
#include <bitset>
// boost includes
#include <boost/array.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

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

        //! \cond SKIPDOX
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
        //! \endcond

        //! Simplex with double values
        template<unsigned DIM>
        struct DSimplex
        {
            typedef typename Simplex<DIM,double>::Type Type;
        };

        //! Simplex with vectors
        template<unsigned SDIM, unsigned MDIM>
        struct VecSimplex
        {
            typedef typename base::Vector<SDIM,double>::Type VecDim;
            typedef typename Simplex<MDIM,VecDim>::Type      Type;
        };

        //----------------------------------------------------------------------
        namespace detail_{
            template<unsigned DIM> struct SimplexOrient;

            //! Signed length of a 1-simplex (= line)
            template<>
            struct SimplexOrient<1>
            {
                static double apply( const VecSimplex<1,1>::Type& lin )
                {
                    return lin[1][0] - lin[0][0];
                }
            };

            //! Signed area of a 2-simplex (= triangle)
            template<>
            struct SimplexOrient<2>
            {
                static double apply( const VecSimplex<2,2>::Type& tri )
                {
                    return
                        (tri[1][0] - tri[0][0]) * (tri[2][1] - tri[0][1]) -
                        (tri[2][0] - tri[0][0]) * (tri[1][1] - tri[0][1]);
                }
            };

            //! Signed volume of a 3-simplex (= tetrahedron)
            template<>
            struct SimplexOrient<3>
            {
                static double apply( const VecSimplex<3,3>::Type& tet )
                {
                    const base::Vector<3>::Type aux =
                        base::crossProduct( tet[1] - tet[0], tet[2] - tet[0] );
                    return
                        base::dotProduct( aux, tet[3] - tet[0] );
                }
            };
            
        }

        //! Return orientation of the simplex
        template<unsigned DIM>
        double simplexOrient( const typename VecSimplex<DIM,DIM>::Type& simplex )
        {
            return detail_::SimplexOrient<DIM>::apply( simplex );
        }

        //----------------------------------------------------------------------
        namespace detail_{
            template<unsigned DIM> struct SurfSimplexSize;

            //! Size of a point = 1
            template<>
            struct SurfSimplexSize<1>
            {
                static double apply( const VecSimplex<1,0>::Type& pt )
                {
                    return 1.0;
                }
            };

            //! Length of a line segment
            template<>
            struct SurfSimplexSize<2>
            {
                static double apply( const VecSimplex<2,1>::Type& lin )
                {
                    return base::norm( lin[1] - lin[0] );
                }
            };

            //! Area of embedded triangle
            template<>
            struct SurfSimplexSize<3>
            {
                static double apply( const VecSimplex<3,2>::Type& tri )
                {
                    const base::Vector<3>::Type aux =
                        base::crossProduct( tri[1] - tri[0], tri[2] - tri[0] );
                    return base::norm( aux );
                }
            };

        }

        //! Return size of the embedded surface simplex
        template<unsigned DIM>
        double surfSimplexSize( const typename VecSimplex<DIM,DIM-1>::Type& surfSimplex )
        {
            return detail_::SurfSimplexSize<DIM>::apply( surfSimplex );
        }
        
        //----------------------------------------------------------------------
        /** Edge with comparison operator, allows unique identification of
         *  vertices and cut points.
         */
        class Edge
        {
        public:
            //! Construct from the line between two vertices
            Edge( const unsigned v1, const unsigned v2 )
                : vertices_( std::make_pair( (v1 < v2? v1 : v2),
                                             (v1 < v2? v2 : v1) ) ) { }

            //! Use build-in comparison for unique sorting
            bool operator<( const Edge& other ) const
            {
                return vertices_ < other.vertices_;
            }

            //! @name Accessor
            //{
            unsigned first() const  { return vertices_.first; }
            unsigned second() const { return vertices_.second; }
            //}

        private:
            //! (v1,v2) denotes an oriented edge
            const std::pair<unsigned,unsigned> vertices_;
        };


        //----------------------------------------------------------------------
        //! @name Dimension-dependent algorithms
        //@{
        template<unsigned SDIM>
        void bisect( const USimplex<1>::Type& indexSimplex,
                     const DSimplex<1>::Type&   distances,
                     std::vector<typename base::Vector<SDIM,double>::Type>& nodes,
                     std::map<Edge,unsigned>& uniqueNodes,
                     std::vector<USimplex<0>::Type>& surface,
                     std::vector<USimplex<1>::Type>& volumeIn, 
                     std::vector<USimplex<1>::Type>& volumeOut );

        template<unsigned SDIM>
        void marchingTri( const USimplex<2>::Type& indexSimplex,
                          const DSimplex<2>::Type&   distances,
                          std::vector<typename base::Vector<SDIM,double>::Type>& nodes,
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
            template<unsigned DIM,unsigned SHAPEDIM> struct ApplyMarchingSimplex;

            //------------------------------------------------------------------
            // DIM = 0 ---> placeholder
            template<unsigned DIM>
            struct ApplyMarchingSimplex<DIM,0>
            {
                static void apply( const USimplex<0>::Type& indexSimplex,
                                   const DSimplex<0>::Type&   distances,
                                   std::vector<typename base::Vector<DIM,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<USimplex<0>::Type>& surface,
                                   std::vector<USimplex<0>::Type>& volumeIn, 
                                   std::vector<USimplex<0>::Type>& volumeOut )
                {
                    return;
                }
            };

            //------------------------------------------------------------------
            // DIM = 1
            template<unsigned DIM>
            struct ApplyMarchingSimplex<DIM,1>
            {
                static void apply( const USimplex<1>::Type& indexSimplex,
                                   const DSimplex<1>::Type&   distances,
                                   std::vector<typename base::Vector<DIM,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<USimplex<0>::Type>& surface,
                                   std::vector<USimplex<1>::Type>& volumeIn, 
                                   std::vector<USimplex<1>::Type>& volumeOut )
                {
                    bisect<DIM>( indexSimplex, distances, nodes, uniqueNodes,
                                 surface, volumeIn, volumeOut );
                }
            };

            //------------------------------------------------------------------
            // DIM = 2
            template<unsigned DIM>
            struct ApplyMarchingSimplex<DIM,2>
            {
                static void apply( const USimplex<2>::Type& indexSimplex,
                                   const DSimplex<2>::Type&   distances,
                                   std::vector<typename base::Vector<DIM,double>::Type>& nodes,
                                   std::map<Edge,unsigned>& uniqueNodes,
                                   std::vector<USimplex<1>::Type>& surface,
                                   std::vector<USimplex<2>::Type>& volumeIn, 
                                   std::vector<USimplex<2>::Type>& volumeOut )
                {
                    marchingTri<DIM>( indexSimplex, distances, nodes, uniqueNodes,
                                      surface, volumeIn, volumeOut );
                }
            };

            //------------------------------------------------------------------
            // DIM = 3
            template<>
            struct ApplyMarchingSimplex<3,3>
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

            //------------------------------------------------------------------
            //! Interpolate the given signed distances
            template<base::Shape SHAPE, unsigned DEGREE>
            struct InterpolatedDistance
            {
                typedef typename base::LagrangeShapeFun<DEGREE,SHAPE> LagrangeFun;
                typedef typename LagrangeFun::VecDim                  VecDim;
                
                static double evaluate( const VecDim& xi,
                                        const boost::array<double,LagrangeFun::numFun>&
                                        signedDistances )
                {
                    double result = 0.;
                    
                    typename LagrangeFun::FunArray phi;
                    LagrangeFun lagFun;
                    lagFun.fun( xi, phi );
                    for ( unsigned f = 0; f < LagrangeFun::numFun; f++ )
                        result += phi[f] * signedDistances[f];

                    return result;
                }

                static VecDim gradient( const VecDim& xi,
                                        const boost::array<double,LagrangeFun::numFun>&
                                        signedDistances )
                {
                    VecDim result = base::constantVector<base::ShapeDim<SHAPE>::value>( 0. );
                    
                    typename LagrangeFun::GradArray dPhi;
                    LagrangeFun lagFun;
                    lagFun.gradient( xi, dPhi );
                    for ( unsigned f = 0; f < LagrangeFun::numFun; f++ )
                        result += signedDistances[f] * dPhi[f];
                    
                    return result;
                }

            };

        } // namespace detail_

    }
}

//------------------------------------------------------------------------------
// Include implementation files
#include <base/cut/bisect.ipp>
#include <base/cut/marchingTri.ipp>
#include <base/cut/marchingTet.ipp>



#endif
