//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   HierarchicOrder.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_hierarchicorder_hpp
#define base_mesh_hierarchicorder_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
// base includes
#include <base/numbers.hpp>
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace mesh{

        //----------------------------------------------------------------------
        namespace detail_{

            //! Specialised construction methods
            template<unsigned DEGREE, unsigned DIM, bool ISHYPERCUBE>
            struct BuildHierarchicLookUpTable;

            template<unsigned DEGREE>
            struct BuildHierarchicLookUpTable<DEGREE,1,true>;

            template<unsigned DEGREE>
            struct BuildHierarchicLookUpTable<DEGREE,2,true>;

            template<unsigned DEGREE>
            struct BuildHierarchicLookUpTable<DEGREE,3,true>;

            template<unsigned DEGREE, unsigned DIM>
            struct BuildHierarchicLookUpTable<DEGREE,DIM,false>
            {
                static const unsigned numNodes = base::Binomial<DEGREE+DIM,
                                                                DEGREE>::value;
                
                //--------------------------------------------------------------
                //static void apply( boost::array<unsigned,numTotal> & lookUpTable )
                static boost::array<unsigned,numNodes> apply(  )
                {
                    boost::array<unsigned,numNodes> lookUpTable;
                    for ( unsigned i = 0; i < numNodes; i++ )
                        lookUpTable[i] = i;
                    return lookUpTable;
                }

            };

        }

        //----------------------------------------------------------------------
        /** Look-up table for conversion from lexicographic to hierarchic order.
         *  The Lagrange elements with the geometry shape LINE, QUAD, and HEX
         *  can be equipped with shape functions in lexicographic order. This is
         *  done in order to make use of the tensor-product structure. But
         *  when generating and handling degrees of freedom, it has to be
         *  differentiated between shape functions associated with a VERTEX,
         *  an EDGE, a FACE or the CELL interior. Therefore, this object
         *  provides a simple array \f$ H \f$ such that
         *  \f[
         *        h(i) = H[ i ]
         *  \f]
         *  where \f$ i \f$ is the lexicographic (linear) index of a shape
         *  function and \f$ h(i) \f$ the corresponding index in a hierarchic
         *  ordering.
         *  \tparam SHAPE    Type of geometry shape
         *  \tparam DEGREE   Polynomial degree of the shape functions
         */
        template<base::Shape SHAPE, unsigned DEGREE>
        class HierarchicOrder
        {
        public:
            static const unsigned shapeDim = base::ShapeDim<SHAPE>::value;
            
            //! Lookup Table for external access            
            static const unsigned numNodes =
                ( SHAPE == base::HyperCubeShape<shapeDim>::value ?
                  base::MToTheN<DEGREE+1,shapeDim>::value :
                  base::Binomial<DEGREE+shapeDim,DEGREE>::value );

            //! Access operation
            static unsigned apply( const unsigned i )
            {
                return lookUpTable_[i];
            }

        private:
            typedef boost::array<unsigned,numNodes> LookUpTable;
            
            static const boost::array<unsigned,numNodes> lookUpTable_;
        };

        //! Initialise lookUp table here, as it is a static member variable
        template<base::Shape SHAPE, unsigned DEGREE>
        const typename HierarchicOrder<SHAPE,DEGREE>::LookUpTable
        HierarchicOrder<SHAPE,DEGREE>::lookUpTable_
        = detail_::BuildHierarchicLookUpTable<
            DEGREE, base::ShapeDim<SHAPE>::value,
            (SHAPE ==
             base::HyperCubeShape<base::ShapeDim<SHAPE>::value>::value ) >::apply();
    
    }
}

//------------------------------------------------------------------------------
/** Construct hierarchic ordering for a line element. 
 *  Here, the ordering goes from (K=degree) the lexicographic
 *  version
 *  <pre>
 *        *-----o-----o-----o-----o ... o-----*
 *        0     1     2     3     4    K-1    K
 *  </pre>
 *  to the hierarchic version
 *  <pre>
 *        *-----o-----o-----o-----o ... o-----*
 *        0     2     3     4     5     K     1
 *  </pre>
 */
template<unsigned DEGREE>
struct base::mesh::detail_::BuildHierarchicLookUpTable<DEGREE,1,true>
{
    //! Number of total nodes along the LINE
    static const int numTotal = DEGREE + 1;

    //--------------------------------------------------------------------------
    static boost::array<unsigned,numTotal> apply()
    {
        boost::array<unsigned,numTotal> lookUpTable;
            
        // First vertex
        lookUpTable[0] = 0;

        // In-between nodes along the only edge
        for ( int i = 1; i < numTotal-1; i++ ) lookUpTable[i] = i+1;

        // Second vertex
        if ( DEGREE > 0 ) lookUpTable[ numTotal-1 ] = 1;

        return lookUpTable;
    }
};

//------------------------------------------------------------------------------
/** Construct hierarchic ordering for a quadrilateral element
 *  Here, the ordering goes from (K=degree) the lexicographic
 *  version
 *  <pre>
 *         
 *          K(K+1)                      K(K+2)
 *              *-----o-----o ... o-----* 
 *              |                       |
 *              |                       |
 *   (K-1)(K+1) 0     x     x ... x     0 K^2+K-1
 *              �                       �
 *              
 *              �                       �
 *         2K+2 0     x     x ... x     0 3K+2
 *              |                       |
 *              |    K+2   K+3          |
 *          K+1 0     x     x ... x     0 2K+1
 *              |                       |
 *              |                       |
 *              *-----o-----o ... o-----*
 *             0      1     2    K-1     K
 *
 *  </pre>
 *  to the hierarchic version
 *  <pre>
 *                           <--
 *             3     3K   3K-1   2K+2    2
 *              *-----o-----o ... o-----* 
 *              |         [2]=(2,3)     |
 *              |  K^2+K+2      K(K+2)  |
 *         3K+1 0     x     x ... x     0 2K+1   
 *              �                       �        
 *   [3]=(3,0)                                  [1]=(1,2)    
 *              �   5K-1   5K           �        
 *         4K-2 0     x     x ... x     0 K+4
 *              |                       |        
 *              |    4K   4K+1  5K-2    |        
 *         4K-1 0     x     x ... x     0 K+3    
 *              |                       |        
 *              |       [0]=(0,1)       |
 *              *-----o-----o ... o-----*
 *             0      4     5    K+2     1
 *                       -->
 *  </pre>
 */
template<unsigned DEGREE>
struct base::mesh::detail_::BuildHierarchicLookUpTable<DEGREE,2,true>
{
    //! Number of total nodes in a QUAD
    static const unsigned numTotal = (DEGREE+1) * (DEGREE+1);

    //! For legibility
    static const int K = DEGREE;
    
    //--------------------------------------------------------------------------
    static boost::array<unsigned,numTotal> apply()
    {
        boost::array<unsigned,numTotal> lookUpTable;
        
        lookUpTable.assign( base::invalidInt );

        // Array of vertex indices
        const boost::array<int,4> vertices = {{ 0, K, K*(K+2), K*(K+1) }};

        // Edges which form a ring
        const boost::array<std::pair<unsigned,unsigned>,4> edges =
            {{ std::make_pair( 0, 1 ), std::make_pair( 1, 2 ),
               std::make_pair( 2, 3 ), std::make_pair( 3, 0 ) }};

        // Global counter
        unsigned ctr = 0;

        // Insert vertices first
        for ( unsigned i = 0; i < 4; i ++ ) {
            lookUpTable[ static_cast<unsigned>(vertices[i]) ] = ctr++;
        }

        //----------------------------------------------------------------------
        // Go through edges
        for ( unsigned i = 0; i < 4; i ++ ) {

            // Begin and end vertices of edge
            const int ev1 = vertices[ edges[i].first  ];
            const int ev2 = vertices[ edges[i].second ];
            // Increment along edge
            const int delta = (ev2 - ev1) / K;

            for ( int n = 0; n < K-1; n ++ ) {
                
                const int j = (ev1+delta) + n * delta;
                lookUpTable[ static_cast<unsigned>(j) ] = ctr++;
            }

        }

        //----------------------------------------------------------------------
        // Go through face

        // Vertices spanning the face
        const int fv1 = vertices[0];
        const int fv2 = vertices[1];
        const int fv3 = vertices[3];

        // Increments per direction of the face
        const int delta1 = (fv2 - fv1) / K;
        const int delta2 = (fv3 - fv1) / K;

        for ( int n2 = 0; n2 < K-1; n2++ ) {
            for ( int n1 = 0; n1 < K-1; n1++ ) {
                    
                const int i = (fv1+delta1+delta2) + n1 * delta1 + n2 * delta2;
                lookUpTable[ static_cast<unsigned>(i) ] = ctr++;
            }
        }
        //----------------------------------------------------------------------
        return lookUpTable;
    }
};

//------------------------------------------------------------------------------
/** Construct hierarchic ordering for a hexahedral element.
 *  Here, the ordering goes from (K=DEGREE) the lexicographic
 *  version
 *  <pre>
 *
 *      
 *  K(K+1)^2 o----o----o----o......* (K+1)^3 - (K+1)
 *           �\                     \   
 *           � 0    x    x    x      0  
 *           �  \                     \ 
 *  2(K+1)^2 0   0     x    x    x     0
 *           |    \                     \
 *           | x   `                     `
 *           |      `                     `
 *   (K+1)^2 0   x   *----o----o----o......* (K+1)^3-1
 *           |       �                     �
 *           | x     �                     �
 *           |       �                     �
 *         0 *   x   0    x    x    x      0
 *            \      |                     |
 *           1 0     |                     |
 *              \    |                     |    
 *             2 0   0    x    x    x      0
 *                \  |                     |
 *                 ` |                     |
 *                . `|                     |
 *                   *----o----o----o......*
 *                  K    2K+1 3K+2          K(K+2)
 *      
 *  </pre>
 *  to the hierarchic version
 *  <pre>
 *
 *      
 *         4 o----o----o----o......* 7
 *           �\           [7]       \   
 *           � 0    x    x    x      0  
 *           �  \                     \ 
 *           0   0     x    x    x     0
 *        [8]|    \[4]     {1}       [6]\
 *           | x   `                     `
 *           |  {2} `            [5]      `
 *           0   x   *----o----o----o......* 6
 *           |       �                     �
 *           | x     �        {3}          �
 *           |       �                     �
 *         0 *   x   0    x    x    x      0
 *            \      |[9]              [10]|
 *           8 0     |                     |
 *              \    |                     |
 *       [0]   9 0   0    x    x    x      0
 *                \  |                     |
 *                 ` |                     |
 *                . `|                     |
 *                   *----o----o----o......*
 *                  1    K+7  K+8           2     
 *                               [1]
 *  </pre>
 *  For clarity, the following ordering is used
 *  <pre>
 *         V = { 0, 1, 2, 3, 4, 5, 6, 7 }
 *         E = { (0,1), (1,2), (2,3), (3,0),
 *               (4,5), (5,6), (6,7), (7,4),
 *               (0,4), (1,5), (2,6), (3,7 ) }
 *         F = { (0,1,2,3), (4,5,6,7),
 *               (0,1,5,4), (1,2,6,5),
 *               (2,3,7,6), (3,0,4,7) }
 *  </pre>
 */
template<unsigned DEGREE>
struct base::mesh::detail_::BuildHierarchicLookUpTable<DEGREE,3,true>
{
    //! Number of total nodes in a QUAD
    static const unsigned numTotal = (DEGREE+1) * (DEGREE+1) * (DEGREE+1);

    //! For legibilty
    static const int K = DEGREE;

    //--------------------------------------------------------------------------
    static boost::array<unsigned,numTotal> apply()
    {
        boost::array<unsigned,numTotal> lookUpTable;
            
        lookUpTable.assign( base::invalidInt );

        // Array of vertex indices
        const unsigned dZ = K * (K+1) * (K+1);
        const boost::array<int,8> vertices =
            {{ 0,  K,    K*(K+2),    K*(K+1),
               dZ, K+dZ, K*(K+2)+dZ, K*(K+1)+dZ }};

        // Edges which form a ring
        const boost::array<std::pair<unsigned,unsigned>,12> edges =
            {{ std::make_pair( 0, 1 ), std::make_pair( 1, 2 ),
               std::make_pair( 2, 3 ), std::make_pair( 3, 0 ),
               std::make_pair( 4, 5 ), std::make_pair( 5, 6 ),
               std::make_pair( 6, 7 ), std::make_pair( 7, 4 ),
               std::make_pair( 0, 4 ), std::make_pair( 1, 5 ),
               std::make_pair( 2, 6 ), std::make_pair( 3, 7 ) }};

        // Faces defined by 3 vertices which span them
        const boost::array< boost::tuple<unsigned,unsigned,unsigned>, 6 > faces =
            {{ boost::make_tuple( 0, 1, 3 ), boost::make_tuple( 4, 5, 7 ),
               boost::make_tuple( 0, 1, 4 ), boost::make_tuple( 1, 2, 5 ),
               boost::make_tuple( 2, 3, 6 ), boost::make_tuple( 3, 0, 7 ) }};

        // Global counter
        unsigned ctr = 0;

        // Insert vertices first
        for ( unsigned i = 0; i < 8; i ++ ) {
            lookUpTable[ vertices[i] ] = ctr++;
        }

        //----------------------------------------------------------------------
        // Go through edges
        for ( unsigned i = 0; i < 12; i ++ ) {

            // Begin and end vertices of edge
            const int ev1 = vertices[ edges[i].first  ];
            const int ev2 = vertices[ edges[i].second ];
            // Increment along edge
            const int delta = (ev2 - ev1) / K;

            for ( int n = 0; n < K-1; n ++ ) {
                
                const unsigned j = (ev1+delta) + n * delta;
                lookUpTable[ j ] = ctr++;
            }

        }

        //----------------------------------------------------------------------
        // Go through faces
        for ( unsigned i = 0; i < 6; i ++ ) {

            // Vertices spanning the face
            const int fv1 = vertices[ faces[i].get<0>() ];
            const int fv2 = vertices[ faces[i].get<1>() ];
            const int fv3 = vertices[ faces[i].get<2>() ];

            // Increments per direction of the face
            const int delta1 = (fv2 - fv1) / K;
            const int delta2 = (fv3 - fv1) / K;

            for ( int n2 = 0; n2 < K-1; n2++ ) {
                for ( int n1 = 0; n1 < K-1; n1++ ) {
                    
                    const unsigned i = (fv1+delta1+delta2) + n1 * delta1 + n2 * delta2;
                    lookUpTable[ i ] = ctr++;
                }
            }
        }

        //----------------------------------------------------------------------
        // Go through cell interior
        const int cv1 = vertices[0];
        const int cv2 = vertices[1];
        const int cv3 = vertices[3];
        const int cv4 = vertices[4];

        const int delta1 = (cv2 - cv1) / K;
        const int delta2 = (cv3 - cv1) / K;
        const int delta3 = (cv4 - cv1) / K;

        for ( int n3 = 0; n3 < K-1; n3++ ) {
            for ( int n2 = 0; n2 < K-1; n2++ ) {
                for ( int n1 = 0; n1 < K-1; n1++ ) {

                    const unsigned i = (cv1+delta1+delta2+delta3) +
                        n1 * delta1 + n2 * delta2 + n3 * delta3;

                    lookUpTable[i] = ctr++;
                }
            }
        }
        //----------------------------------------------------------------------
        return lookUpTable;
    }
};

#endif
