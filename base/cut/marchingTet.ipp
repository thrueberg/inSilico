//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   marchingTet.ipp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        namespace detail_{

            //------------------------------------------------------------------
            /** Given an index i from (0,1,2,3), generate indices j, k, and el
             *  such that i-j-k-el forms a proper tetrahedron.
             *  \param[in]  i      Given first index of the vertex permutation
             *  \param[out] other  Array of the remaining three indices
             */
            void permutateTet( const unsigned i,
                               boost::array<unsigned,3>& other )
            {
                const unsigned inc = (2*i + 1)%4;
                other[0] = (i       +inc)%4;
                other[1] = (other[0]+inc)%4;
                other[2] = (other[1]+inc)%4;
                return;
            }

            //------------------------------------------------------------------
            /** Case of one node inside of the domain and 3 nodes outside.
             *  Note that this function is also used for the reverse case of
             *  precisely one node outside of the domain. In that case, the
             *  meanings of 'out' and 'in' have to be reversed in the following
             *  and the produced surface triangle needs to be re-oriented.
             *  \param[in]  I            Vertex ID of the vertex inside 
             *  \param[in]  iS           Vertex numbers of the tet to cut
             *  \param[in]  distances    Vertex distance values of that tet
             *  \param      nodes        Coordinates of vertices and cut points
             *  \param      uniqueNodes  Node uniqueness map
             *  \param[out] inSimplex    Connectivity of the inside tetrahedron
             *  \param[out] outSimplices Array of tetrahedra of the outer volume
             *  \param[out] surfSimplex  Connectivity of the surface triangle
             */
            void cutTet1( const unsigned I,
                          const USimplex<3>::Type& iS,
                          const DSimplex<3>::Type& distances,
                          std::vector<base::Vector<3,double>::Type>& nodes,
                          std::map<base::cut::Edge,unsigned>&        uniqueNodes,
                          USimplex<3>::Type&                 inSimplex,
                          boost::array<USimplex<3>::Type,3>& outSimplices,
                          USimplex<2>::Type&                 surfSimplex )
            {

                // indices of the three nodes outside
                boost::array<unsigned,3> out;
                detail_::permutateTet( I, out );
                
                // generate the intersection nodes
                boost::array<unsigned,3> cut;
                for ( unsigned c = 0; c < 3; c++ )
                    cut[c] = 
                        detail_::intersect<3>( iS[I], iS[ out[c] ],
                                               distances[I],    distances[    out[c] ],
                                               nodes, uniqueNodes );
                
                // inside is only one triangle
                inSimplex = USimplex<3>::create( iS[I], cut[0], cut[1], cut[2] );
                //        = {{ iS[I], cut[0], cut[1], cut[2] }};
                
                // outside are two triangles (choice seems arbitrary)
                outSimplices[0] =
                    USimplex<3>::create( iS[ out[0] ], iS[ out[2] ], iS[ out[1] ], cut[0] );
                outSimplices[1] =
                    USimplex<3>::create( iS[ out[1] ], iS[ out[2] ], cut[2],       cut[0] );
                outSimplices[2] =
                    USimplex<3>::create( cut[1],       iS[ out[1] ], cut[2],       cut[0] );
                //outSimplices[0] = {{ iS[ out[0] ], iS[ out[2] ], iS[ out[1] ], cut[0] }};
                //outSimplices[1] = {{ iS[ out[1] ], iS[ out[2] ], cut[2],       cut[0] }};
                //outSimplices[2] = {{ cut[1],       iS[ out[1] ], cut[2],       cut[0] }};
                
                // surface simplex
                surfSimplex = USimplex<2>::create( cut[0], cut[1], cut[2] );
                //surfSimplex = {{ cut[0], cut[1], cut[2] }};

                return;

            }

            
            //------------------------------------------------------------------
            /** Case of two nodes inside of the domain and two nodes outside.
             *  \param[in]  I1, I2       Vertex IDs of the vertices inside 
             *  \param[in]  iS           Vertex numbers of the tet to cut
             *  \param[in]  distances    Vertex distance values of that tet
             *  \param      nodes        Coordinates of vertices and cut points
             *  \param      uniqueNodes  Node uniqueness map
             *  \param[out] inSimplices  Array of tetrahedra of the inner volume  
             *  \param[out] outSimplices Array of tetrahedra of the outer volume
             *  \param[out] surfSimplex1, surfSimplex2
             *                           Connectivity of the surface triangles
             */
            void cutTet2( const unsigned I1, const unsigned I2,
                          const USimplex<3>::Type&          iS,
                          const DSimplex<3>::Type&          distances,
                          std::vector<base::Vector<3,double>::Type>& nodes,
                          std::map<base::cut::Edge,unsigned>&        uniqueNodes,
                          boost::array<USimplex<3>::Type,3>& inSimplices,
                          boost::array<USimplex<3>::Type,3>& outSimplices,
                          USimplex<2>::Type&                 surfSimplex1,
                          USimplex<2>::Type&                 surfSimplex2 )
            {
                // indices of the three nodes outside
                unsigned O1, O2;
                boost::array<unsigned,3> aux;

                // get other indices for a proper triangle
                detail_::permutateTet( I1, aux );
                // rotate these indices until the first equals the in2-index

                while ( aux[0] != I2 ) {
                    const unsigned tmp = aux[0];
                    aux[0] = aux[1];
                    aux[1] = aux[2];
                    aux[2] = tmp;
                }
                // assign remaining to the out indices
                O1 = aux[1];
                O2 = aux[2];
                // I1-I2-O1-O2 is now a proper TET

                // generate the intersection nodes
                const unsigned C1 =
                    detail_::intersect<3>( iS[I1], iS[O1],
                                           distances[I1],    distances[O1],
                                           nodes, uniqueNodes );
                const unsigned C2 =
                    detail_::intersect<3>( iS[I1], iS[O2],
                                           distances[I1],    distances[O2],
                                           nodes, uniqueNodes );
                const unsigned C3 =
                    detail_::intersect<3>( iS[I2], iS[O1],
                                           distances[I2],    distances[O1],
                                           nodes, uniqueNodes );
                const unsigned C4 =
                    detail_::intersect<3>( iS[I2], iS[O2],
                                           distances[I2],    distances[O2],
                                           nodes, uniqueNodes );

                // three inside tets
                inSimplices[0] = USimplex<3>::create( iS[I2], C4, C3, iS[I1] );
                inSimplices[1] = USimplex<3>::create( C3,     C4, C2, iS[I1] );
                inSimplices[2] = USimplex<3>::create( C1,     C3, C2, iS[I1] );
                
                // three outside tets
                outSimplices[0] = USimplex<3>::create( iS[O2], C4, C2, iS[O1] );
                outSimplices[1] = USimplex<3>::create( C2,     C4, C3, iS[O1] );
                outSimplices[2] = USimplex<3>::create( C1,     C2, C3, iS[O1] );
                
                // surface simplex
                surfSimplex1 = USimplex<2>::create( C1, C2, C3 );
                surfSimplex2 = USimplex<2>::create( C2, C4, C3 );

                // // three inside tets
                // inSimplices[0] = {{ iS[I2], C4, C3, iS[I1] }};
                // inSimplices[1] = {{ C3,     C4, C2, iS[I1] }};
                // inSimplices[2] = {{ C1,     C3, C2, iS[I1] }};
                // 
                // // three outside tets
                // outSimplices[0] = {{ iS[O2], C4, C2, iS[O1] }};
                // outSimplices[1] = {{ C2,     C4, C3, iS[O1] }};
                // outSimplices[2] = {{ C1,     C2, C3, iS[O1] }};
                // 
                // // surface simplex
                // surfSimplex1 = {{ C1, C2, C3 }};
                // surfSimplex2 = {{ C2, C4, C3 }};

                return;
            }

        } // namespace detail_
    }
}

//------------------------------------------------------------------------------
/** Given the four values of a signed distance function of a tetrahedron,
 *  recover a surface triangulation of the implicit surface and create
 *  tetrahedra of the volumes in- and outside of that surface.
 *  A first simplification is taken by assuming that the signed distances are
 *  never exactly zero. This can easily achieved by shifting them away from
 *  zero by a numerical epsilon. Then, by symmetry only five cases are possible:
 *
 *    1.  All nodes are inside
 *    2.  All nodes are outside
 *    3.  One node is inside
 *    4.  Two nodes are inside
 *    5.  Three nodes are inside
 *
 *  The cases 1 and 2 are trivial and the cases 3 and 5 are complementary to
 *  each other with the meaning of in- and outside reversed and the surface
 *  triangles re-oriented. Therefore, two main cases are identified:
 *
 *    1.  The implicit surface separates one node from the three others
 *    2.  The implicit surface separates two nodes from the other two
 *
 *  In the first case, the tetrahedron is split into one tetrahedron on one
 *  side of surface, a surface triangle and a wedge on the other side which
 *  itself is decomposed into three tetrahedra.
 *  \image html tet1.png
 *  In the second case, the surface is a quadrilateral divided into two
 *  triangles and separates two wedges which lie on either side of it. These
 *  wedges are themselves divided into three tetrahedra.
 *  \image html tet2.png
 * 
 *  \param[in]     indexSimplex   Vertex IDs of the tetrahedron
 *  \param[in]     distances      Values of the signed distance function
 *  \param[in,out] nodes          Coordinates of the vertices and cut points
 *  \param[in,out] uniqueNodes    Map for avoiding node duplications
 *  \param[out]    surface        Connectivity of the surface elements
 *  \param[out]    volumeIn       Connectivity of the volume inside the domain
 *  \param[out]    volumeOut      Connectivity of the volume outside the domain
 */
void base::cut::marchingTet( const base::cut::USimplex<3>::Type& indexSimplex,
                             const base::cut::DSimplex<3>::Type& distances,
                             std::vector<base::Vector<3,double>::Type>&  nodes,
                             std::map<base::cut::Edge,unsigned>&         uniqueNodes,
                             std::vector<base::cut::USimplex<2>::Type>& surface,
                             std::vector<base::cut::USimplex<3>::Type>& volumeIn, 
                             std::vector<base::cut::USimplex<3>::Type>& volumeOut )
{
    // flags checking the values of the distances
    std::bitset<4> flags;
    for ( unsigned d = 0; d < 4; d++ )
        flags[d] = (distances[d] >= 0.);

    //--------------------------------------------------------------------------
    // case I: all nodes are inside
    if      ( flags.count() == 4 ) {
        volumeIn.push_back( indexSimplex );
    }
    //--------------------------------------------------------------------------
    // case II: all nodes are outside
    else if ( flags.count() == 0 ) {
        volumeOut.push_back( indexSimplex );
    }
    //--------------------------------------------------------------------------
    // case III: one node is inside
    else if ( flags.count() == 1 ) {
        // index of the one vertex which is inside
        const unsigned in =
            (flags.test(0) ? 0 : (flags.test(1) ? 1 : (flags.test(2) ? 2 : 3)));
        // indices of the outside vertices, such that in-out1-out2-out3 is proper
        // 
        base::cut::USimplex<3>::Type                   inSimplex;
        boost::array<base::cut::USimplex<3>::Type,3>  outSimplices;
        base::cut::USimplex<2>::Type                 surfSimplex;
        detail_::cutTet1( in, indexSimplex, distances, nodes, uniqueNodes,
                          inSimplex, outSimplices, surfSimplex );
        
        volumeIn.push_back(   inSimplex );

        for ( unsigned os = 0; os < 3; os++ )
            volumeOut.push_back( outSimplices[os] );
        
        surface.push_back(  surfSimplex );

    }
    //--------------------------------------------------------------------------
    // case IV:  two nodes are inside
    else if ( flags.count() == 2 ) {
        // get indices of the first and second vertex inside the domain
        unsigned in1,  in2;
        bool foundFirst = false;
        for ( unsigned d = 0; d < 4; d ++ ) {
            if ( flags.test(d) and not foundFirst ) { in1 = d; foundFirst = true; }
            if ( flags.test(d) and     foundFirst )   in2 = d;
        }

        boost::array<base::cut::USimplex<3>::Type,3> inSimplices, outSimplices;
        base::cut::USimplex<2>::Type                surfSimplex1, surfSimplex2;
        detail_::cutTet2( in1, in2, indexSimplex, distances, nodes, uniqueNodes,
                          inSimplices, outSimplices, surfSimplex1, surfSimplex2 );

        for ( unsigned is = 0; is < 3; is++ )  volumeIn.push_back(  inSimplices[is] );
        for ( unsigned os = 0; os < 3; os++ ) volumeOut.push_back( outSimplices[os] );
        
        surface.push_back(  surfSimplex1 );
        surface.push_back(  surfSimplex2 );

    }
    //--------------------------------------------------------------------------
    // case V:  three nodes are inside
    else {
        
        // index of the one vertex which is inside
        const unsigned out =
            (flags.test(0) == false ? 0 : (flags.test(1) == false ? 1 :
                                          (flags.test(2) == false  ? 2 : 3)));

        base::cut::USimplex<3>::Type                  outSimplex;
        boost::array<base::cut::USimplex<3>::Type,3>   inSimplices;
        base::cut::USimplex<2>::Type                 surfSimplex;
        detail_::cutTet1( out, indexSimplex, distances, nodes, uniqueNodes,
                          outSimplex, inSimplices, surfSimplex );
        
        volumeOut.push_back( outSimplex );

        for ( unsigned is = 0; is < 3; is ++ ) 
            volumeIn.push_back(   inSimplices[is] );

        detail_::reverseSimplex<2>( surfSimplex );        
        surface.push_back(  surfSimplex );

    }
    return;
}

