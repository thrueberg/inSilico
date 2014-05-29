//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   marchingTri.ipp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        namespace detail_{

            //------------------------------------------------------------------
            /** Cut-element situation for the triangle.
             *
             *                      I
             *                      x          d>0
             *                     / \
             *                  C1/   \C2
             *            -------+-----+-----  d=0
             *                  /    ´  \
             *                 /   ´     \
             *                /  ´        \ 
             *               x-------------x
             *              O1             O2  d<0
             *
             *    Inside triangle:  { I, C1, C2 }
             *   Outside triangles: { C1, O1, C2 } and { O1, O2, C2 }
             *   Surface line:      { C1, C2 }
             *
             *   Note that there is a complemenatary case in which one vertex
             *   is outside and two are inside. This function can be called in
             *   exactly the same way with the meanings of 'in' and 'out' reversed.
             *   But in that case the surface line needs to be reversed !!
             */
            template<unsigned SDIM>
            void cutTriangle( const unsigned I,
                              const USimplex<2>::Type& indexSimplex,
                              const DSimplex<2>::Type& distances,
                              std::vector<typename base::Vector<SDIM,double>::Type>& nodes,
                              std::map<base::cut::Edge,unsigned>&        uniqueNodes,
                              USimplex<2>::Type&  inSimplex,
                              USimplex<2>::Type& outSimplex1,
                              USimplex<2>::Type& outSimplex2,
                              USimplex<1>::Type& surfSimplex )
            {
                // indices of the two nodes outside
                const unsigned O1 = (I + 1) % 3;
                const unsigned O2 = (I + 2) % 3;
                
                // generate the intersection nodes
                const unsigned C1 =
                    detail_::intersect<SDIM>( indexSimplex[I], indexSimplex[O1],
                                              distances[I],    distances[O1],
                                              nodes, uniqueNodes );
                const unsigned C2 =
                    detail_::intersect<SDIM>( indexSimplex[I], indexSimplex[O2],
                                              distances[I],    distances[O2],
                                              nodes, uniqueNodes );
                
                // inside is only one triangle
                inSimplex = USimplex<2>::create( indexSimplex[I], C1, C2 );
                // outside are two triangles (choice seems arbitrary)
                outSimplex1 = USimplex<2>::create( C1, indexSimplex[O1], C2 );
                outSimplex2 = USimplex<2>::create( indexSimplex[O1], indexSimplex[O2], C2 );
                // surface simplex
                surfSimplex = USimplex<1>::create( C1, C2 ); 

                return;
            }
            
        }
    }
}

//------------------------------------------------------------------------------
/** Given the three values of a signed distance function of a triangle,
 *  this function generates the line element along the implicit surface and
 *  the triangles which form the volumes in- and out-side of the domain.
 *  Under the assumption that no signed distance value is precisely zero, the
 *  following four cases occur:
 *
 *    1. All nodes are inside
 *    2. All nodes are outside 
 *    3. One node is inside
 *    4. Two  nodes are inside
 *
 *  Simple observation shows that cases 3 and 4 are complemenatary to each other
 *  with the meanings of inside and outside reversed and the surface line
 *  re-oriented. Since the first two cases are trivial, it remains to handle
 *  the non-trivial case in which the implicit surface separates one node from
 *  the other two. This yields one triangle formed by the lonely node and the
 *  two cut points, a line between these cut points and a quadrilateral on the
 *  other side with the cut points and the two other nodes. This quadrilateral
 *  is divided into two triangles.
 *  \image html tri.png
 *
 *                             
 *  \param[in]     indexSimplex   Vertex IDs of the tetrahedron
 *  \param[in]     distances      Values of the signed distance function
 *  \param[in,out] nodes          Coordinates of the vertices and cut points
 *  \param[in,out] uniqueNodes    Map for avoiding node duplications
 *  \param[out]    surface        Connectivity of the surface elements
 *  \param[out]    volumeIn       Connectivity of the volume inside the domain
 *  \param[out]    volumeOut      Connectivity of the volume outside the domain
 */
template<unsigned SDIM>
void base::cut::marchingTri( const base::cut::USimplex<2>::Type& indexSimplex,
                             const base::cut::DSimplex<2>::Type&   distances,
                             std::vector<typename base::Vector<SDIM,double>::Type>&  nodes,
                             std::map<base::cut::Edge,unsigned>&         uniqueNodes,
                             std::vector<base::cut::USimplex<1>::Type>& surface,
                             std::vector<base::cut::USimplex<2>::Type>& volumeIn, 
                             std::vector<base::cut::USimplex<2>::Type>& volumeOut )
{
    // flags checking the values of the distances
    std::bitset<3> flags;
    for ( unsigned d = 0; d < 3; d++ )
        flags[d] = (distances[d] >= 0.);

    // case I: all nodes are inside
    if      ( flags.count() == 3 ) {
        volumeIn.push_back( indexSimplex );
    }
    // case II: all nodes are outside
    else if ( flags.count() == 0 ) {
        volumeOut.push_back( indexSimplex );
    }
    // case III: one node is inside
    else if ( flags.count() == 1 ) {
        
        // index of the node which is inside
        const unsigned in = ( flags.test(0) ? 0 : (flags.test(1) ? 1 : 2 ) );

        USimplex<2>::Type inSimplex, outSimplex1, outSimplex2;
        USimplex<1>::Type surfSimplex;
        detail_::cutTriangle<SDIM>( in, indexSimplex, distances, nodes, uniqueNodes,
                                    inSimplex, outSimplex1, outSimplex2, surfSimplex );
        
        volumeIn.push_back(   inSimplex );
        volumeOut.push_back( outSimplex1 );
        volumeOut.push_back( outSimplex2 );
        surface.push_back(   surfSimplex );
    }
    // case IV: one node is outside
    else {

        // index of the node which is outside
        const unsigned out = ( flags.test(0) == false ? 0 :
                              (flags.test(1) == false ? 1 : 2 ) );

        USimplex<2>::Type outSimplex, inSimplex1, inSimplex2;
        USimplex<1>::Type surfSimplex;
        detail_::cutTriangle<SDIM>( out, indexSimplex, distances, nodes, uniqueNodes,
                                    outSimplex, inSimplex1, inSimplex2, surfSimplex );
        
        volumeIn.push_back(  inSimplex1 );
        volumeIn.push_back(  inSimplex2 );
        
        volumeOut.push_back( outSimplex );

        // reverse surface simplex
        detail_::reverseSimplex<1>( surfSimplex );
        surface.push_back( surfSimplex );

    }

    return;
}
