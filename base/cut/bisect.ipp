//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   bisect.ipp
//! @author Thomas Rueberg
//! @date   2013

//------------------------------------------------------------------------------
/** Cut a line element based on a signed distance function.
 *  Given the signed distances of the vertices of a line element, this function
 *  computes the cut point location by linear interpolation along the line,
 *  generates to line elements (one inside, one outside) and identifies the new
 *  point as the surface portion of the cut.
 *  \param[in]     indexSimplex   Vertex IDs of the line element
 *  \param[in]     distances      Values of the signed distance function
 *  \param[in,out] nodes          Coordinates of the vertices and cut points
 *  \param[in,out] uniqueNodes    Map for avoiding node duplications
 *  \param[out]    surface        Connectivity of the surface elements
 *  \param[out]    volumeIn       Connectivity of the volume inside the domain
 *  \param[out]    volumeOut      Connectivity of the volume outside the domain
 */
void base::cut::bisect( const base::cut::USimplex<1>::Type& indexSimplex,
                        const base::cut::DSimplex<1>::Type&   distances,
                        std::vector<base::Vector<1,double>::Type>& nodes,
                        std::map<base::cut::Edge,unsigned>& uniqueNodes,
                        std::vector<base::cut::USimplex<0>::Type>& surface,
                        std::vector<base::cut::USimplex<1>::Type>& volumeIn, 
                        std::vector<base::cut::USimplex<1>::Type>& volumeOut )
{
    // flags checking the values of the distances
    std::bitset<2> flags;
    for ( unsigned d = 0; d < 2; d++ )
        flags[d] = (distances[d] >= 0.);
    
    // case I: all nodes are inside (>= )
    if (flags.count() == 2 ) {
        volumeIn.push_back( indexSimplex );
    }
    // case II: all nodes are outside ( < )
    else if (flags.count() == 0) {
        volumeOut.push_back( indexSimplex );
    }
    // case III: the line is cut
    else {
        // the nodeal indices
        const unsigned i1 = indexSimplex[0];
        const unsigned i2 = indexSimplex[1];
        // the nodal distances
        const double d1 = distances[0];
        const double d2 = distances[1];
        // perform intersection along the line
        const unsigned numNewNode =
            detail_::intersect<1>( i1, i2, d1, d2, nodes, uniqueNodes );
        
        // create volume simplices
        base::cut::USimplex<1>::Type in  = USimplex<1>::create( i1, numNewNode );
        //                               = {{ i1, numNewNode }};
        base::cut::USimplex<1>::Type out = USimplex<1>::create( numNewNode, i2 );
        //                               = {{ numNewNode, i2 }};
        volumeIn.push_back( in );
        volumeOut.push_back( out );
        // create surface simplex
        base::cut::USimplex<0>::Type surf = USimplex<0>::create( numNewNode );
        //                                = {{ numNewNode }};
        surface.push_back( surf );
    }

    return;
}

