//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   createBoundaryFromUnstructured.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_createboundaryfromstructuredface_hpp
#define base_mesh_createboundaryfromstructuredface_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
// boost includes
#include <boost/utility.hpp>
// base includes
#include <base/shape.hpp>
#include <base/types.hpp>
// base/aux
#include <base/aux/SortArray.hpp>
// base/mesh includes
#include <base/mesh/ElementFaces.hpp>
#include <base/mesh/FaceIterator.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        namespace detail_{

            template<unsigned DIM>
            class FaceNumberLookUp
            {
            public:
                static unsigned apply( const unsigned direction,
                                       const unsigned side )
                {
                    return faceNumbers_[ 2*direction + side ];
                }
                
            private:
                // see base/mesh/ElementFaces for explanation
                static const boost::array<unsigned,2*DIM> faceNumbers_;
            };

            // Line element: faces are 0, 1
            template<> const boost::array<unsigned,2>
            FaceNumberLookUp<1>::faceNumbers_ = {{ 0, 1 }};

            // Quad element
            template<> const boost::array<unsigned,4>
            FaceNumberLookUp<2>::faceNumbers_ = {{ 3, 1, 0, 2 }};

            // Hex element
            template<> const boost::array<unsigned,6>
            FaceNumberLookUp<3>::faceNumbers_ = {{ 5, 3, 2, 4, 0, 1 }};
            
        } // namespace detail_

        //----------------------------------------------------------------------
        template<unsigned DIM>
        void createBoundaryFromStructuredFace(
            const typename base::MultiIndex<DIM>::Type& gridSizes,
            const unsigned direction,
            const unsigned side,
            std::vector< std::pair<std::size_t,unsigned> >&
            boundaryElementContainer );
    }
}

//------------------------------------------------------------------------------
/** 
 */
template<unsigned DIM>
void base::mesh::createBoundaryFromStructuredFace(
    const typename base::MultiIndex<DIM>::Type& gridSizes,
    const unsigned direction,
    const unsigned side,
    std::vector< std::pair<std::size_t,unsigned> >&
    boundaryElementContainer )
{
    VERIFY_MSG( (direction < DIM),    "Direction is nonsense ");
    VERIFY_MSG( (side==0) or (side==1), "Side ID is nonsense" );

    // fixed number
    const unsigned fixed = ( side == 0 ? 0 : gridSizes[ direction ]-1 );

    // number of element faces
    typedef typename base::MultiIndex<DIM-1> FaceMI;
    typedef typename FaceMI::Type            FaceMIT;
    FaceMIT numFacesM;

    unsigned ctr = 0;
    for ( unsigned d = 0; d < DIM; d++ )
        if ( d != direction ) numFacesM[ctr++] = gridSizes[d];
    const unsigned numFaces = FaceMI::length( numFacesM );

    // face number
    const unsigned faceNum =
        detail_::FaceNumberLookUp<DIM>::apply( direction, side );

    // go through all faces
    for ( unsigned f = 0; f < numFaces; f ++ ) {

        // multi-index on the grid face
        const FaceMIT fM = FaceMI::wrap( f, numFacesM );

        // construct domain face index
        typename base::MultiIndex<DIM>::Type domIndexM;
        unsigned ctr = 0;
        for ( unsigned d = 0; d < DIM; d++ ) {
            if ( d == direction ) domIndexM[d] = fixed;
            else                  domIndexM[d] = fM[ctr++];
        }

        // linearise
        const std::size_t domIndex = base::MultiIndex<DIM>::unwrap( domIndexM,
                                                                    gridSizes );

        // store pair of element number and face number
        boundaryElementContainer.push_back( std::make_pair( domIndex, faceNum ) );
    }
        
    return;
}

#endif
