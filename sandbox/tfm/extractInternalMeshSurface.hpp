//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   extractInternalMeshSurface.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef tfm_extractinternalmeshsurface_hpp
#define tfm_extractinternalmeshsurface_hpp

//------------------------------------------------------------------------------
// std  includes
#include <map>
#include <vector>
#include <iterator>
// boost includes
#include <boost/utility.hpp>
// base includes
#include <base/shape.hpp>
#include <base/types.hpp>
// base/aux
#include <base/auxi/SortArray.hpp>
// base/mesh includes
#include <base/mesh/ElementFaces.hpp>
#include <base/mesh/FaceIterator.hpp>

//------------------------------------------------------------------------------
namespace tfm{

    // by coordinates 
    template<typename EITER, typename FILTER>
    void extractInternalMeshSurface(
        EITER first, EITER last, FILTER filter, 
        std::vector< std::pair<std::size_t,unsigned> >&
        boundaryElementContainer );

    // by indices
    template<typename EITER, typename FILTER>
    void extractInternalMeshSurface2(
        EITER first, EITER last, FILTER filter, 
        std::vector< std::pair<std::size_t,unsigned> >&
        boundaryElementContainer );
}

//------------------------------------------------------------------------------
/**
 *  \tparam EITER                Type of element iterator
 *  \tparam FILTER 
 *  \param[in]  first, last      Range of elements
 *  \param[out] boundaryElements Container of (elemNum,faceNum) pairs
 */
template<typename EITER, typename FILTER>
void tfm::extractInternalMeshSurface(
    EITER first, EITER last, FILTER filter,
    std::vector< std::pair<std::size_t,unsigned> >& boundaryElements )
{
    // Template parameter: type of element
    typedef typename base::TypeReduction<
        typename std::iterator_traits<EITER>::value_type>::Type Element;

    // Deduction of geometric entities
    static const base::Shape elemShape = Element::shape;
    static const base::NFace surface   = base::ShapeSurface<elemShape>::value;

    // Iterator over all faces of the element range
    typedef base::mesh::FaceIterator<EITER,surface> FaceIter;

    // Type of face
    typedef typename FaceIter::Face Face;

    // Geometry element is a LagrangeElement
    typedef base::fe::LagrangeElement<Element::shape,
                                      Element::GeomFun::degree> GeomElement;

    // extract face
    typedef base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;

    // Coordinate types
    typedef typename Element::Node::VecDim    GlobalVecDim;
    typedef typename Element::GeomFun::VecDim LocalVecDim;

    // Support points of the geometry shape function
    boost::array<LocalVecDim,Element::GeomFun::numFun> supportPoints;
    Element::GeomFun::supportPoints( supportPoints );

    // Temporary map for storage of the boundary faces
    typedef std::map< Face, std::pair<std::size_t, unsigned> > BoundaryMap;
    BoundaryMap boundaryMap;

    // Go through all faces of the range
    FaceIter faceIter = FaceIter( first );
    FaceIter faceEnd  = FaceIter( last  );

    for ( ; faceIter != faceEnd; ++faceIter ) {

        // number of face
        const unsigned nf = faceIter.faceNum();

        // Get indicdes from face extraction
        std::vector<unsigned> faceIndices;
        FaceExtraction::apply( nf, faceIndices );


        // Geometry filter
        bool passesFilter = true;
        for ( std::size_t f = 0; f < faceIndices.size(); f++ ) {
            
            // parameter coordinate inherited from the domain element
            const LocalVecDim xi = supportPoints[ faceIndices[f] ];

            const GlobalVecDim x =
                base::Geometry<Element>()( *(faceIter.elementIterator()), xi );
        
            // Filter
            if ( not filter( x ) ) passesFilter = false;
        }

        if ( passesFilter ) {
        
            // Get face be dereferencing the iterator
            Face faceSorted = *faceIter;

            // Sort face array in order to have a unique key
            base::auxi::SortArray<Face>::apply( faceSorted );

            // Try to find sorted face in the map
            typename BoundaryMap::iterator check = boundaryMap.find( faceSorted );

            // If not found, insert new pair of Element pointer and face number
            if ( check == boundaryMap.end() ) {
                const std::size_t elemID = (*(faceIter.elementIterator() )) -> getID();
                boundaryMap.insert( std::make_pair( faceSorted,
                                                    std::make_pair( elemID, nf ) ) );
            }
        }
    }

    // Pass element,faceNumber pairs to storage
    typename BoundaryMap::iterator bMapIter = boundaryMap.begin();
    typename BoundaryMap::iterator bMapEnd  = boundaryMap.end();
    for( ; bMapIter != bMapEnd; ++bMapIter )
        boundaryElements.push_back( bMapIter -> second );
        

    return;
}

//------------------------------------------------------------------------------
/**
 *  \tparam EITER                Type of element iterator
 *  \tparam FILTER 
 *  \param[in]  first, last      Range of elements
 *  \param[out] boundaryElements Container of (elemNum,faceNum) pairs
 */
template<typename EITER, typename FILTER>
void tfm::extractInternalMeshSurface2(
    EITER first, EITER last, FILTER filter,
    std::vector< std::pair<std::size_t,unsigned> >& boundaryElements )
{
    // Template parameter: type of element
    typedef typename base::TypeReduction<
        typename std::iterator_traits<EITER>::value_type>::Type Element;

    // Deduction of geometric entities
    static const base::Shape elemShape = Element::shape;
    static const base::NFace surface   = base::ShapeSurface<elemShape>::value;

    // Iterator over all faces of the element range
    typedef base::mesh::FaceIterator<EITER,surface> FaceIter;

    // Type of face
    typedef typename FaceIter::Face Face;

    // Geometry element is a LagrangeElement
    typedef base::fe::LagrangeElement<Element::shape,
                                      Element::GeomFun::degree> GeomElement;

    // extract face
    typedef base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;

    // Coordinate types
    typedef typename Element::Node::VecDim    GlobalVecDim;
    typedef typename Element::GeomFun::VecDim LocalVecDim;

    // Support points of the geometry shape function
    boost::array<LocalVecDim,Element::GeomFun::numFun> supportPoints;
    Element::GeomFun::supportPoints( supportPoints );

    // Temporary map for storage of the boundary faces
    typedef std::map< Face, std::pair<std::size_t, unsigned> > BoundaryMap;
    BoundaryMap boundaryMap;

    // Go through all faces of the range
    FaceIter faceIter = FaceIter( first );
    FaceIter faceEnd  = FaceIter( last  );

    for ( ; faceIter != faceEnd; ++faceIter ) {

        // number of face
        const unsigned nf = faceIter.faceNum();

        // Get indicdes from face extraction
        std::vector<unsigned> faceIndices;
        FaceExtraction::apply( nf, faceIndices );


        // Geometry filter
        std::vector<std::size_t> face;
        for ( std::size_t f = 0; f < faceIndices.size(); f++ ) {

            const std::size_t nodeID =
                (*(faceIter.elementIterator()) ) -> nodePtr( faceIndices[f] ) -> getID();

            face.push_back( nodeID );
        }

        // Filter
        const bool passesFilter = filter( face );

        if ( passesFilter ) {
        
            // Get face be dereferencing the iterator
            Face faceSorted = *faceIter;

            // Sort face array in order to have a unique key
            base::auxi::SortArray<Face>::apply( faceSorted );

            // Try to find sorted face in the map
            typename BoundaryMap::iterator check = boundaryMap.find( faceSorted );

            // If not found, insert new pair of Element pointer and face number
            if ( check == boundaryMap.end() ) {
                const std::size_t elemID = (*(faceIter.elementIterator() )) -> getID();
                boundaryMap.insert( std::make_pair( faceSorted,
                                                    std::make_pair( elemID, nf ) ) );
            }
        }
    }

    // Pass element,faceNumber pairs to storage
    typename BoundaryMap::iterator bMapIter = boundaryMap.begin();
    typename BoundaryMap::iterator bMapEnd  = boundaryMap.end();
    for( ; bMapIter != bMapEnd; ++bMapIter )
        boundaryElements.push_back( bMapIter -> second );
        

    return;
}

#endif
