//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ExtractMeshFaces.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_extractmeshfaces_hpp
#define base_mesh_extractmeshfaces_hpp

//------------------------------------------------------------------------------
// std  includes
#include <set>
#include <map>
#include <vector>                               
// base includes
#include <base/shape.hpp>
// base/aux
#include <base/aux/SortArray.hpp>
// base/mesh includes
#include <base/mesh/FaceIterator.hpp>
#include <base/mesh/ElementFaces.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<base::Shape SHAPE, base::NFace NFACE>
        class ExtractMeshFaces;
    }
}

//------------------------------------------------------------------------------
/** Function to extract faces from a range of elements
 *  \tparam SHAPE Shape type to extract from
 *  \tparam NFACE Type of face to extract.
 */
template<base::Shape SHAPE,
         base::NFace NFACE = base::ShapeSurface<SHAPE>::value >
class base::mesh::ExtractMeshFaces
{
public:
    //! @name Template parameter
    //@{
    static const base::Shape elemShape = SHAPE;
    static const base::NFace nFace     = NFACE;
    //@}

    //! Type of face
    typedef typename base::mesh::ElementFaces<elemShape,nFace>::Face Face;

    //--------------------------------------------------------------------------
    //! Get unique faces from the range of elements, i.e. the skeleton
    template<typename EITER>
    static void unique( EITER first, EITER last, 
                        std::vector<Face> & faces )
    {
        //! Sanity check
        static const base::Shape iterShape =
            base::aux::TypeReduction<typename EITER::value_type>::Type::shape;
        STATIC_ASSERT_MSG( (iterShape == SHAPE), "Shapes do not match!" );
        
        //! Iterator over the faces
        typedef base::mesh::FaceIterator<EITER,NFACE> FaceIter;

        //! Storage for unique amount (consider boost::flat_set)
        std::set<Face> uniqueSortedFaces;

        //! Begin and end iterators
        FaceIter faceIter = FaceIter( first );
        FaceIter faceEnd  = FaceIter( last  );

        for ( ; faceIter != faceEnd; ++faceIter ) {

            //! Dereference in order to get the face
            Face face = *faceIter;

            //! Use a copy which will be sorted in order to create a UID
            Face sorted = face;
            base::aux::SortArray<Face>::apply( sorted );

            //! Try to insert into set of sorted faces
            std::pair< typename std::set<Face>::iterator,
                       bool> check = uniqueSortedFaces.insert( sorted );

            //! If insertion was successful, a new face has been encountered 
            if ( check.second ) faces.push_back( face );

        }

    }

    //--------------------------------------------------------------------------
    //! Get required faces of every element and store the unqiue ones
    template<typename EITER>
    static void boundary( EITER first, EITER last,
                          std::vector<typename base::mesh::FaceIterator<
                              EITER,NFACE>::Face> & faces )
    {
        //! Sanity checks
        static const base::Shape iterShape =
            base::aux::TypeReduction<typename EITER::value_type>::Type::shape;
        STATIC_ASSERT_MSG( (iterShape == SHAPE), "Shapes do not match!" );
        STATIC_ASSERT_MSG( (base::ShapeDim<SHAPE>::value == NFACE+1),
                           "Co-dimension of boundary is not 1" );
        
        //! Iterator over the faces
        typedef base::mesh::FaceIterator<EITER,NFACE> FaceIter;

        //! Storage for unique amount (consider boost::flat_set)
        //! Key is the sorted face, value is the original face
        typedef std::map< Face, Face > FaceMap;
        FaceMap aux;

        //! Begin and end iterators
        FaceIter faceIter = FaceIter( first );
        FaceIter faceEnd  = FaceIter( last  );

        for ( ; faceIter != faceEnd; ++faceIter ) {

            //! Dereference in order to get the face
            Face face = *faceIter;

            //! Use a copy which will be sorted in order to create a UID
            Face sorted = face;
            base::aux::SortArray<Face>::apply( sorted );

            //! Check if face already exists
            typename FaceMap::iterator check = aux.find( sorted );

            if ( check == aux.end() )
                aux.insert( std::make_pair( sorted, face ) );
            else
                aux.erase( check );
            
        }

        typename FaceMap::iterator auxIter = aux.begin();
        typename FaceMap::iterator auxEnd  = aux.end();
        for ( ; auxIter != auxEnd; ++auxIter ) {

            faces.push_back( (auxIter -> second) );
            
        }

    }

    
    
};


#endif
