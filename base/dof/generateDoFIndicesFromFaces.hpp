//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   generateDoFIndicesFromFaces.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_generatedofindicesfromfaces_hpp
#define base_dof_generatedofindicesfromfaces_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
// boost includes
#include <boost/array.hpp>
// base includs
#include <base/types.hpp>
// base/mesh includes
#include <base/mesh/FaceIterator.hpp>
// base/fe includes
#include <base/fe/Policies.hpp>
// base/aux includes
#include <base/aux/SortArray.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename FEBASIS,typename FACEITER,
                 typename DOFINDEXARRAY>
        void generateDoFIndicesFromFaces( FACEITER first, FACEITER last,
                                          DOFINDEXARRAY& doFIndexArray,
                                          std::size_t& numDoFs );
    }
}


//--------------------------------------------------------------------------
/** Generation of degree of freedom indices for a range of faces.
 *  Going through a range of faces, this function assigns each of them unique
 *  degree-of-freedom numbers. The unqiueness is optional according to the
 *  choice of the FE Basis (e.g. discontinuous shape functions are possible).
 *  The uniqueness is guaranteed by temporary maps between the already
 *  considered faces and pairs of element ID and local face number.
 *  If continuity is required, this map is queried if the current face has
 *  already been considered. If so, its indices are copied. If not, new
 *  indices are generated.
 *  \note The index container doFIndexArray has to be resized prior to
 *        calling this function.
 *
 *  ### Example ###
 *  As an example for illustration of the steps involved in this function,
 *  consider this mesh
 *  \code{.txt}
 *
 *       6-------7-------8
 *       |       |       |
 *       |  [2]  |  [3]  |
 *       |       |       |
 *       3-------4-------5
 *       |       |       |
 *       |  [0]  |  [1]  |
 *       |       |       |
 *       0-------1-------2
 *
 *  \endcode
 *  which implies the connectivity table
 *
 *    | %Element | Vertices |
 *    | :------: | :------- |
 *    | 0        | 0-1-4-3  |
 *    | 1        | 1-2-5-4  |
 *    | 2        | 3-4-7-6  |
 *    | 3        | 4-5-8-7  |
 *
 *  The task is now to generate the DoF-Indices for fully quadratic, i.e.,
 *  \f$ Q_2 \f$ shape functions. This means that every element has
 *  one DoF per vertex, one DoF per edge and one DoF at the cell center;
 *  in total nine DoFs. The doFIndexArray, passed to this function will
 *  thus be a 2D array of dimensions 4 x 9.
 *  
 *  In the first step, this function is called to generate indices for the
 *  9 vertex DoFs. Finally, this leads to the vertex DoF numbering
 *  \code{.txt}
 *
 *       7-------6-------8
 *       |       |       |
 *       |  [2]  |  [3]  |
 *       |       |       |
 *       3-------2-------5
 *       |       |       |
 *       |  [0]  |  [1]  |
 *       |       |       |
 *       0-------1-------4
 *
 *  \endcode
 *  This numbering comes from going element-wise through the mesh (hidden in the
 *  base::mesh::FaceIterator) and consecutively numbering the newly appearing
 *  vertex DoFs. Therefore, the first element receives new DoF indices 0,1,2,3
 *  and the second element inherits 1,2 and gets the new indices 4 and 5.
 *
 *  The next step is to generate the edge DoFs, beginning with the number 9.
 *  Thus this function is called with FaceIterators for the base::EDGE and the
 *  value of numDoFs is already 9 (result of the previous call). The picture of
 *  the edges DoFs is:
 *  \code{.txt}
 *
 *       X--17---X--20---X
 *       |       |       |
 *      18  [2] 16  [3]  19
 *       |       |       |
 *       X--11---X--15---X
 *       |       |       |
 *      12  [0]  10 [1]  14
 *       |       |       |
 *       X---9---X--13---X
 *
 *  \endcode
 *  Here again, the first element receives four new DoFs, 9, 10, 11, 12, and the
 *  second element inherits the DoF 10 and gets the new ones 13, 14 and 15. For
 *  instance, DoF number 10 lies on the edge connecting vertices 1 and 2. This
 *  edge gets the label [1,2] and the data (0, 1): it is the local edge number
 *  1 of element number 0. Going through the edges of element 1, the local edge
 *  number 3 connects vertices 2 and 1. The ordered (!) label of this edge is
 *  again [1,2]. A lookup in a map reveals, that has already been equipped with
 *  DoF indices as edge 1 of element 0.
 *
 *  Finally, the 9th DoF of every element is associated with the element interior
 *  and therefore disconnected from any neighbor elements. No continuity check is
 *  requrired and every element obtains one additional DoF. This function is
 *  called with a base::mesh::FaceIterator for base::FACE and numDoFs = 20.
 *
 *  The final picture of DoF indices is
 *  \code{.txt}
 *
 *       7--17---6--20---8
 *       |       |       |
 *      18  23  16  24   19
 *       |       |       |
 *       3--11---2--15---5
 *       |       |       |
 *      12  21   10 22   14
 *       |       |       |
 *       0---9---1--13---4
 *
 *  \endcode
 *  and the array doFIndexArray has the entries
 *
 *    | %Element | DoF indices              |
 *    | :------: | :----------------------- |
 *    | 0        | 0-1-2-3;  9-10-11-12; 21 |
 *    | 1        | 1-4-5-2; 13-14-15-10; 22 |
 *    | 2        | 3-2-6-7; 11-16-17-18; 23 |
 *    | 3        | 2-5-8-6; 15-19-20-16; 24 |
 *
 *
 *  \tparam FEBASIS        Type of finite element basis
 *  \tparam FACEITER       Type of face iterator
 *  \tparam DOFINDEXARRAY  Type of array (2D) of elemnet dof indices
 *  \param[in]      first,last    Range of faces
 *  \param[in,out]  doFIndexArray Array with the DoF numbers
 *  \param[in,out]  numDoFs       Total number of unique DoFs
 */
template<typename FEBASIS,typename FACEITER, typename DOFINDEXARRAY>
void base::dof::generateDoFIndicesFromFaces( FACEITER first, FACEITER last,
                                             DOFINDEXARRAY& doFIndexArray,
                                             std::size_t& numDoFs )
{
    // deduce element type
    typedef typename FACEITER::ElementIterator ElementIterator;
    typedef typename
        base::TypeReduction<typename
                            std::iterator_traits<ElementIterator>::value_type>::Type
        Element;

    // Types needed for unique storage
    typedef typename FACEITER::Face Face;
    typedef std::map<Face, std::pair<std::size_t,unsigned> > FaceMap;
    FaceMap     faceMap;

    // Type of face to generate dofs from
    static const base::NFace nFace = FACEITER::nFace;

    // Positions within the element's dof array
    typedef base::fe::DoFArrayPositions<typename FEBASIS::FiniteElement,nFace> DAP;

    // check the co-dimension between face and element
    static const unsigned coDimension =
        base::ShapeDim<Element::shape>::value -
        static_cast<unsigned>( nFace );

    // only do continuity checks if shape functions are continuous
    // and the face has a lower dimension then the element
    const bool continuityCheck =
        ( FEBASIS::continuity > -1 ) and ( coDimension > 0 );

    //----------------------------------------------------------------------
    // Go through all faces
    for( FACEITER fIter = first; fIter != last; ++fIter ) {

        // Element pointer
        Element* elementPtr = *(fIter.elementIterator());

        // Element ID
        const std::size_t elementID = elementPtr -> getID();

        //------------------------------------------------------------------
        // First Section, try to find existing dofs (only for continuous approx)
        if ( continuityCheck ) {

            // Temporary storage of element's sub-faces
            std::vector<Face>        elementFaces;
                
            // Extract vertices first
            typedef base::mesh::ExtractElementVertices<Element> ExtractVertices;
            boost::array<std::size_t,ExtractVertices::numVertices> vertexIDs;
            ExtractVertices::apply( elementPtr, vertexIDs );
            
            // Extract all edges from element
            base::mesh::ExtractElementFaces<Element::shape,nFace>::
                apply( vertexIDs, elementFaces );

            // Go through all faces
            const int begin  = DAP::begin;
            const int stride = DAP::stride;

            // Go through all sub faces
            for ( unsigned f = 0; f < elementFaces.size(); f++ ) {

                // Create a sorted face as unique ID
                Face sortedFace = elementFaces[f];
                base::aux::SortArray<Face>::apply( sortedFace );
            
                // Try to inser new edge into the edge map
                const std::pair<
                    typename FaceMap::iterator,
                    bool>
                    check =
                    faceMap.insert( std::make_pair( sortedFace,
                                                    std::make_pair(elementID,f) ) );

                // Face does already exist in map
                if ( not check.second ) {

                    // ID of that element
                    const std::size_t otherElem = ( (check.first) -> second ).first;

                    // Face number of this face on the other element
                    const unsigned faceNum = ((check.first) -> second).second;
                
                    // Go through all dofs per edge
                    for ( int fDof = 0; fDof < stride; fDof++ ) {
                        doFIndexArray[ elementID     ][ begin + f       * stride + fDof] =
                            doFIndexArray[ otherElem ][ begin + faceNum * stride + fDof];
                    }
                }
                else {

                    // Go through all dofs per edge
                    for ( int fDof = 0; fDof < stride; fDof++ ) {
                        doFIndexArray[ elementID ][ begin + f*stride + fDof ] =
                            numDoFs++;

                    }
                }
                    
                    
            }
                
        }
        else { // No continuity check performed: simply enter new numbers

            // Go through all faces
            const int begin  = DAP::begin;
            const int   end  = DAP::end;

            for ( int fDoF = begin; fDoF < end; fDoF++ ) {
                doFIndexArray[ elementID ][ fDoF ] = numDoFs++;
            }
                
        } // if-else continuity

        
    } // end loop over element range

    return;
}

#endif
