//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   dof/IndexMap.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_indexmap_hpp
#define base_dof_indexmap_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <set>
#include <map>
#include <algorithm>
// boost includes
#include <boost/array.hpp>
#include <boost/utility.hpp>
// base includs
#include <base/types.hpp>
// base mesh includes
#include <base/mesh/FaceIterator.hpp>
#include <base/mesh/Structured.hpp>
// base/fe includes
#include <base/fe/Policies.hpp>
// base/aux includes
#include <base/aux/SortArray.hpp>
// base/dof includes
#include <base/dof/copyConnectivity.hpp>
#include <base/dof/generateDoFIndicesFromFaces.hpp>
#include <base/dof/generateDoFIndicesStructured.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename FEBASIS>
        class IndexMap;

        namespace detail_{

            //------------------------------------------------------------------
            //! Call policy in order to avoid to call
            //! e.g. the generation of Cell dofs for a 2D element
            template<typename FEBASIS,base::NFace NFACE,bool CALLME>
            struct FaceDoFGeneration;

            //! Call Generation of DoFs from face
            template<typename FEBASIS,base::NFace NFACE>
            struct FaceDoFGeneration<FEBASIS,NFACE,true>
            {
                template<typename EITER, typename DOFINDEXARRAY>
                static void apply( EITER first, EITER last,
                                   DOFINDEXARRAY& doFIndexArray,
                                   std::size_t& numDoFs )
                {
                    typedef base::mesh::FaceIterator<EITER,NFACE> FaceIter;

                    base::dof::generateDoFIndicesFromFaces<FEBASIS>( FaceIter( first ),
                                                                     FaceIter( last ),
                                                                     doFIndexArray,
                                                                     numDoFs );
                    return;
                }
            };

            //! Do Nothing
            template<typename FEBASIS,base::NFace NFACE>
            struct FaceDoFGeneration<FEBASIS,NFACE,false>
            {
                template<typename EITER, typename DOFINDEXARRAY>
                static void apply( EITER first, EITER last,
                                   DOFINDEXARRAY& doFIndexArray,
                                   std::size_t dummy )
                {
                    return;
                }
            };

            //! Call above functor for all possible faces
            template<typename FEBASIS, unsigned DIM, bool ISUNSTRUCTURED>
            struct CallFaceDoFGeneration
            {
                template<typename EITER, typename DOFINDEXARRAY>
                static std::size_t apply( EITER first, EITER last,
                                          DOFINDEXARRAY& doFIndexArray )
                {
                    return 0; // do nothing
                }

            };
            
            template<typename FEBASIS, unsigned DIM>
            struct CallFaceDoFGeneration<FEBASIS,DIM,true>
            {
                template<typename EITER, typename DOFINDEXARRAY>
                static std::size_t apply( EITER first, EITER last,
                                          DOFINDEXARRAY& doFIndexArray )
                {
                    std::size_t numDoFs = 0;

                    FaceDoFGeneration<FEBASIS,base::VERTEX,
                                      true>::apply( first, last,
                                                    doFIndexArray, numDoFs);
                    
                    FaceDoFGeneration<FEBASIS, base::EDGE,
                                      (DIM >= 1)>::apply( first, last,
                                                          doFIndexArray, numDoFs );
                    
                    FaceDoFGeneration<FEBASIS,base::FACE,
                                      (DIM >= 2)>::apply( first, last,
                                                          doFIndexArray, numDoFs );
                    
                    FaceDoFGeneration<FEBASIS,base::CELL,
                                      (DIM >= 3)>::apply( first, last,
                                                          doFIndexArray, numDoFs );
                    
                    return numDoFs;
                }
            };

            //------------------------------------------------------------------
            template<typename FEBASIS, bool ISSTRUCTURED>
            struct CallStructuredDoFIndexGeneration
            {
                template<typename GRID, typename DOFINDEXARRAY>
                static std::size_t apply( const GRID& grid,
                                          DOFINDEXARRAY& doFIndexArray )
                {
                    return 0; // do nothing
                }
            };

            template<typename FEBASIS>
            struct CallStructuredDoFIndexGeneration<FEBASIS,true>
            {
                template<typename GRID, typename DOFINDEXARRAY>
                static std::size_t apply( const GRID& grid,
                                          DOFINDEXARRAY& doFIndexArray )
                {
                    return 
                        base::dof::generateDoFIndicesStructured<FEBASIS,
                                                                GRID::dim>( doFIndexArray,
                                                                            grid.gridSizes() );
                }
            };

        } // namesapce detail_


    }
}

//------------------------------------------------------------------------------
/** Generation and storage of a look-up table for DoF object indices.
 *  For every element in a given range, this container generates unique index
 *  numbers corresponding to the element's n-faces (i.e., VERTEX, EDGE, .. ).
 *  In the generation process, temporary maps are used in order to guarantee the
 *  uniqueness for the n-Faces which are shared among elements.
 *  \tparam FEBASIS  Description of the used finite element basis
 */
template<typename FEBASIS>
class base::dof::IndexMap
    : boost::noncopyable
{
public:
    //! Template parameter: the FE Basis struct
    typedef FEBASIS FEBasis;

    //! Type of finite element
    typedef typename FEBasis::FiniteElement FiniteElement;

    //! Is continuous (better, maybe, degree of continuity)
    static const bool isContinuous = (FEBasis::continuity >= 0 );

    //! Total number of dofs per element
    static const unsigned numTotalDoFs   = FiniteElement::numTotalDoFs;

    //! Storage of DoF indices per element
    typedef boost::array<std::size_t, numTotalDoFs > ElementIndexArray;

    //--------------------------------------------------------------------------
    //! Constructor to initialise the number of dofs stored
    IndexMap()
        : numDoFs_( 0 )
    {
        // empty
    }

    //--------------------------------------------------------------------------
    //! Access total number of dofs
    std::size_t numDoFs() const { return numDoFs_; }

    //--------------------------------------------------------------------------
    //! Degree of freedom lookup for a given element
    template<typename OUTITER>
    void lookUpElementDoFIndices( OUTITER iter,
                                  const std::size_t elemID ) const
    {
        // Copy the DoF pointers to the provided iterator
        std::copy( elementDoFIndices_[elemID].begin(),
                   elementDoFIndices_[elemID].end(),
                   iter );
    }

    //--------------------------------------------------------------------------
    /** Generation of degrees of freedom indices for a given mesh.
     *  This function generates in the isoparametric case a DoF index for every
     *  node of the mesh or delegates the task to generateDoFIndicesFromFaces,
     *  which does the job for every considered face in the mesh.
     *  Note that this function call is indirected in order to avoid, e.g., 
     *  the generation of cell-dof indices for 2D elements. Only if the
     *  co-dimension between the element shape and the face is larger than or
     *  equal to zero, the other member function is called.
     *  \tparam    MESH   Type of mesh
     *  \param[in] mesh   Mesh 
     */
    template<typename MESH>
    void generateDoFIndices( const MESH& mesh )
    {
        // Number of new elements
        const std::size_t numNewElements = std::distance( mesh.elementsBegin(),
                                                          mesh.elementsEnd() );

        // Current contents
        const std::size_t numCurrentElements = elementDoFIndices_.size();

        // resize storage
        {
            // invalid element
            ElementIndexArray invalid; invalid.assign( base::invalidInt );

            elementDoFIndices_.resize( numCurrentElements + numNewElements,
                                       invalid );
        }

        // deduce element type
        typedef typename MESH::Element Element;
        
        // check for iso-parametric case
        const bool isoParametric =
            ( boost::is_same<typename FEBasis::FEFun,typename Element::GeomFun>::value
              and isContinuous );


        // Check if MESH is base::mesh::Structured
        const bool isStructured =
            boost::is_same< typename base::TypeReduction<MESH>::Type,
                            base::mesh::Structured<Element> >::value;
              

        if ( isoParametric ) {
            // simply copy the node IDs
            numDoFs_ +=
                detail_::copyConnectivity( mesh, elementDoFIndices_  );
        }
        else if ( isStructured ) {
            // gnerate via multi-indices
            numDoFs_ += 
                detail_::CallStructuredDoFIndexGeneration<
                    FEBasis, isStructured>::apply( mesh, elementDoFIndices_ );
        }
        else {
            // generate from faces
            typedef typename MESH::ElementPtrConstIter EIter;
            EIter first = mesh.elementsBegin();
            EIter  last = mesh.elementsEnd();
            
            // use spatial dimension for check of co-dimension
            static const unsigned dim = base::ShapeDim<Element::shape>::value;

            // Compile-time conditional call of face DoF generation
            numDoFs_ += detail_::CallFaceDoFGeneration<
                FEBasis, dim, not isStructured>::apply( first, last,
                                                        elementDoFIndices_ );

        }
    }

public:
    //--------------------------------------------------------------------------
    /** Generate the connectivity pattern among DoF objects.
     *  The purpose of this function is only for illustrations. It generates
     *  for every DoF-ID a vector containing the DoF-IDs of the directly
     *  connected DoFs. A connection is assumed if DoFs are located on the same
     *  element, because in that case their shape functions have an overlapping
     *  support.
     *  \param[in] sparsity  Vector of ID-vectors representing the pattern
     */
    void generateSparsityPattern( std::vector<std::vector<std::size_t> >&
                                  sparsity ) const
    {
        // for every dof create a set of dofs sharing the same element
        std::vector< std::set< std::size_t > > doFConnectivity( numDoFs_ );
        {
            for ( std::size_t e = 0; e < elementDoFIndices_.size(); e++ ) {
                for ( unsigned d1 = 0; d1 < numTotalDoFs; d1++ ) {
                    for ( unsigned d2 = 0; d2 < numTotalDoFs; d2++ ) {
                        
                        const std::size_t thisDoF      = elementDoFIndices_[e][d1];
                        const std::size_t connectedDoF = elementDoFIndices_[e][d2];
                        doFConnectivity[ thisDoF ].insert( connectedDoF );
                        
                    }
                }
            }
        }

        // copy back to input container
        sparsity.resize( doFConnectivity.size() );
        for ( std::size_t d = 0; d < doFConnectivity.size(); d++ ) {
            sparsity[d].resize( doFConnectivity[d].size() );
            std::copy( doFConnectivity[d].begin(), doFConnectivity[d].end(),
                       sparsity[d].begin() );

        }
    }
    

private:
    //! Store for every element, an array of its degree-of-freedom indices
    std::vector<ElementIndexArray> elementDoFIndices_; 

    //! Counter for numbering
    std::size_t numDoFs_;
};

#endif
