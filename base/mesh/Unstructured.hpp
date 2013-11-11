//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   mesh/Unstructured.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_mesh_unstructured_hpp
#define base_mesh_unstructured_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/verify.hpp>
#include <base/types.hpp>
// base/mesh includes
#include <base/fe/Field.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename ELEMENT, bool SHAREDNODES>
        class Unstructured;
        
    }
}

//------------------------------------------------------------------------------
/** Representation of an unstructured mesh.
 *  Makes use of the implementation in base::fe::Field with geometry as the
 *  field and nodes as the coefficients. Provides an interface with number of
 *  nodes and elements for the allocation process.
 *  \tparam ELEMENT Type of element of the mesh
 *  \tparam SHAREDNODES Node-sharing
 */
template<typename ELEMENT,
         bool SHAREDNODES = false>
class base::mesh::Unstructured
    : public base::fe::Field<ELEMENT,typename ELEMENT::Node,SHAREDNODES>
{
public:
    //! Template parameter: type of element
    typedef ELEMENT Element;

    //! Deduce type of node
    typedef typename Element::Node Node;
    
    //! Type of container
    typedef base::fe::Field<Element,Node,SHAREDNODES> Mesh;

    //! Main allocation call
    void allocate( const std::size_t nNodes, const std::size_t nElements )
    {
        Mesh::addCoefficients_( nNodes );
        Mesh::addElements_(     nElements );
    }

    //! @name Iterator access typedefs
    //@{
    typedef typename Mesh::CoeffPtrIter       NodePtrIter;
    typedef typename Mesh::CoeffPtrConstIter  NodePtrConstIter;
    //@}

    //! @name Iterator access functions
    //@{
    NodePtrIter      nodesBegin()        { return Mesh::coefficientsBegin(); }
    NodePtrIter      nodesEnd()          { return Mesh::coefficientsEnd();   }
    NodePtrConstIter nodesBegin()  const { return Mesh::coefficientsBegin(); }
    NodePtrConstIter nodesEnd()    const { return Mesh::coefficientsEnd();   }
    //@}

    //! @name Random access
    //@{
    Node* nodePtr( const std::size_t c) const { return Mesh::coefficientPtr(c); }
    //@}

    //--------------------------------------------------------------------------
    /** Append another mesh to this one.
     *  Given ranges of iterators for nodes and elements, make a deep copy of
     *  all the data. 
     *  \note There will be no check of double nodes. The copied mesh will be
     *        simply added to the containers of this mesh.
     *  \note The mesh to be copied can be destroyed as the copy is truely
     *        \em deep.
     *  \tparam MESH Type of mesh to copy from
     *  \param[in] other The other mesh to copy from
     */
    template<typename NITER, typename EITER>
    void append( NITER firstOtherNode,    NITER lastOtherNode,
                 EITER firstOtherElement, EITER lastOtherElement )
    {
        typedef typename std::iterator_traits<NITER>::value_type OtherNodePtr;
        typedef typename std::iterator_traits<EITER>::value_type OtherElementPtr;

        typedef typename base::TypeReduction<OtherNodePtr>::Type    OtherNode;
        typedef typename base::TypeReduction<OtherElementPtr>::Type OtherElement;
        
        // sanity checks:
        {
            // spatial dimensions
            STATIC_ASSERT_MSG( OtherNode::dim == Node::dim,
                               "Dimensions do not fit" );
            // geometric shapes of elements
            STATIC_ASSERT_MSG( OtherElement::shape == Element::shape,
                               "Shapes do not fit" );
            // geometry representation function of the element
            typedef base::TypeEquality<typename OtherElement::GeomFun,
                                       typename Element::GeomFun> TE;
        }
        
        // sizes of the old mesh
        const std::size_t numOldNodes    = std::distance( nodesBegin(),
                                                          nodesEnd() );
        const std::size_t numOldElements = std::distance( Mesh::elementsBegin(),
                                                          Mesh::elementsEnd() );

        // sizes of the mesh to append
        const std::size_t numNewNodes    = std::distance( firstOtherNode,
                                                          lastOtherNode );
        const std::size_t numNewElements = std::distance( firstOtherElement,
                                                          lastOtherElement );

        // allocate this mesh to hold more nodes
        Mesh::addCoefficients_( numNewNodes );
        Mesh::addElements_(     numNewElements );

        // copy nodes
        NodePtrIter newNodeIter = nodesBegin();
        std::advance( newNodeIter, numOldNodes );
        NITER copyNodeIter = firstOtherNode;
        NITER copyNodeEnd  = lastOtherNode;
        for ( ; copyNodeIter != copyNodeEnd; ++copyNodeIter, ++newNodeIter ) {
            // plain copy
            (*newNodeIter) -> deepCopy( *copyNodeIter );
            // shift ID
            (*newNodeIter) -> setID( (*newNodeIter)->getID() + numOldNodes );
        }

        // Have an iterator point to the begin of the new nodes
        newNodeIter = (this -> nodesBegin() );
        std::advance( newNodeIter, numOldNodes );

        // copy elements
        typename Mesh::ElementPtrIter newElementIter = Mesh::elementsBegin();
        std::advance( newElementIter, numOldElements );
        EITER copyElementIter = firstOtherElement;
        EITER copyElementEnd  = lastOtherElement;
        for ( ; copyElementIter != copyElementEnd;
              ++copyElementIter, ++newElementIter ) {
            // plain copy
            (*newElementIter) -> deepCopy( *copyElementIter, newNodeIter );
            // shift ID
            (*newElementIter) -> setID( (*newElementIter)->getID() +
                                        numOldElements );
        }

        return;
    }

    //--------------------------------------------------------------------------
    /** Append another mesh to this one.
     *  Given a reference to a mesh (with the same geometry representation),
     *  make a deep copy of its data and append them to this mesh.
     *  The implementation is delegated to append<NITER,EITER>() above.
     *  \tparam MESH Type of mesh to copy from
     *  \param[in] other The other mesh to copy from
     */
    template<typename MESH>
    void append( const MESH& other )
    {
        this -> append( other.nodesBegin(),    other.nodesEnd(),
                        other.elementsBegin(), other.elementsEnd() );
    }


};

#endif
