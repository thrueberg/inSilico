//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Unstructured.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_unstructured_hpp
#define base_unstructured_hpp

//------------------------------------------------------------------------------
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

};

#endif
