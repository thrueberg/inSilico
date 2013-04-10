//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   mesh/Element.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_mesh_element_hpp
#define base_mesh_element_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/utility.hpp>
#include <boost/array.hpp>
// base includes
#include <base/verify.hpp>
#include <base/shape.hpp>
#include <base/LagrangeShapeFun.hpp>
// base/fe includes
#include <base/fe/Element.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename NODE, typename GEOMFUN,
                 typename BASISELEMENT = base::fe::Element<NODE,GEOMFUN> >
        class Element;

        //! Simple case of a linear Lagrange element
        template<typename NODE, base::Shape SHAPE>
        class LinearElement
            : public Element<NODE, base::LagrangeShapeFun<1,SHAPE> >
        {};

    }
}

//------------------------------------------------------------------------------
/** Element of a mesh.
 *  \tparam NODE         Type of node the elements has pointers to
 *  \tparam GEOMFUN      Type of shape function used for geometry representation
 *  \tparam BASISELEMENT Representation of a field (here geometry) on one element
 */
template<typename NODE,typename GEOMFUN, typename BASISELEMENT>
class base::mesh::Element
    : public BASISELEMENT
{
public:

    //! @name Template parameters
    //@{
    typedef NODE             Node;
    typedef GEOMFUN          GeomFun;
    typedef BASISELEMENT     BasisElement;
    //@}

    //! Deduce number of nodes from geometry function
    static const unsigned numNodes = BasisElement::numCoeff;

    //! @name Iterators for DoF access
    //@{
    typedef typename BasisElement::CoeffPtrIter      NodePtrIter;
    typedef typename BasisElement::CoeffPtrConstIter NodePtrConstIter;

    NodePtrIter      nodesBegin()       { return BasisElement::coefficientsBegin(); }
    NodePtrIter      nodesEnd()         { return BasisElement::coefficientsEnd();   }
    NodePtrConstIter nodesBegin() const { return BasisElement::coefficientsBegin(); }
    NodePtrConstIter nodesEnd()   const { return BasisElement::coefficientsEnd();   }
    //@}

    //! Random access to node pointers
    Node* nodePtr( const unsigned n ) const
    {
        return BasisElement::coefficientPtr( n );
    }

    //! Access to shape function
    const GeomFun & geomFun() const { return BasisElement::shapeFun(); }

};


#endif
