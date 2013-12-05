//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Structured.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_structured_hpp
#define base_structured_hpp

//------------------------------------------------------------------------------
// std includes
#include <iterator>
// base includes
#include <base/MultiIndex.hpp>
// base/mesh includes
#include <base/fe/Field.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename ELEMENT>
        class Structured;

        //----------------------------------------------------------------------
        namespace detail_{

            //! Number of nodes of a b-spline grid
            template<unsigned DIM, typename SFUN1D>
            struct NumBSplineGridNodes
            {
                typedef typename base::MultiIndex<DIM>::Type MIT;

                static const unsigned continuity = SFUN1D::continuity;
                
                static const unsigned defect =
                    SFUN1D::degree - SFUN1D::continuity - 1;

                /** For given degree \f$ p \f$ and regualrity \f$ k \f$, the
                 *  total number of grid nodes (i.e. geometry coefficients) is
                 *  \f[
                 *      (p-r+1) n + r
                 *  \f]
                 *   where \f$ n \f$ is the multi-index of number of elements
                 *   per direction.
                 */
                static MIT apply( const MIT & gridSizes )
                {
                    const MIT result =
                        (defect+1) * gridSizes + continuity + 1;
                    return result;
                }
               
            };

            //! Number of nodes of a Lagrangian grid
            template<unsigned DIM, typename SFUN1D>
            struct NumLagrangianGridNodes
            {
                typedef typename base::MultiIndex<DIM>::Type MIT;

                static MIT apply( const MIT & gridSizes )
                {
                    const MIT result =
                        SFUN1D::degree * gridSizes + 1;
                        
                    return result;
                }
               
            };
            
        }
        //----------------------------------------------------------------------
        
    }
}

//------------------------------------------------------------------------------
/** Representation of a structured grid.
 *  This object provides the interface for a structured grid to the mesh
 *  Container. Allocation is carried out by means of the grid dimensions (i.e.,
 *  number of elements per Cartesian direction). An extra interface is provided
 *  to access node and element pointers via a multi-index.
 *  \tparam ELEMENT Type of element of the grid
 */
template<typename ELEMENT>
class base::mesh::Structured
    : public base::fe::Field<ELEMENT,typename ELEMENT::Node>
{
public:
    //! Template parameter: type of element
    typedef ELEMENT Element;

    //! Deduce type of node
    typedef typename Element::Node Node;

    //! Mesh type
    typedef base::fe::Field<Element,Node> Mesh;

    //! Polynomial degree of the geometry shape function dictates node number
    static const unsigned degree = Element::GeomFun::degree;

    //! Manifold dimension of element shape = spatial dimension of the grid
    static const unsigned dim = base::ShapeDim<Element::shape>::value;

    //! Sanity check concerning the element shape
    STATIC_ASSERT_MSG( (Element::shape == base::HyperCubeShape<dim>::value),
                       "Element does not have the shape of a hyper-cube" );

    typedef typename Element::GeomFun::ShapeFun1D ShapeFun1D;
    
    //! Check if shape function is Lagrangian (assume spline else)
    static const bool isLagrangeFun = boost::is_same< ShapeFun1D,
                                                      base::sfun::Lagrange1D<degree> >::value;

    //! Structure providing the number of grid nodes
    typedef typename
    base::IfElse< isLagrangeFun,
                  detail_::NumLagrangianGridNodes<dim,ShapeFun1D>,
                  detail_::NumBSplineGridNodes<dim,ShapeFun1D> >::Type NumGridNodes;

    //! Multi-index structure for element and node access
    typedef base::MultiIndex<dim>     MultiIndex;
    typedef typename MultiIndex::Type MultiIndexType;

    //! Main allocation call
    void allocate( const MultiIndexType & elementDim )
    {
        // Number of grid nodes
        const MultiIndexType nodeDim = NumGridNodes::apply( elementDim);
        
        const std::size_t nNodes = MultiIndex::length( nodeDim );
        Mesh::addCoefficients_(    nNodes );

        const std::size_t nElements = MultiIndex::length( elementDim );
        Mesh::addElements_( nElements );

        // Store the grid dimensions
        gridSizes_ = elementDim;
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


    //! @name Direct access functions via multi-indices
    //@{
    Node* nodePtr( const MultiIndexType & mi ) const
    {
        const std::size_t linIndex = MultiIndex::unwrap( mi, gridSizes_ + degree);
        return Mesh::coefficientPtr( linIndex );
    }

    Element* elementPtr( const MultiIndexType & mi ) const
    {
        const std::size_t linIndex = MultiIndex::unwrap( mi, gridSizes_ );
        return Mesh::elementPtr( linIndex );
    }

    //! Repeat this, since compilation otherwise fails :(
    Element* elementPtr( const std::size_t& index ) const
    {
        return Mesh::elementPtr( index );
    }
    //@}

    //! Let user know about the grid sizes
    MultiIndexType gridSizes() const { return gridSizes_; }

private:
    //! Dimension of the grid
    MultiIndexType gridSizes_;
    
};

#endif
