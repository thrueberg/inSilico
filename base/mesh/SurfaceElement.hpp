//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   SurfaceElement.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_surfaceelement_hpp
#define base_mesh_surfaceelement_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/verify.hpp>
#include <base/shape.hpp>
#include <base/linearAlgebra.hpp>
#include <base/LagrangeShapeFun.hpp>
// base/mesh includes
#include <base/mesh/ElementFaces.hpp>
#include <base/mesh/Element.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        template<typename ELEMENT,
                 base::Shape SHAPE = base::FaceShape<ELEMENT::shape>::value>
        class SurfaceElement;


        namespace detail_{

            template<typename SELEMENT, base::Shape shape>
            struct LocalDomainCoordinate
            {
                typedef typename SELEMENT::DomainCoordinate DomainCoordinate;
                typedef typename SELEMENT::GeomFun::VecDim  SurfaceCoordinate;
                
                static DomainCoordinate apply( const SELEMENT* sep, 
                                               const SurfaceCoordinate& eta )
                {
                    DomainCoordinate xi = DomainCoordinate::Constant( 0. );
                    typename SELEMENT::GeomFun::FunArray funValues;
                    (sep -> geomFun()).fun( eta, funValues );

                    typename SELEMENT::ParamConstIter pIter = sep -> parametricBegin();
                    typename SELEMENT::ParamConstIter pEnd  = sep -> parametricEnd();
                    for ( unsigned p = 0; pIter != pEnd; ++pIter, p++ )
                        xi += (*pIter) * funValues[p];

                    return xi;
                }
                
            };

            template<typename SELEMENT>
            struct LocalDomainCoordinate<SELEMENT,base::POINT>
            {
                typedef typename SELEMENT::DomainCoordinate DomainCoordinate;
                typedef typename SELEMENT::GeomFun::VecDim  SurfaceCoordinate;
                
                static DomainCoordinate apply( const SELEMENT* sep, 
                                               const SurfaceCoordinate& eta )
                {
                    typename SELEMENT::ParamConstIter pIter = sep -> parametricBegin();
                    DomainCoordinate xi = *pIter;
                    return xi;
                }
                
            };

                
        }
    }
}

//------------------------------------------------------------------------------
/** Representation of a surface element as a parametric restriction of a
 *  volume element.
 *  In addition to a normal base::mesh::Element, this class stores a pointer to
 *  the domain element it is associated with and the parametric coordinates
 *  of that domain element which define this surface element.
 *
 *  \tparam ELEMENT  Type of volume element
 */
template<typename ELEMENT, base::Shape SHAPE>
class base::mesh::SurfaceElement
    : public base::mesh::Element<
    typename ELEMENT::Node,
    base::LagrangeShapeFun<ELEMENT::GeomFun::degree,SHAPE> >
{
public:
    //! @name Template parameter
    //@{
    typedef ELEMENT DomainElement;
    static const base::Shape shape = SHAPE;
    //@}

    //! Type of nodes is the same as for the volume
    typedef typename DomainElement::Node Node;

    //! Local dimension of the surface
    static const unsigned dim = base::ShapeDim<shape>::value;

    //! Geometry function: Lagrangian with same degree as volume element
    typedef typename base::LagrangeShapeFun<ELEMENT::GeomFun::degree,
                                            shape>  GeomFun;

    //! Inherits from this element
    typedef base::mesh::Element<Node,GeomFun>       BasisElement;

    //! This type
    typedef base::mesh::SurfaceElement<DomainElement,shape> SelfType;

    //! Local coordinates of the volume element
    typedef typename base::Vector<dim+1>::Type  DomainCoordinate;

    //! Storage of the parametric coordinates of the volume element
    typedef boost::array<DomainCoordinate,BasisElement::numNodes> ParamtricArray;

    //! Constructor which invalidates data
    SurfaceElement( )
        : domainElementPtr_( NULL )
    {
        parametric_.assign( base::invalidVector<dim+1>() );
    }

    //! @name Set and get volume element pointer
    //@{
    void setDomainElementPointer( DomainElement* vep )
    {
        domainElementPtr_ = vep;
    }

    DomainElement* getDomainElementPointer() const
    {
        return domainElementPtr_;
    }
    //@}

    //! @name Access to parametric coordinates
    //@{
    typedef typename ParamtricArray::iterator ParamIter;
    ParamIter parametricBegin() { return parametric_.begin(); }
    ParamIter parametricEnd()   { return parametric_.end();   }

    typedef typename ParamtricArray::const_iterator ParamConstIter;
    ParamConstIter parametricBegin() const { return parametric_.begin(); }
    ParamConstIter parametricEnd()   const { return parametric_.end();   }
    //@}
    
    //! Override the ID function
    std::size_t getID() const { return domainElementPtr_ -> getID(); }

    //--------------------------------------------------------------------------
    /** Computation of a domain-element parametric coordinate.
     *  This element stores in addition to BasisElement the values of
     *  domain parametric coordinates of its corresponding domain element.
     *  By evaluating the linear combination
     *  \f[
     *       \xi(\eta) = \sum \phi^K(\eta) \xi^K
     *  \f]
     *  the parametric coordinate of the domain element at a local
     *  coordinate \f$ \eta \f$ of this surface element is computed.
     *  \param[in]  eta Local surface parameter coordinate
     *  \return         Value of the corresponding domain coordinate
     */
    DomainCoordinate
    localDomainCoordinate( const typename GeomFun::VecDim& eta ) const
    {
        // Delegate to a policy in order to avoid function evaluation
        // in case of domain-dim = 1 ---> surface-dim = 0
        return detail_::LocalDomainCoordinate<SelfType,
                                              shape>::apply( this, eta );
    }
        
    /** Make a deep copy of a given surface element.
     *  Delegate call to BasisElement and copy the local data.
     *  \tparam SELEMENT Surface element type to copy from
     *  \tparam CITER    Coefficient iterator for passing the pointers
     *  \param[in] other     Other element to copy from
     *  \param[in] coeffIter Iterator pointing to begin of coefficient array
     */
    template<typename SELEMENT,typename CITER>
    void deepCopy( const SELEMENT* other, const CITER coeffIter )
    {
        // Deep copy of basis element
        BasisElement::template deepCopy<SELEMENT,CITER>( other, coeffIter );

        // Pass pointer to domain element
        domainElementPtr_ = other -> getDomainElementPointer();

        // copy parametric coordinates one-by-one
        std::copy( other -> parametricBegin(), other -> parametricEnd(),
                   parametric_.begin() );
        
    }
    
private:
    DomainElement* domainElementPtr_; //!< Pointer to connected volume element
    ParamtricArray       parametric_; //!< Storage of parametric coordinates
};


#endif
