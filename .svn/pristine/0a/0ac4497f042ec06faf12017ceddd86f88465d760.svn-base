//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   fe/Element.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_fe_element_hpp
#define base_fe_element_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/utility.hpp>
#include <boost/array.hpp>
// base includes
#include <base/verify.hpp>
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace fe{

        template<typename COEFF, typename SFUN>
        class Element;
    }
}

//------------------------------------------------------------------------------
/** Finite element for a generic field.
 *  Considering a generic field \f$ f \f$, scalar- or tensor-valued, and a
 *  function of the \f$ d \f$-dimensional coordinate
 *  \f$ x \in \Omega \subset R^d\f$. By \f$ f^h \f$ we denote the Finite %Element
 *  representation of this field. This can be either an interpolation,
 *  a given representation of e.g. the geometry or the approximation of the
 *  solution field. In any case, \f$ f^h \f$ is represented by its piece-wise
 *  interpolation over elements \f$ \tau_e \f$  which form a decomposition of
 *  \f$ \Omega \f$
 *  \f[
 *         \Omega \approx \Omega^h = \cup_{e} \tau_e
 *  \f]
 *  The restriction of \f$ f^h \f$ to an element \f$ \tau_e \f$ reads
 *  \f[
 *         f^h (x)|_{\tau_e} = \sum_{K} f^K \phi^K(x)
 *  \f]
 *  Here, \f$ f^K \f$ are the (generalised) nodal values of \f$ f^h \f$ and
 *  \f$ \phi^K \f$ are the interpolation functions.
 *  
 *  \tparam COEFF  Type of object holding the field coefficients
 *                 (the \f$ f^K \f$ are held by OBJ)
 *  \tparam SFUN   Type of shape function used for field approximation
 *                 (the \f$ \phi^K \f$ are of type SFUN)
 */
template<typename COEFF, typename SFUN>
class base::fe::Element : boost::noncopyable
{
public:
    //! @name Template parameters
    //@{
    typedef COEFF         Coefficient;
    typedef SFUN          ShapeFun;
    //@}

    //! Geometric shape of the element
    static const base::Shape shape = ShapeFun::shape;

    //! Deduce number of coefficients from shape function
    static const unsigned numCoeff = ShapeFun::numFun;

    //! Manifold dimension determined by the shape 
    static const unsigned dim = base::ShapeDim<shape>::value;

    //! Coefficient storage container
    typedef boost::array<Coefficient*, numCoeff> CoeffPtrArray;
    
    //! Default constructor
    Element()
        : shapeFun_( ShapeFun() )
    {
        coefficients_.assign( NULL );
    }

    //! @name Iterators for coefficient access
    //@{
    typedef typename CoeffPtrArray::iterator       CoeffPtrIter;
    typedef typename CoeffPtrArray::const_iterator CoeffPtrConstIter;

    CoeffPtrIter      coefficientsBegin()       { return coefficients_.begin(); }
    CoeffPtrIter      coefficientsEnd()         { return coefficients_.end();   }
    CoeffPtrConstIter coefficientsBegin() const { return coefficients_.begin(); }
    CoeffPtrConstIter coefficientsEnd()   const { return coefficients_.end();   }
    //@}

    //! Random access to coefficient
    Coefficient* coefficientPtr( const unsigned which ) const
    {
        return coefficients_[ which ];
    }

    //! Access to shape function
    ShapeFun & shapeFun() { return shapeFun_; }

    //! const Access to shape function
    const ShapeFun & shapeFun() const { return shapeFun_; }

    //! @name ID handling
    //@{
    void   setID( const std::size_t id ) { id_ = id; }
    std::size_t getID( ) const { return id_; }
    //@}

    //--------------------------------------------------------------------------
    /** Make a deep copy of a given element.
     *  Pass private data of given element to this one and set
     *  the coefficient pointers of this element to the right location
     *  provided by an iterator pointing to the begin of the new coefficients.
     *  \tparam ELEMENT  Type of element to copy from
     *  \tparam CITER    Type of iterator over coefficients
     *  \param[in] other        Element to copy from
     *  \param[in] coeffIter    Iterator pointing to the coefficients array
     */
    template<typename ELEMENT,typename CITER>
    void deepCopy( const ELEMENT* other,
                   const CITER coeffIter )
    {
        // copy the element ID
        id_ = other -> getID();

        // copy shape function ??

        // pass pointers to coefficients (e.g. nodes)
        CoeffPtrIter cIter = coefficients_.begin();
        CoeffPtrIter cEnd  = coefficients_.end();
        typename ELEMENT::CoeffPtrConstIter cOther = other -> coefficientsBegin();
        for ( ; cIter != cEnd; ++cIter, ++cOther ) {

            // get pointer to the right position in the coefficients container
            CITER thisIter = coeffIter;
            std::advance( thisIter, (*cOther) -> getID() );

            *cIter = *thisIter;
        }
    }
    
private:
    ShapeFun        shapeFun_;     //!< Finite element shape function
    CoeffPtrArray   coefficients_; //!< Field representation coefficients
    std::size_t     id_;           //!< Global ID of this element
};

#endif
