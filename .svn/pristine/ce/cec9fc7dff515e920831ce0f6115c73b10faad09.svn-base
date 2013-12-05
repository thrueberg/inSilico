//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   dof/Element.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_element_hpp
#define base_dof_element_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/utility.hpp>
#include <boost/array.hpp>
// base includes
#include <base/verify.hpp>
#include <base/shape.hpp>
// base/fe includes
#include <base/fe/Element.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename DOF, typename FEFUN,
                 typename BASISELEMENT = base::fe::Element<DOF,FEFUN> >
        class Element;
    }
}

//------------------------------------------------------------------------------
/** Finite element for a field variable.
 *  \tparam DOF          Type of degree of freedom
 *  \tparam FEFUN        Type of shape function used for field approximation
 *  \tparam BASISELEMENT Basis element which holds the field representation for
 *                       one element
 */
template<typename DOF, typename FEFUN, typename BASISELEMENT>
class base::dof::Element : public BASISELEMENT
{
public:
    //! @name Template parameters
    //@{
    typedef DOF           DegreeOfFreedom;
    typedef FEFUN         FEFun;
    typedef BASISELEMENT  BasisElement;
    //@}

    //! Deduce number of nodes from geometry function
    static const unsigned numDoFs = BasisElement::numCoeff;

    //! @name Iterators for DoF access
    //@{
    typedef typename BasisElement::CoeffPtrIter      DoFPtrIter;
    typedef typename BasisElement::CoeffPtrConstIter DoFPtrConstIter;

    DoFPtrIter      doFsBegin()       { return BasisElement::coefficientsBegin(); }
    DoFPtrIter      doFsEnd()         { return BasisElement::coefficientsEnd();   }
    DoFPtrConstIter doFsBegin() const { return BasisElement::coefficientsBegin(); }
    DoFPtrConstIter doFsEnd()   const { return BasisElement::coefficientsEnd();   }
    //@}

    //! Access to shape function
    const FEFun & fEFun() const { return BasisElement::shapeFun(); }
};

#endif
