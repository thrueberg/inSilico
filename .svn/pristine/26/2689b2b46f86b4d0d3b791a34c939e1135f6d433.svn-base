//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   dof/Field.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_dof_field_hpp
#define base_dof_field_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <algorithm>
// boost includes
#include <boost/utility.hpp>
// base/fe includes
#include <base/fe/Field.hpp>
// base/dof includes
#include <base/dof/Element.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename ELEMENT>
        class Field;

    }
}

//------------------------------------------------------------------------------
/** Representation for FE solution field.
 *  \tparam ELEMENT Type of element representing the DoF interpolation
 */
template<typename ELEMENT>
class base::dof::Field
    : public base::fe::Field<ELEMENT,typename ELEMENT::DegreeOfFreedom>
{
public:
    typedef ELEMENT Element;
    
    //! DoF type defined in basis
    typedef typename Element::DegreeOfFreedom DegreeOfFreedom;

    //! Basis type as container
    typedef typename base::fe::Field<Element,DegreeOfFreedom> Container;
    
    //! @name Container and iterator definitions
    //@{
    typedef typename Container::CoeffPtrIter          DoFPtrIter;
    typedef typename Container::CoeffPtrConstIter     DoFPtrConstIter;
    //@}

    //! @name Iterator access functions
    //@{
    DoFPtrIter      doFsBegin()       { return Container::coefficientsBegin(); }
    DoFPtrIter      doFsEnd()         { return Container::coefficientsEnd();   }
    DoFPtrConstIter doFsBegin() const { return Container::coefficientsBegin(); }
    DoFPtrConstIter doFsEnd()   const { return Container::coefficientsEnd();   }
    //@}

    //! @name Random access
    //@{
    DegreeOfFreedom* doFPtr( const std::size_t& d) const
    {
        return Container::coefficientPtr( d );
    }
    //@}

    //! @name Add storage of new dofs and elements
    //@{
    void addDoFs( const std::size_t numNewDofs )
    {
        Container::addCoefficients_( numNewDofs );
    }

    void addElements( const std::size_t numNewElements )
    {
        Container::addElements_( numNewElements );
    }
    //@}

};

#endif
