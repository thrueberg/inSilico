//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   base/Field.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_field_hpp
#define base_field_hpp

//------------------------------------------------------------------------------
#include <base/shape.hpp>
#include <base/dof/DegreeOfFreedom.hpp>
#include <base/dof/Element.hpp>
#include <base/dof/Field.hpp>


//------------------------------------------------------------------------------
namespace base{

    namespace detail_{

        template<typename FEBASIS, unsigned DOFSIZE, unsigned NHIST>
        struct FieldTraits
        {
            typedef base::dof::DegreeOfFreedom<DOFSIZE,NHIST>    DegreeOfFreedom;
            typedef typename FEBASIS::FEFun                      FieldFun;
            typedef base::dof::Element<DegreeOfFreedom,FieldFun> Element;
            typedef base::dof::Field<Element>                    Type;
        };
    }

    //--------------------------------------------------------------------------
    /** Convenience type for easy field type binding.
     *  In a standard FE Analysis the representation of an FE field, e.g. the
     *  solution field, is uniquely defined by the chosen basis and the type
     *  of degree of freedom. Here, the latter is usually defined in terms of
     *  the size of the individual DoF, i.e. its number of components, and the
     *  size of the history storage as needed in time-dependent problems.
     *  Given the type of FE Basis, the size of the degree of freedom and,
     *  optinally, the number of history terms to store, this object inherits
     *  from the corresponding base::dof::Field class for a convenient type
     *  binding.
     *  \tparam FEBASIS Type of FE Basis
     *  \tparam DOFSIZE Size of individual degree of freedom
     *  \tparam NHIST   Number of history terms to store
     */
    template<typename FEBASIS, unsigned DOFSIZE, unsigned NHIST = 0>
    struct Field
        : public detail_::FieldTraits<FEBASIS,DOFSIZE,NHIST>::Type
    {
        // empty, all is inherited
    };
    
}
#endif
