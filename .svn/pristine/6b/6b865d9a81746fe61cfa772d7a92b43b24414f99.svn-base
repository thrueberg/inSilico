//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   base/Structured.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_structured_hpp
#define base_structured_hpp

//------------------------------------------------------------------------------
#include <base/shape.hpp>
#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Structured.hpp>
#include <base/BSplineShapeFun.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace detail_{
        
        //! Derive type of mesh from dimensions and geometry degree
        template<unsigned    MANIFOLDDIM,
                 unsigned    GDEG,
                 unsigned    GLOBALDIM>
        struct StructuredTraits
        {
            typedef base::mesh::Node<GLOBALDIM>                 Node;
            typedef base::BSplineShapeFun<MANIFOLDDIM,GDEG>     SFun;
            typedef base::mesh::Element<Node,SFun>              Element;
            typedef base::mesh::Structured<Element>             Type;
        };

    }

    //--------------------------------------------------------------------------
    /** Convenience class for easy structured mesh type binding.
     *  For most applications based on a structured mesh, this interface
     *  will be sufficient and the application only needs the code
     *  \code{.cpp}
     *  const unsigned   dim     = someDim;   // spatial dimension
     *  const unsigned   geomDeg = value;   // Polynomial degree of geometry
     *  typedef base::Structured<dim,geomDeg> Grid;
     *  \endcode
     *  to define the type of the Grid (i.e., structured mesh).
     *  \tparam GLOBALDIM   Spatial (embedding) dimension of the problem
     *  \tparam GDEG        Polynomial degree of geometry represenation
     *  \tparam MANIFOLDDIM Local (manifold) dimension of the problem
     */
    template<unsigned GLOBALDIM,
             unsigned GDEG,
             unsigned MANIFOLDDIM=GLOBALDIM>
    class Structured
        : public detail_::StructuredTraits<GLOBALDIM,GDEG,MANIFOLDDIM>::Type
    {
        // empty, all is inherited
    };

    
}

#endif
