//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   base/Unstructured.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_unstructured_hpp
#define base_unstructured_hpp

//------------------------------------------------------------------------------
#include <base/shape.hpp>
#include <base/mesh/Node.hpp>
#include <base/mesh/Element.hpp>
#include <base/mesh/Unstructured.hpp>
#include <base/LagrangeShapeFun.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace detail_{
        
        //! Derive type of mesh from shape, shape-function and dimension
        template<base::Shape SHAPE,
                 unsigned    GDEG,
                 unsigned    DIM, 
                 template<unsigned DEG, base::Shape> class SFUN>
        struct UnstructuredTraits
        {
            typedef base::mesh::Node<DIM>                 Node;
            typedef SFUN<GDEG,SHAPE>                      SFun;
            typedef base::mesh::Element<Node,SFun>        Element;
            typedef base::mesh::Unstructured<Element>     Type;
        };

    }

    //--------------------------------------------------------------------------
    /** Convenience class for easy mesh type binding.
     *  For most applications based on an unstructured mesh, this interface
     *  will be sufficient. Moreover, the chosen default values for the template
     *  parameters DIM and SFUN also correspond to standard appliations.
     *  In this simple case, the application only needs the code
     *  \code{.cpp}
     *  const base::Shape shape = someShape; // choose from base/shape.hpp
     *  const unsigned    geomDeg = value;   // Polynomial degree of geometry
     *  typedef base::Unstructured<shape,geomDeg> Mesh;
     *  \endcode
     *  to define the type Mesh.
     *  \tparam SHAPE  Geometric shape of element
     *  \tparam GDEG   Polynomial degree of geometry represenation
     *  \tparam DIM    Dimension of the embedding space
     *  \tparam SFUN   Type of shape function for the geometry representation
     */
    template<base::Shape SHAPE,
             unsigned GDEG,
             unsigned DIM=base::ShapeDim<SHAPE>::value,
             template<unsigned DEG,base::Shape> class SFUN=base::LagrangeShapeFun>
    class Unstructured
        : public detail_::UnstructuredTraits<SHAPE,GDEG,DIM,SFUN>::Type
    {
        // empty, all is inherited
    };

    
}

#endif
