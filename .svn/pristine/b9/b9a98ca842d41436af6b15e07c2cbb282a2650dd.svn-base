//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ElementTraits.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_xdmf_elementtraits_hpp
#define base_io_xdmf_elementtraits_hpp

//------------------------------------------------------------------------------
// std includes
#include <string>
// base includes
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace xdmf{

            template<base::Shape SHAPE, unsigned NNODES>
            struct ElementName;

            template<base::Shape SHAPE, unsigned NNODES>
            struct ElementNumOutputNodes
            {
                static const unsigned value = NNODES;
            };

            //! \cond
            //------------------------------------------------------------------
            // LINE elements
            template<>
            struct ElementName<base::LINE,2>
            {
                static std::string value() { return "Polyline"; }
            };

            template<>
            struct ElementName<base::LINE,3>
            {
                static std::string value() { return "Edge_3"; }
            };

            //------------------------------------------------------------------
            // TRIANGLE
            template<>
            struct ElementName<base::TRI,3>
            {
                static std::string value() { return "Triangle"; }
            };
            
            template<>
            struct ElementName<base::TRI,6>
            {
                static std::string value() { return "Tri_6"; }
            };

            //------------------------------------------------------------------
            // QUADRILATERAL
            template<>
            struct ElementName<base::QUAD,4>
            {
                static std::string value() { return "Quadrilateral"; }
            };

            template<>
            struct ElementName<base::QUAD,9>
            {
                static std::string value() { return "Quad_8"; }
            };

            template<>
            struct ElementNumOutputNodes<base::QUAD,9>
            {
                static const unsigned value = 8;
            };
            
            //------------------------------------------------------------------
            // TETRAHEDRON
            template<>
            struct ElementName<base::TET,4>
            {
                static std::string value() { return "Tetrahedron"; }
            };

            template<>
            struct ElementName<base::TET,10>
            {
                static std::string value() { return "Tet_10"; } 
            };

            //------------------------------------------------------------------
            // HEXAHEDRON
            template<>
            struct ElementName<base::HEX,8>
            {
                static std::string value() { return "Hexahedron"; }
            };

            template<>
            struct ElementName<base::HEX,27>
            {
                static std::string value() { return "Hex_20"; }
            };

            template<>
            struct ElementNumOutputNodes<base::HEX,27>
            {
                static const unsigned value = 20;                
            };
            //! \endcond
        }
    }
}

#endif
