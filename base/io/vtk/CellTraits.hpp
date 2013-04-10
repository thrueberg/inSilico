//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   CellTraits.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_vtk_celltraits_hpp
#define base_io_vtk_celltraits_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace vtk{
        
            template<base::Shape SHAPE, unsigned NNODES>
            struct CellType;

            template<base::Shape SHAPE, unsigned NNODES>
            struct CellNumOutputNodes
            {
                static const unsigned value = NNODES;
            };

            //------------------------------------------------------------------
            // LINE ELEMENTS
            template<>
            struct CellType<base::LINE,2> { static const unsigned value =  3; };

            template<>
            struct CellType<base::LINE,3> { static const unsigned value = 21; };

            //------------------------------------------------------------------
            // TRIANGLE ELEMENTS
            template<>
            struct CellType<base::TRI,3> { static const unsigned value =  5; };

            template<>
            struct CellType<base::TRI,6> { static const unsigned value = 22; };

            //------------------------------------------------------------------
            // QUADRILATERAL ELEMENTS
            template<>
            struct CellType<base::QUAD,4> { static const unsigned value =  9; };

            template<>
            struct CellType<base::QUAD,9> { static const unsigned value = 23; };

            template<>
            struct CellNumOutputNodes<base::QUAD,9> { static const unsigned value = 8; };

            //------------------------------------------------------------------
            // TETRAHEDRAL ELEMENTS
            template<>
            struct CellType<base::TET, 4> { static const unsigned value = 10; };

            template<>
            struct CellType<base::TET,10> { static const unsigned value = 24; };

            //------------------------------------------------------------------
            // TETRAHEDRAL ELEMENTS
            template<>
            struct CellType<base::HEX, 8> { static const unsigned value = 12; };
            
            template<>
            struct CellType<base::HEX,27> { static const unsigned value = 25; };

            template<>
            struct CellNumOutputNodes<base::HEX,27> { static const unsigned value = 20; };

        }
    }
}


#endif

