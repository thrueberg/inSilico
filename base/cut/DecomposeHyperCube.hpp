//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   DecomposeHyperCube.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_decomposehypercube_hpp
#define base_cut_decomposehypercube_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/meta.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        //----------------------------------------------------------------------
        /** Decomposition of a DIM-Hypercube into simplex elements.
         *  Consider the task to decompose the hypercube as created in
         *  base::cut::HyperCube into simplex elements. In 1D this task is
         *  trivial and 2D it yields 2 triangles with a cut of the square along
         *  one of the diagonals. Here, the choice is to cut along the 1-2
         *  diagonal. The decomposition becomes non-trivial in 3D where several
         *  variants are available to decompose a cube into tetrahedra. In order
         *  to ensure that the cuts along the cube's faces are coincident for
         *  neighbouring cubes. This is not possible for the minimal
         *  decomposition into 5 tetrahedra, but for 6 it is possible.
         *
         *  \image html cube.png
         *
         *  \tparam DIM Spatial dimension of the hypercube
         */
        template<unsigned DIM>
        class DecomposeHyperCube
        {
        public:
            //! Number of simplices in the decomposition
            static const unsigned numSimplices = base::Factorial<DIM>::value;
            //! A simplex
            typedef boost::array<unsigned,DIM+1>              Simplex;
            //! The decomposition table
            typedef boost::array<Simplex,numSimplices>        SimplexTable;

            //! Return from a specified simplex a certain node number
            static unsigned apply( const unsigned numSimplex,
                                   const unsigned numNode )
            {
                return simplexTable_[numSimplex][numNode];
            }

        private:
            //! Store the decomposition
            static const SimplexTable simplexTable_;
        };

        //----------------------------------------------------------------------
        // Line element has no sub-simplices
        template<>
        const DecomposeHyperCube<1>::SimplexTable
        DecomposeHyperCube<1>::simplexTable_ = {{
                {{ 0, 1 }}
            }};

        //----------------------------------------------------------------------
        // Quadrilateral element has two sub-triangles
        template<>
        const DecomposeHyperCube<2>::SimplexTable
        DecomposeHyperCube<2>::simplexTable_ = {{
                {{ 0, 1, 2 }},
                {{ 3, 2, 1 }}
            }};

        //----------------------------------------------------------------------
        // Hexahedral element has six sub-tetrahedra
        template<>
        const DecomposeHyperCube<3>::SimplexTable
        DecomposeHyperCube<3>::simplexTable_ = {{
                {{ 0, 1, 2, 4 }},  // first  wedge 0-1-2=4-5-6
                {{ 1, 2, 4, 5 }},
                {{ 2, 4, 5, 6 }},
                {{ 1, 2, 5, 3 }},  // second wedge 1-3-2=5-7-6
                {{ 2, 5, 3, 6 }},
                {{ 5, 3, 6, 7 }}
            }};
        
    }
}

#endif
