//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   HyperCube.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_hypercube_hpp
#define base_cut_hypercube_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/meta.hpp>
#include <base/linearAlgebra.hpp>

namespace base{
    namespace cut{

        namespace detail_{

            // Recursively generate the vertices of a hypercube
            template<unsigned DIM>
            struct GenerateVertices
            {
                static const unsigned size = base::MToTheN<2,DIM>::value;
                
                static boost::array<typename base::Vector<DIM,double>::Type,size> apply()
                {
                    // result array
                    boost::array<typename base::Vector<DIM,double>::Type,size> result;
                    // result from lower dimension
                    boost::array<typename base::Vector<DIM-1,double>::Type,
                                 size/2> lowerDimPoints
                    = GenerateVertices<DIM-1>::apply();

                    // construct results
                    unsigned ctr = 0;
                    for ( unsigned nOuter = 0; nOuter < 2; nOuter++ ) {
                        for ( unsigned nInner = 0; nInner < size/2; nInner ++ ) {

                            // copy from lower dimension
                            for ( unsigned d = 0; d < DIM-1; d++ )
                                result[ctr][d] = lowerDimPoints[nInner][d];
                            // add (0/1) to this dimension's entry
                            result[ctr][DIM-1] = 
                                static_cast<double>( nOuter );
                            ctr++;
                        }
                    }
                    
                    return result;
                }
            };

            // Special case for DIM=1
            template<>
            struct GenerateVertices<1>
            {
                static boost::array<base::Vector<1,double>::Type,2> apply()
                {
                    boost::array<base::Vector<1,double>::Type,2> result;
                    result[0][0] = 0.;
                    result[1][0] = 1.;
                    return result;
                }
            };
            
        }

        //--------------------------------------------------------------------------
        /** Store and provide access to the vertices of a DIM-hypercube.
         *  In lexicographic order, the vertices of a hypercube are:
         *  - DIM = 1
         *    | Number |  Vertex |
         *    |:------:|:-------:|
         *    |   0    |    0    |
         *    |   1    |    1    |
         *
         *  - DIM = 2
         *    | Number |  Vertex |
         *    |:------:|:-------:|
         *    |   0    |  (0,0)  |
         *    |   1    |  (1,0)  |
         *    |   2    |  (0,1)  |
         *    |   3    |  (1,1)  |
         *
         *  - DIM = 3
         *    | Number |  Vertex |
         *    |:------:|:-------:|
         *    |   0    | (0,0,0) |
         *    |   1    | (1,0,0) |
         *    |   2    | (0,1,0) |
         *    |   3    | (1,1,0) |
         *    |   4    | (0,0,1) |
         *    |   5    | (1,0,1) |
         *    |   6    | (0,1,1) |
         *    |   7    | (1,1,1) |
         *
         *  See: http://en.wikipedia.org/wiki/Hypercube
         *  \tparam DIM Spatial dimension of the hypercube
         */
        template<unsigned DIM>
        class HyperCube
        {
        public:
            //! Type of vertex
            typedef typename base::Vector<DIM,double>::Type VecDim;

            //! Number of vertices
            static const unsigned numVertices = base::MToTheN<2,DIM>::value;

            //! Type of array for storage
            typedef boost::array<VecDim,numVertices> VertexArray;

            //! Return the \a i -th vertex
            static VecDim vertex( const unsigned i ) { return vertexArray_[i];}
        private:
            static const VertexArray vertexArray_; //!< The vertices
        };

        //! Initialise via generic helper template
        template<unsigned DIM>
        const typename HyperCube<DIM>::VertexArray HyperCube<DIM>::vertexArray_ =
            detail_::GenerateVertices<DIM>::apply();
    }
}
#endif
