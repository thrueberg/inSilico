//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ElementFilter.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_io_gp_elementfilter_hpp
#define base_io_gp_elementfilter_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
// boost includes

// base includes
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace io{
        namespace gp{

            //------------------------------------------------------------------
            /** Provide an element connectivity filter for line plots.
             *  Using GnuPlot for a connectivity plot, the element edges shall
             *  be plotted. This filter gives the order in which the points
             *  are to be plotted in order to look like a wireframe of the
             *  element. Therefore, some edges or points need to be detached
             *  from the other lines. This is indicated by inserting '-1' into
             *  the filter list.
             *  \tparam SHAPE
             *  \tparam NUMNODES
             */
            template<base::Shape SHAPE,unsigned NUMNODES>
            class ElementFilter
            {
            public:
                typedef std::vector<int>                      FilterVector;
                typedef typename FilterVector::const_iterator FilterIter;
                
                ElementFilter();

                FilterIter begin() const { return filter_.begin(); }
                FilterIter   end() const { return filter_.end(); }

            private:
                std::vector<int> filter_;
            };

            //------------------------------------------------------------------
            //! \cond SKIPDOX

            //------------------------------------------------------------------
            // SHAPE==LINE
            
            // <POINT,1>
            template<>
            ElementFilter<base::POINT,1>::ElementFilter()
            {
                filter_.push_back( 0 ); 
            }

            //------------------------------------------------------------------
            // SHAPE==LINE
            
            // <LINE,2>
            template<>
            ElementFilter<base::LINE,2>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 1 );
            }

            // <LINE,3>
            template<>
            ElementFilter<base::LINE,3>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 2 );
                filter_.push_back( 1 );
            }
            
            //------------------------------------------------------------------
            // SHAPE==TRI
            
            // <TRI,3>
            template<>
            ElementFilter<base::TRI,3>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 1 );
                filter_.push_back( 2 ); filter_.push_back( 0 );
            }
            
            // <TRI,6>
            template<>
            ElementFilter<base::TRI,6>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 3 );
                filter_.push_back( 1 ); filter_.push_back( 4 );
                filter_.push_back( 2 ); filter_.push_back( 5 );
                filter_.push_back( 0 );
            }

            //------------------------------------------------------------------
            // SHAPE==QUAD
            
            // <QUAD,4>
            template<>
            ElementFilter<base::QUAD,4>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 1 );
                filter_.push_back( 2 ); filter_.push_back( 3 );
                filter_.push_back( 0 );
            }
            
            // <QUAD,9>
            template<>
            ElementFilter<base::QUAD,9>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 4 );
                filter_.push_back( 1 ); filter_.push_back( 5 );
                filter_.push_back( 2 ); filter_.push_back( 6 );
                filter_.push_back( 3 ); filter_.push_back( 7 );
                filter_.push_back( 0 ); filter_.push_back( -1 );
                filter_.push_back( 8 );
            }

            //------------------------------------------------------------------
            // SHAPE==TET
            
            // <TET,4>
            template<>
            ElementFilter<base::TET,4>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 1 );
                filter_.push_back( 2 ); filter_.push_back( 0 );
                filter_.push_back( 3 ); filter_.push_back( 1 );
                filter_.push_back( -1 );
                filter_.push_back( 3 ); filter_.push_back( 2 );
            }
            
            // <TET,10>
            template<>
            ElementFilter<base::TET,10>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 4 );
                filter_.push_back( 1 ); filter_.push_back( 5 );
                filter_.push_back( 2 ); filter_.push_back( 6 );
                filter_.push_back( 0 ); filter_.push_back( 7 );
                filter_.push_back( 3 ); filter_.push_back( 8 );
                filter_.push_back( 1 );
                filter_.push_back(-1 );
                filter_.push_back( 3 ); filter_.push_back( 9 );
                filter_.push_back( 2 );
            }

            //------------------------------------------------------------------
            // SHAPE==QUAD
            
            // <HEX,8>
            template<>
            ElementFilter<base::HEX,8>::ElementFilter()
            {
                filter_.push_back( 0 ); filter_.push_back( 1 );
                filter_.push_back( 2 ); filter_.push_back( 3 );
                filter_.push_back( 0 ); filter_.push_back( 4 );
                filter_.push_back( 5 ); filter_.push_back( 1 );
                filter_.push_back( -1 );
                filter_.push_back( 5 ); filter_.push_back( 6 );
                filter_.push_back( 2 ); filter_.push_back(-1 );
                filter_.push_back( 6 ); filter_.push_back( 7 );
                filter_.push_back( 3 ); filter_.push_back(-1 );
                filter_.push_back( 7 ); filter_.push_back( 4 );
            }
            
            // <HEX,27>
            template<>
            ElementFilter<base::HEX,27>::ElementFilter()
            {
                filter_.push_back(  0 ); filter_.push_back(  8 );
                filter_.push_back(  1 ); filter_.push_back(  9 );
                filter_.push_back(  2 ); filter_.push_back( 10 );
                filter_.push_back(  3 ); filter_.push_back( 11 );
                filter_.push_back(  0 ); filter_.push_back( 16 );
                filter_.push_back(  4 ); filter_.push_back( 12 );
                filter_.push_back(  5 ); filter_.push_back( 17 );
                filter_.push_back( -1 );
                filter_.push_back(  5 ); filter_.push_back( 13 );
                filter_.push_back(  6 ); filter_.push_back( 18 );
                filter_.push_back(  2 ); filter_.push_back( -1 );
                filter_.push_back(  6 ); filter_.push_back( 14 );
                filter_.push_back(  7 ); filter_.push_back( 19 );
                filter_.push_back(  3 ); filter_.push_back( -1 );
                filter_.push_back(  7 ); filter_.push_back( 15 );
                filter_.push_back(  4 ); filter_.push_back( -1 );
                for ( int i = 20; i < 27; i++ ) {
                    filter_.push_back( i ); filter_.push_back( -1 );
                }
            }
            //! \endcond
            //
            
        }
    }
}

#endif

