//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   meshGeneration.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_meshgeneration_h
#define tools_meshgeneration_h

// std includes
#include <iostream>
#include <vector>
// boost includes
#include <boost/array.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace meshGeneration{

        //----------------------------------------------------------------------
        typedef boost::array<double,3> Point;

        void writePoint( const Point& p, std::ostream& out )
        {
            for ( unsigned d = 0; d < 3; d++ )
                out << p[d] << " ";
            return;
        }


        void writePoints( const std::vector<Point>& points, std::ostream& out )
        {
            std::vector<Point>::const_iterator pIter = points.begin();
            std::vector<Point>::const_iterator pEnd  = points.end();
            for ( ; pIter != pEnd; ++pIter ) {
                writePoint( *pIter, out );
                out << '\n';
            }
    
            return;
        }

        //----------------------------------------------------------------------
        typedef std::vector<std::size_t> Element;

        void writeElement( const Element& e, std::ostream& out )
        {
            for ( unsigned s = 0; s < e.size(); s++ )
                out << e[s] << " ";
            return;
        }

        void writeElements( const std::vector<Element>& elements, std::ostream& out )
        {
            std::vector<Element>::const_iterator eIter = elements.begin();
            std::vector<Element>::const_iterator eEnd  = elements.end();
            for ( ; eIter != eEnd; ++eIter ) {
                writeElement( *eIter, out );
                out << '\n';
            }
    
            return;
        }

        //----------------------------------------------------------------------
        void writeSMFHeader( const std::string elementShape,
                             const unsigned    elementNumPoints,
                             const std::size_t numNodes,
                             const std::size_t numElements,
                             std::ostream& out )
        {
            out << "! elementShape "     << elementShape     << '\n'
                << "! elementNumPoints " << elementNumPoints << '\n'
                << numNodes << "  " << numElements << '\n';
        }

        //----------------------------------------------------------------------
        void writeSMFComment( const std::string& message,
                              std::ostream& out )
        {
            out << "# " << message << '\n';
        }
        
        
    }
}


        
//------------------------------------------------------------------------------
#endif
