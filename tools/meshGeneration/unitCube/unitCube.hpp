//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   unitCube.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef tools_meshgeneration_unitcube_h
#define tools_meshgeneration_unitcube_h

// std includes
#include <iostream>
#include <vector>
// boost includes
#include <boost/lexical_cast.hpp>
// tools includes
#include <tools/meshGeneration/meshGeneration.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace meshGeneration{
        namespace unitCube{

            //------------------------------------------------------------------
            // Get dimensions from the user
            bool userInput( const int argc, char* argv[], unsigned& dim,
                            unsigned& n1, unsigned& n2, unsigned& n3,
                            unsigned& e1, unsigned& e2, unsigned& e3 )
            {
                // Usage message
                if ( (argc < 2) or (argc > 4) ) {
                    std::cout << "Usage:  " << argv[0]
                              << " N1  [N2  [N3] ]\n";
                    return false;
                }

                // Induced spatial dimension
                dim = boost::lexical_cast<unsigned>( argc-1 );

                // number of elements per direction
                e1 =          boost::lexical_cast<unsigned>( argv[1] );
                e2 = (dim>1 ? boost::lexical_cast<unsigned>( argv[2] ) : 1);
                e3 = (dim>2 ? boost::lexical_cast<unsigned>( argv[3] ) : 1);

                // number of points per direction
                n1 = e1+1;
                n2 = (dim > 1 ? e2 + 1 : 1);
                n3 = (dim > 2 ? e3 + 1 : 1);

                return true;
            }

            //------------------------------------------------------------------
            // Generate equi-distant node grid
            void generatePoints( const unsigned n1, const unsigned n2, const unsigned n3,
                                 const unsigned e1, const unsigned e2, const unsigned e3,
                                 std::vector<tools::meshGeneration::Point> & points )
            {
                // spacing
                const double h1 = 1.0 / static_cast<double>( e1 );
                const double h2 = 1.0 / static_cast<double>( e2 );
                const double h3 = 1.0 / static_cast<double>( e3 );

                //--------------------------------------------------------------
                // node coordinates
                for ( unsigned i3 = 0; i3 < n3; i3++ ) {
                    for ( unsigned i2 = 0; i2 < n2; i2++ ) {
                        for ( unsigned i1 = 0; i1 < n1; i1++ ) {

                            const double x = h1 * i1;
                            const double y = h2 * i2;
                            const double z = h3 * i3;
                            
                            tools::meshGeneration::Point point = {{ x, y, z }};
                            points.push_back( point );
                        }
                    }
                }
                return;
            }


        } // namespace unitCube
    }
}
#endif

