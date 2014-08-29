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
// base includes
#include <base/shape.hpp>
#include <base/cut/DecomposeHyperCube.hpp>
#include <base/mesh/HierarchicOrder.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace meshGeneration{
        namespace unitCube{

            bool userInput( const int argc, char* argv[], unsigned& dim,
                            unsigned& e1, unsigned& e2, unsigned& e3 );

            void generateGrid( const unsigned n1, const unsigned n2, const unsigned n3,
                               const double   h1, const double   h2, const double   h3,
                               std::vector<tools::meshGeneration::Point> & points );

            template<unsigned DIM,unsigned DEGREE> struct TriangulateCube;
            template<unsigned DIM,unsigned DEGREE> struct HierarchicCube;

            template<unsigned DIM, bool MAKESIMPLEX, unsigned DEGREE>
            struct Cube
                : base::IfElse<MAKESIMPLEX,
                               TriangulateCube<DIM,DEGREE>,
                               HierarchicCube<DIM,DEGREE> >::Type
            {};

            template<unsigned DIM, bool MAKESIMPLEX, unsigned DEGREE=1>
            struct SMF;

            template<bool MAKESIMPLEX,unsigned DEGREE>
            void unitCubeSMF( const unsigned dim, 
                              const unsigned e1, const unsigned e2, const unsigned e3,
                              std::ostream& out );
                              
        }
    }
}

//------------------------------------------------------------------------------
// Get dimensions from the user
bool tools::meshGeneration::unitCube::userInput(
    const int argc, char* argv[], unsigned& dim,
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

    return true;
}

//------------------------------------------------------------------
// Generate equi-distant node grid
void tools::meshGeneration::unitCube::generateGrid(
    const unsigned n1, const unsigned n2, const unsigned n3,
    const double   h1, const double   h2, const double   h3,
    std::vector<tools::meshGeneration::Point> & points )
{
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

//------------------------------------------------------------------------------
//! Make mesh of simplex elements in cube
template<unsigned DIM,unsigned DEGREE>
struct tools::meshGeneration::unitCube::TriangulateCube
{
    STATIC_ASSERT_MSG( DEGREE==1, "Higher-order simplex not implemented" );
    
    typedef base::cut::DecomposeHyperCube<DIM> DHC;

    typedef base::mesh::HierarchicOrder<base::HyperCubeShape<DIM>::value,DEGREE>
    HO;

    static const unsigned numNodes = base::Binomial<DEGREE+DIM,DEGREE>::value ;
        
    static void apply( const tools::meshGeneration::Element& cube,
                       std::vector<tools::meshGeneration::Element>& elements )
    {
        for ( unsigned s = 0; s < DHC::numSimplices; s++ ) {
            tools::meshGeneration::Element simplex;
            for ( unsigned v = 0; v < DIM+1; v++ ) {
                const std::size_t index = HO::apply( DHC::apply( s, v ) );
                simplex.push_back( cube[index] );
            } 
            elements.push_back( simplex );
        }
        return;
    }
};

//------------------------------------------------------------------------------
//! Generate hierarchic element from tensor-product ordering
template<unsigned DIM,unsigned DEGREE>
struct tools::meshGeneration::unitCube::HierarchicCube
{
    typedef base::mesh::HierarchicOrder<base::HyperCubeShape<DIM>::value,DEGREE>
    HO;

    static const unsigned numNodes = HO::numNodes;
        
    static void apply( const tools::meshGeneration::Element& cube,
                       std::vector<tools::meshGeneration::Element>& elements )
    {
        tools::meshGeneration::Element orderedCube;
        orderedCube.resize( numNodes );
        for ( unsigned v = 0; v < numNodes; v++ ) {
            orderedCube[ HO::apply( v ) ] = cube[v];
        }
        elements.push_back( orderedCube );
            
        return;
    }
};

           
//------------------------------------------------------------------------------
// Generate an SMF file on a unit-cube
template<unsigned DIM, bool MAKESIMPLEX, unsigned DEGREE>
struct tools::meshGeneration::unitCube::SMF
{
    // shape of the elements (static_casts are owed to a INTEL compiler bug)
    static const base::Shape shape = (MAKESIMPLEX ?
                                      static_cast<base::Shape>(
                                          base::SimplexShape<DIM>::value) :
                                      static_cast<base::Shape>(
                                          base::HyperCubeShape<DIM>::value) );

    typedef unitCube::Cube<DIM,MAKESIMPLEX,DEGREE> Element;
    
    // number of points per element
    static const unsigned elementNumPoints = Element::numNodes;

    // call with mesh parameters and stream
    static void apply( const unsigned e1, const unsigned e2, const unsigned e3, 
                       std::ostream& out )
    {
        // name of the element shape
        const std::string shapeName = base::ShapeName<shape>::apply();

        // number of points per direction
        const unsigned n1 = DEGREE*e1 + 1;
        const unsigned n2 = (DIM > 1 ? DEGREE*e2 + 1 : 1);
        const unsigned n3 = (DIM > 2 ? DEGREE*e3 + 1 : 1);
    
        // number of nodes of the mesh
        const unsigned numNodes = n1 * n2 * n3;

        // number of elements of the mesh
        const unsigned numElements =
            e1 * e2 * e3 * (MAKESIMPLEX ? (DIM==3 ? 6 : DIM) : 1);

        // write the header
        tools::meshGeneration::writeSMFHeader( shapeName,
                                               elementNumPoints,
                                               numNodes, numElements, out );

        // spacing
        const double h1 = 1.0 / static_cast<double>( DEGREE * e1 );
        const double h2 = 1.0 / static_cast<double>( DEGREE * e2 );
        const double h3 = 1.0 / static_cast<double>( DEGREE * e3 );
        
        // generate points and write them
        std::vector<tools::meshGeneration::Point> points;
        unitCube::generateGrid( n1, n2, n3, h1, h2, h3, points );
        tools::meshGeneration::writePoints( points, out );

        //--------------------------------------------------------------------------
        // connectivity
        std::vector<tools::meshGeneration::Element> elements;
    
        for ( unsigned i3 = 0; i3 < e3; i3++ ) {
            for ( unsigned i2 = 0; i2 < e2; i2++ ) {
                for ( unsigned i1 = 0; i1 < e1; i1++ ) {

                    tools::meshGeneration::Element cube;
                
                    // serialise index of lower left corner
                    const unsigned i = DEGREE * i1 +
                        (DIM > 1 ? DEGREE * i2 * n1      : 0) +
                        (DIM > 2 ? DEGREE * i3 * n1 * n2 : 0 );

                    // line element
                    if ( DIM == 1 ) {
                        for ( unsigned d = 0; d <= DEGREE; d++ )
                            cube.push_back( i + d );
                    }
                    else if (DIM==2) {
                        // quad
                        for ( unsigned d2 = 0; d2 <= DEGREE; d2++ ){
                            for ( unsigned d1 = 0; d1 <= DEGREE; d1++ ){
                                cube.push_back( i + d2*n1 + d1 );
                            }
                        }
                    }
                    else {
                        // hex
                        for ( unsigned d3 = 0; d3 <= DEGREE; d3++ ){
                            for ( unsigned d2 = 0; d2 <= DEGREE; d2++ ){
                                for ( unsigned d1 = 0; d1 <= DEGREE; d1++ ){
                                    cube.push_back( i + d3*(n1*n2) + d2*n1 + d1 );
                                }
                            }
                        }
                    }

                    // generate elements
                    Element::apply( cube, elements );
                    
                }
            }
        }

        // write to stream
        tools::meshGeneration::writeElements( elements, out );

        return;
    }

};

//------------------------------------------------------------------------------
// Convenience function to call the instantiated templates
template<bool MAKESIMPLEX, unsigned DEGREE>
void tools::meshGeneration::unitCube::unitCubeSMF(
    const unsigned dim, 
    const unsigned e1, const unsigned e2, const unsigned e3,
    std::ostream& out )
{
    namespace unitCube = tools::meshGeneration::unitCube;

    if ( dim == 1 ) {
        unitCube::SMF<1,MAKESIMPLEX,DEGREE>::apply( e1, e2, e3, std::cout );
    }
    else if ( dim == 2 ) {
        unitCube::SMF<2,MAKESIMPLEX,DEGREE>::apply( e1, e2, e3, std::cout );
    }
    else if ( dim == 3 ) {
        unitCube::SMF<3,MAKESIMPLEX,DEGREE>::apply( e1, e2, e3, std::cout );
    }
    

    return;
}

#endif

