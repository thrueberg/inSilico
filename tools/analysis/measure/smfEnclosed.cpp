//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smfEnclosed.cpp
//! @author Thomas Rueberg
//! @date   2014

//------------------------------------------------------------------------------
// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>      
// boost includes
#include <boost/lexical_cast.hpp>
// base includes
#include <base/verify.hpp>
#include <base/numbers.hpp>
#include <base/Unstructured.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/kernel/Measure.hpp>
#include <base/Quadrature.hpp>
// tools includes
#include <tools/converter/smf2xx/SMFHeader.hpp>
#include <tools/converter/smf2xx/Conversion.hpp>
// surf includes
#include <surf/Moments.hpp> 

//------------------------------------------------------------------------------
namespace tools{
    namespace analysis{
        namespace enclosed{
        
            //------------------------------------------------------------------
            template<base::Shape SHAPE,unsigned DEGREE>
            struct Analysis
            {
                static void apply( std::istream& smfIn,
                                   std::ostream& out )
                {
                    // Attributes of the mesh
                    static const base::Shape shape   = SHAPE;
                    static const unsigned degree     = DEGREE;
                    static const unsigned    dim     = base::ShapeDim<shape>::value+1;
                    static const unsigned quadOrder  = 5;

                    VERIFY_MSG( (base::ShapeDim<shape>::value < 3),
                                "Works only for surface meshes" );

                    // Mesh type and object
                    typedef base::Unstructured<shape,degree,dim>  Mesh;
                    Mesh mesh;

                    // SMF input
                    base::io::smf::readMesh( smfIn, mesh );

                    // Trivial field binding
                    typedef base::asmb::FieldBinder<Mesh> Binder;
                    Binder binder( mesh );
                    typedef typename Binder::template TupleBinder<>::Type STB;

                    // Quadrature rule
                    base::Quadrature<quadOrder,shape> quadrature;

                    // check closed-ness of surface
                    const double indicator =
                        surf::isClosed<STB>( quadrature, binder );
                    
                    
                    // surface area
                    double area = 0.;
                    base::asmb::simplyIntegrate<STB>(
                        quadrature, area, binder,
                        base::kernel::Measure<typename STB::Tuple>() );
                    
                    // enclosed volume
                    const double volume =
                        surf::enclosedVolume<STB>( quadrature, binder );

                    std::cout << std::setprecision(12);
                    std::cout << "Is closed?      = " << indicator << " (0=yes)\n"
                              << "Surface area    = " << area << "\n"
                              << "Enclosed volume = " << volume << "\n";

                }
            };

            //------------------------------------------------------------------
            // Exclude certain combinations
            struct IDoNotWork
            {
                static void apply( std::istream& smfIn,
                                   std::ostream& out )
                {
                    VERIFY_MSG( false, "Works only for surface meshes" );

                }
            };

            template<unsigned DEGREE>
            struct Analysis<base::POINT,DEGREE> : public IDoNotWork { };
            template<unsigned DEGREE>
            struct Analysis<base::TET,DEGREE> : public IDoNotWork { };
            template<unsigned DEGREE>
            struct Analysis<base::HEX,DEGREE> : public IDoNotWork { };
            
        }
    }
}
//------------------------------------------------------------------------------
/** Read smf formatted file, compute the bounding box and write to stdout
 */
int main( int argc, char * argv[] )
{
    namespace smf2xx   = tools::converter::smf2xx;
    namespace enclosed = tools::analysis::enclosed;
    
    // Sanity check of the number of input arguments
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0]
                  << " fileIn.smf \n\n";
        return 0;
    }

    // Name of smf input file, its basename and the data output file name
    const std::string smfFileIn  = boost::lexical_cast<std::string>( argv[1] );

    // Element attributes
    base::Shape elementShape;
    unsigned    elementNumPoints;
    
    {
        // extract data from header
        std::ifstream smf( smfFileIn.c_str() );
        smf2xx::readSMFHeader( smf, elementShape, elementNumPoints );
        smf.close();
    }
    
    // Input and output file streams
    std::ifstream smfIn(   smfFileIn.c_str() );

    // Call generic conversion helper
    smf2xx::Conversion< enclosed::Analysis >::apply( elementShape,
                                                     elementNumPoints,
                                                     smfIn, std::cout );

    // Close the streams
    smfIn.close();
    
    return 0;
}
