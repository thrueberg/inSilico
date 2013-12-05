//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   smf2xx/Conversion.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef smf2xx_conversion_hpp
#define smf2xx_conversion_hpp

//------------------------------------------------------------------------------
// std includes
#include <iostream>
// base includes
#include <base/verify.hpp>
#include <base/shape.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace smf2xx{

            template< template<base::Shape,unsigned> class CONVERTER>
            struct Conversion;
        }
    }
}

//------------------------------------------------------------------------------
/** Generic conversion application in function of geometry shape and degree.
 *  Given a template type CONVERTER with apply function, it is used for some
 *  conversion from the data in the smf format to some other format
 *  (implicit in the CONVERTER type).
 *  \tparam CONVERTER   Type of converter
 */
template< template<base::Shape,unsigned> class CONVERTER >
struct tools::converter::smf2xx::Conversion
{
    /** Apply converter object in function of shape and numPoints
     *  \param[in]  shape      Geometric shape of element
     *  \param[in]  numPoints  Number of geometry nodes of element
     *  \param[in]  smf        Input stream (smf format)
     *  \param[out] out        Output stream with converted data
     */
    static void apply( const base::Shape shape, const unsigned numPoints,
                       std::istream& smf, std::ostream& out )
    {
        //! Apply converter depending on the combination (shape,numPoints)
        switch( shape ) {
        
        case base::LINE:
            // Possible LINE element implementations
            switch (numPoints) {
            case 2: CONVERTER<base::LINE,1>::apply( smf, out ); break;
            case 3: CONVERTER<base::LINE,2>::apply( smf, out ); break;
            default:
                VERIFY_MSG( false, "LINE element with this number of points not supported" );
            }
            break;

        case base::TRI:
            // Possible TRI element implementations
            switch (numPoints) {
            case  3: CONVERTER<base::TRI,1>::apply( smf, out ); break;
            case  6: CONVERTER<base::TRI,2>::apply( smf, out ); break;
            default:
                VERIFY_MSG( false, "TRI element with this number of points not supported" );
            }
            break;

        case base::QUAD:
            // Possible QUAD element implementations
            switch (numPoints) {
            case  4: CONVERTER<base::QUAD,1>::apply( smf, out ); break;
            case  9: CONVERTER<base::QUAD,2>::apply( smf, out ); break;
            default:
                VERIFY_MSG( false, "QUAD element with this number of points not supported" );
            }
            break;

        case base::TET:
            // Possible TET element implementations
            switch (numPoints) {
            case  4: CONVERTER<base::TET,1>::apply( smf, out ); break;
            case 10: CONVERTER<base::TET,2>::apply( smf, out ); break;
            default:
                VERIFY_MSG( false, "TET element with this number of points not supported" );
            }
            break;
        
        case base::HEX:
            // Possible HEX element implementations
            switch (numPoints) {
            case  8: CONVERTER<base::HEX,1>::apply( smf, out ); break;
            case 27: CONVERTER<base::HEX,2>::apply( smf, out ); break;
            default:
                VERIFY_MSG( false, "HEX element with this number of points not supported" );
            }
            break;
        
        default:
            VERIFY_MSG( false, "Shape not detected" );
        
        }

    }
};
    
#endif
