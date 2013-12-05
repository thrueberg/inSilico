//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   sgf2xx/Conversion.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef sgf2xx_conversion_hpp
#define sgf2xx_conversion_hpp

//------------------------------------------------------------------------------
// std includes
#include <iostream>
// base includes
#include <base/verify.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace converter{
        namespace sgf2xx{

            template< template<unsigned,unsigned> class CONVERTER>
            struct Conversion;

            namespace detail_{

                template< template<unsigned,unsigned> class CONVERTER,
                          unsigned DIM>
                struct ConversionHelper
                {
                    /** Apply converter object in function of spatial dimension
                     *  \param[in]  degree     Polynomial degree of geometry
                     *  \param[in]  sgf        Input stream (sgf format)
                     *  \param[out] out        Output stream with converted data
                     */
                    static void apply( const unsigned degree,
                                       std::istream& sgf, std::ostream& out )
                    {
                        //! Apply converter depending on the degree
                        switch( degree ) {
                        case 1:
                            CONVERTER<DIM,1>::apply( sgf, out ); break;
                        case 2:
                            CONVERTER<DIM,2>::apply( sgf, out ); break;
                        case 3:
                            CONVERTER<DIM,3>::apply( sgf, out ); break;
                        default:
                            VERIFY_MSG( false,
                                        "Degree not implemented" );
                            
                        }
                        return;
                    }
                };
                
            } // namespace detail_
            
        }
    }
}

//------------------------------------------------------------------------------
/** Generic conversion application in function of spatial dimension and
 *  geometry representation degree.
 *  Given a template type CONVERTER with apply function, it is used for some
 *  conversion from the data in the sgf format to some other format
 *  (implicit in the CONVERTER type).
 *  \tparam CONVERTER   Type of converter
 */
template< template<unsigned,unsigned> class CONVERTER >
struct tools::converter::sgf2xx::Conversion
{
    /** Apply converter object in function of spatial dimension
     *  \param[in]  dim        Spatial dimension
     *  \param[in]  degree     Polynomial degree of geometry
     *  \param[in]  sgf        Input stream (sgf format)
     *  \param[out] out        Output stream with converted data
     */
    static void apply( const unsigned dim, const unsigned degree,
                       std::istream& sgf, std::ostream& out )
    {
        //! Apply converter depending on the dimension
        switch( dim ) {
        
        case 1:
            detail_::ConversionHelper<CONVERTER,1>::apply(degree, sgf, out); break;
        case 2:
            detail_::ConversionHelper<CONVERTER,2>::apply(degree, sgf, out); break;
        case 3:
            detail_::ConversionHelper<CONVERTER,3>::apply(degree, sgf, out); break;
        default:
            VERIFY_MSG( false, "Dimension is wrong" );
        
        }

        return;
    }
};
    
#endif
