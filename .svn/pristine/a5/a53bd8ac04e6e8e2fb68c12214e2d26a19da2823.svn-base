//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   verify.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_verify_hpp
#define base_verify_hpp

//------------------------------------------------------------------------------
// std   includes
#include <cassert>
#include <iostream>
// boost includes
#include <boost/current_function.hpp>
#include <boost/static_assert.hpp>
#include <boost/assert.hpp>

//------------------------------------------------------------------------------
namespace detail_
{
    void assertion_failed( char const * expr, char const * function, 
                           char const * file, long line );

    void assertion_failed_msg( char const * expr, char const * function, 
                               char const * file, long line,
                               char const * message );
} 

//------------------------------------------------------------------------------
/** Macro for permanent assertion.
 *  The macro provided in cassert of the standard library
 *  (see http://www.cplusplus.com/reference/cassert/) is disabled via the
 *  compilation flag `-DNDEBUG` which passes `NDEBUG` as defined to the code.
 *  This flag needs to be set for performance reasons to most external 
 *  libraries, such as Eigen3. But some runtime assertions shall be always
 *  active because they are not time-critical in their evaluation but provide
 *  fundamental insight in case of a runtime error. For this reasone,
 *  this macro is alwasy active and can be used by
 *  \code{.cpp}
 *  VERIFY( expression );
 *  \endcode
 *  In case `expression` evaluates to false, the program aborts with a message
 *  in which file and line the verificiation has failed.
 */
#define VERIFY(expr) \
    ((expr)? ((void)0): \
     ::detail_::assertion_failed(#expr, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__ ))

//------------------------------------------------------------------------------
/** Macro for permanent assertion with user-provided message.
 *  This macro is functionally identical to VERIFY, but has an additional
 *  argument which contains a message from the user of the macro. Therefore,
 *  \code{.cpp}
 *  VERIFY_MSG( expression, message );
 *  \endcode
 *  does the same as `VERIFY( expression )` but in case of a false assertion an
 *  additional message is written to the output.
 */
#define VERIFY_MSG(expr,message) \
    ((expr)? ((void)0): \
     ::detail_::assertion_failed_msg(#expr, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__, message ))

//------------------------------------------------------------------------------
/** Macro for assertion that is just a renaming of the classic c-assert.
 *  This renaming is for beauty reasons only as the upper-case writing clarifies
 *  the macro character of this tool
 *  (see http://www.cplusplus.com/reference/cassert/ for its documentation).
 *  It's usage is clearly
 *  \code{.cpp}
 *  ASSERT( expression );
 *  \endcode
 *  and causes a run-time error if `expression` is false.
 *  \note This macro is disabled with the NDEBUG compilation flag.
 */
#define ASSERT(expr) assert( expr )

//------------------------------------------------------------------------------
/** Macro similar to the ASSERT but with a message.
 *  This tool is just a renaming of boost's assert with message
 *  (see http://www.boost.org/doc/libs/1_54_0/libs/utility/assert.html)
 *  The usage is
 *  \code{.cpp}
 *  ASSERT_MSG( expression, message );
 *  \endcode
 *  and causes a run-time error with the message if the expression is false
 *  \note This macro is disabled with the NDEBUG compilation flag.
 */
#define ASSERT_MSG(expr,message) BOOST_ASSERT_MSG(expr,message)

//------------------------------------------------------------------------------
/** Macro for compile time assertion.
 *  This macro delivers a more verbose output in case of a compilation error.
 *  It is adapted from http://www.boost.org/libs/static_assert/
 *  and used as described therein. 
 */
#ifndef BOOST_NO_STATIC_ASSERT
#  define STATIC_ASSERT_MSG( B, Msg ) static_assert(B, Msg)
#else
#  define STATIC_ASSERT_MSG( B, Msg ) BOOST_STATIC_ASSERT( B )
#endif


//------------------------------------------------------------------------------
/** Output to std::cerr with standard error message
 * 
 *  \param[in]  expr        Boolean expression which is verified to be true.
 *                          Program is stopped if expression is false.
 *  \param[in]  function    Function name in which verification is stated.
 *  \param[in]  file        File in which verification has been called.
 *  \param[in]  line        Line number on which verification is located.
 */
void detail_::assertion_failed( char const * expr, char const * function, 
                                char const * file, long line )
{
    // print error message
    std::cerr << "(EE) Assertion of \"" << expr << "\" in function \"" << function << "\", \n"
              << "(EE) line " << line  << " of file \"" << file << "\" failed! \n";
    abort( );  
}

//------------------------------------------------------------------------------
/** Output to std::cerr with custom error message
 * 
 *  \param[in]  expr        Boolean expression which is verified to be true.
 *                          Program is stopped if expression is false.
 *  \param[in]  function    Function name in which verification is stated.
 *  \param[in]  file        File in which verification has been called.
 *  \param[in]  line        Line number on which verification is located.
 *  \param[in]  message     String containing a message from the programmer
 */
void detail_::assertion_failed_msg( char const * expr, char const * function, 
                                    char const * file, long line,
                                    char const * message )
{
    // print error message
    std::cerr << "(EE) Assertion of \"" << expr << "\" in function \"" << function << "\", \n"
              << "(EE) line " << line  << " of file \"" << file << "\" failed! \n"
              << "     Programmer's message: \n"
              << "     " << message << "\n";
    abort( );  
}

#endif
