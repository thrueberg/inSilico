#include <boost/test/minimal.hpp>

#include <base/linearAlgebra.hpp>
#include <base/auxi/compareNumbers.hpp>

//! \cond SKIPDOX
template<unsigned DIM>
int testVectorEquality( const double a, const double b, const double c )
{
   typedef typename base::Vector<DIM,double>::Type Vec;
   const Vec A = base::constantVector<DIM>( a );
   const Vec B = base::constantVector<DIM>( b );
   const Vec C = base::constantVector<DIM>( c );

   // A == B shall pass
   BOOST_CHECK(     base::auxi::almostEqualVectors<DIM>( A, B ) );
   // A == C shall not pass
   BOOST_CHECK( not base::auxi::almostEqualVectors<DIM>( A, C ) );
   // A == C shall pass with a higher tolerance
   BOOST_CHECK(     base::auxi::almostEqualVectors<DIM>( A, C, std::sqrt(c-a) ) );

   return 0;
}


int test_main( int, char *[] )            
{
    //--------------------------------------------------------------------------
    // test numbers
    const double a = 1.0;
    const double b = a + std::sqrt( std::numeric_limits<double>::epsilon() );
    const double c = a + std::sqrt(b-a);

    // a == b shall pass
    BOOST_CHECK(     base::auxi::almostEqualNumbers( a, b ) );
    // a == c shall not pass
    BOOST_CHECK( not base::auxi::almostEqualNumbers( a, c ) );
    // a == c shall pass with a higher tolerance
    BOOST_CHECK(     base::auxi::almostEqualNumbers( a, c, std::sqrt(c-a) ) );

    //--------------------------------------------------------------------------
    // test vectors
    testVectorEquality<1>( a, b, c );
    testVectorEquality<2>( a, b, c );
    testVectorEquality<3>( a, b, c );

    return 0;
}

//! \endcond
