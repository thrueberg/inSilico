//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   meta.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_meta_hpp
#define base_meta_hpp

// base includes
#include <base/verify.hpp>

//------------------------------------------------------------------------------
namespace base{

    //--------------------------------------------------------------------------
    /**  Compute \f$ M^N \f$.
     *   \tparam M base
     *   \tparam N exponent
     */
    template<unsigned M, unsigned N>
    struct MToTheN
    {
        static const unsigned value = M * MToTheN<M,N-1>::value;
    };
    
    //! \cond SKIPDOX Specialisation for N=0.
    template<unsigned M>
    struct MToTheN<M,0>
    {
        static const unsigned value = 1;
    };
    //! \endcond

    //--------------------------------------------------------------------------
    /** Compute the integral N-th power of floating point number, \f$ x^N \f$.
     *  \tparam N exponent
     */
    template<unsigned N>
    struct Power
    {
        static double apply( const double mantissa )
        {
            return mantissa * Power<N-1>::apply( mantissa );
        }
    };

    //! \cond SKIPDOX Specialisation for N=0
    template<>
    struct Power<0>
    {
        static double apply( const double mantissa ) { return 1.0; }
    };
    //! \endcond

    //--------------------------------------------------------------------------
    /** Compute the factorial of N,  \f$ N! \f$.
     *  A static assertion is triggered when \f$ N \leq 13 \f$ as this would
     *  result beyond the range of ordinary integers.
     *  \sa http://en.wikipedia.org/wiki/Factorial
     *  \tparam N
     */
    template<unsigned N>
    struct Factorial
    {
        STATIC_ASSERT_MSG( (N < 13), "Risk of integer overflow" );   
        static const unsigned value = N * Factorial<N-1>::value;
    };
    
    //! \cond SKIPDOX Specialisation for N=0.
    template<> struct Factorial<0>{ static const unsigned value = 1; };
    //! \endcond

    //--------------------------------------------------------------------------
    /** Compute the binomial coefficient \f$ C(N,K) \f$.
     *  This coefficient is defined for \f$ 0 \leq K \leq K \f$ as
     *  \f[
     *      C(N,K) = \frac{ N! }{ K! (N-K)! }\,.
     *  \f]
     *  Consider Pascal's triangle or 'N choose K'. Here, we extend the
     *  definition to \f$ K > N \f$ by setting the values to zero.
     *  \sa http://en.wikipedia.org/wiki/Binomial_coefficient
     *  \tparam N,K  Non-negative integers 
     */
    template<unsigned N, unsigned K>
    struct Binomial
    {
        // Necessary aux variable in order to avoid the call of Factorial
        // with a negative valued integer converted to a large unsigned 
        static const unsigned diff  = ( K > N ? 0 : N-K );
        static const unsigned value = ( K > N ? 0 :
                                        Factorial<N>::value / 
                                        Factorial<K>::value /
                                        Factorial<diff>::value );
    };

    //--------------------------------------------------------------------------
    /** Compute the N-th integer root of a radicant.
     *  \tparam RADICANT  Radicant
     *  \tparam N         Order of root
     */
    template<unsigned RADICANT,unsigned N,int AUX=-1,bool less = false> 
    struct NthRoot 
    {
        // call with condition (AUX^N <= RADICANT)
        static const unsigned value=
            NthRoot<RADICANT,N,AUX-1,(MToTheN<AUX,N>::value <= RADICANT)>::value; 
    };
    
    //! \cond SKIPDOX
    template<unsigned RADICANT, unsigned N, int AUX> 
    struct NthRoot<RADICANT, N, AUX, true>
    {
        // here AUX^N <= RADICANT is true
        static const unsigned value = AUX+1; 
    };

    template<unsigned RADICANT, unsigned N> 
    struct NthRoot<RADICANT, N, -1, false> 
    {
        static const unsigned value = NthRoot<RADICANT,N,RADICANT,false>::value;
    }; 
    //! \endcond

    //----------------------------------------------------------------------
    /** \class base::IfElse
     *  Compile time conditional.
     *  The following code
     *  \code
     *  IfElse< trueOrFalse, A, B >::Type
     *  \endcode
     *  evaluates to the type \a A in case of \a true and
     *  to \a B in case of \a false.
     *  \tparam COND  Compile time condition to evaluate
     *  \tparam IF    Type for the true-case
     *  \tparam ELSE  Type for the false-case
     */
    template<bool COND, typename IF, typename ELSE> struct IfElse;

    //! \cond SKIPDOX
    template<typename IF, typename ELSE>
    struct IfElse<true, IF,ELSE> { typedef IF   Type; };
        
    template<typename IF, typename ELSE>
    struct IfElse<false,IF,ELSE> { typedef ELSE Type; };
    //! \endcond


}

#endif
