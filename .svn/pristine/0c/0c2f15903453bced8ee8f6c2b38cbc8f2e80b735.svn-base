
//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   funSpace.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_funspace_hpp
#define base_funspace_hpp

//------------------------------------------------------------------------------

namespace base{

    //--------------------------------------------------------------------------
    //! Classify the spaces of typical FE shape functions.
    enum FunSpace{
        LAGRANGE,       //!< Lagrange polynomials
        BSPLINE,        //!< Uniform B-Splines
        HERMITE,        //!< Hermite polynomials
        RAVIARTTHOMAS   //!< Raviart-Thomas vectorial shape functions
        // ... Nédélec, Crouzeix-Raviart, ...
    };

    //--------------------------------------------------------------------------
    //! Dimension of the functions return type.
    template<unsigned DIM, base::FunSpace SPACE>
    struct ValueDimension
    {
        static const unsigned value = 1;
    };

    //\cond SKIPDOX
    //! Raviart-Thomas functions have vector-valued return types
    template<unsigned DIM>
    struct ValueDimension<DIM,RAVIARTTHOMAS>
    {
        static const unsigned value = DIM;
    };

    // Nedelec etc.
    //\endcond

}

#endif
