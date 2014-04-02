//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Lagrange1D.ipp
//! @author Thomas Rueberg
//! @date   2012

namespace base{
    namespace sfun{

        //----------------------------------------------------------------------
        //! @name Constant function
        //@{

        //! \image html lagrange_0.png
        template<> inline
        void Lagrange1D<0>::fun( const VecDim& xi,
                                 FunArray& values ) const
        {
            values[0] = 1.;
        }

        // gradient and Hessian are not implemented on purpose
        // (trigger compile-time error)
        //@}

        //----------------------------------------------------------------------
        //! @name Linear functions
        //@{

        //! \image html lagrange_1.png        
        template<> inline
        void Lagrange1D<1>::fun( const VecDim& xi,
                                 FunArray& values ) const
        {
            values[0] = 1. - xi[0];
            values[1] =      xi[0];
        }

        template<>  inline
        void Lagrange1D<1>:: gradient( const VecDim& xi,
                                       GradArray& values ) const
        {
            values[0][0] = -1.;
            values[1][0] =  1.;
        }

        // Hessian not implement on purpose
        //@}

        //----------------------------------------------------------------------
        //! @name Quadratic functions
        //@{

        //! \image html lagrange_2.png
        template<> inline
        void Lagrange1D<2>::fun( const VecDim & xi,
                                 FunArray& values ) const
        {
            values[0] = (1. - xi[0]) * (1. - 2.*xi[0]);
            values[1] = 4. * xi[0] * (1. - xi[0]);
            values[2] = xi[0] * (2. * xi[0] - 1.);
        }

        template<> inline
        void Lagrange1D<2>::gradient( const VecDim& xi,
                                      GradArray& values ) const
        {
            values[0][0] = 4. * xi[0] - 3.;
            values[1][0] = 4. - 8. * xi[0];
            values[2][0] = 4. * xi[0] - 1.;
        }

        template<> inline
        void Lagrange1D<2>::hessian(  const VecDim& xi,
                                      HessianArray& values ) const
        {
            values[0](0,0) =  4.;
            values[1](0,0) = -8.;
            values[2](0,0) =  4.;
        }
        //@}


        //----------------------------------------------------------------------
        //! @name Cubic functions
        //@{

        //! \image html lagrange_3.png
        template<> inline
        void Lagrange1D<3>::fun( const VecDim& xi,
                                 FunArray& values ) const
        {
            const double zeta0 = 1. - xi[0];
            const double zeta1 =      xi[0];

            values[0] = 0.5 * zeta0 * (3. * zeta0 - 1.) * (3. * zeta0 - 2.);
            values[1] = 4.5 * zeta0 * (3. * zeta0 - 1.) * zeta1;
            values[2] = 4.5 * zeta1 * (3. * zeta1 - 1.) * zeta0;
            values[3] = 0.5 * zeta1 * (3. * zeta1 - 1.) * (3. * zeta1 - 2.);
        }

        template<> inline
        void Lagrange1D<3>::gradient( const VecDim& xi,
                                      GradArray& values ) const
        {
            const double zeta0 = 1. - xi(0);
            const double zeta1 =      xi(0);
            
            values[0][0] =
                -0.5 * ( (3.*zeta0 - 1.)*(3.*zeta0 - 2.) + 3.*zeta0 * (6.*zeta0 - 3.) );
            values[1][0] =
                -4.5 * ( zeta1 * (6.*zeta0 - 1.) - zeta0 * (3.*zeta0 - 1.) );
            values[2][0] =
                4.5 * ( zeta0 * (6.*zeta1 - 1.) - zeta1 * (3.*zeta1 - 1.) );
            values[3][0] =
                0.5 * ( (3.*zeta1 - 1.)*(3.*zeta1 - 2.) + 3.*zeta1 * (6.*zeta1 - 3.) );
        }

        template<> inline
        void Lagrange1D<3>::hessian(  const VecDim& xi,
                                      HessianArray& values ) const
        {
            values[0](0,0) = -9. * (3. * xi(0) - 2.);
            values[1](0,0) = -9. * (9. * xi(0) - 4.);
            values[2](0,0) =  9. * (3. * xi(0) - 1.);
            values[3](0,0) =  9. * (9. * xi(0) - 5.);
        }
        //@}

    }
}
