//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   LagrangeTetrahedron.ipp
//! @author Thomas Rueberg
//! @date   2012

namespace base{
    namespace sfun{

        //----------------------------------------------------------------------
        //! @name Constant function
        //@{
        template<> inline
        void LagrangeTetrahedron<0>::fun( const VecDim& xi,
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
        template<> inline
        void LagrangeTetrahedron<1>::fun( const VecDim& xi,
                                               FunArray& values ) const
        {
            values[0] = 1. - xi[0] - xi[1] - xi[2];
            values[1] = xi[0];
            values[2] = xi[1];
            values[3] = xi[2];
        }

        template<> inline
        void LagrangeTetrahedron<1>::gradient( const VecDim& xi,
                                               GradArray& values ) const
        {
            values[0][0] = -1.; values[0][1] = -1.; values[0][2] = -1.;
            values[1][0] =  1.; values[1][1] =  0.; values[1][2] =  0.;
            values[2][0] =  0.; values[2][1] =  1.; values[2][2] =  0.;
            values[3][0] =  0.; values[3][1] =  0.; values[3][2] =  1.;
        }
        
        template<> inline
        void LagrangeTetrahedron<1>::supportPoints( boost::array<VecDim,
                                                                 numFun> & values )
        {
            values[0][0] = 0.; values[0][1] = 0.; values[0][2] = 0.;
            values[1][0] = 1.; values[1][1] = 0.; values[1][2] = 0.;
            values[2][0] = 0.; values[2][1] = 1.; values[2][2] = 0.;
            values[3][0] = 0.; values[3][1] = 0.; values[3][2] = 1.;
        }
        //@}

        //----------------------------------------------------------------------
        //! @name Quadratic functions
        //@{
        template<> inline
        void LagrangeTetrahedron<2>::fun( const VecDim& xi,
                                               FunArray& values ) const
        {
            // Baricentric coordinates
            const double z1 = xi[0];
            const double z2 = xi[1];
            const double z3 = xi[2];
            const double z0 = 1. - z1 - z2 - z3;

            values[0] = z0 * (2. * z0 - 1.);
            values[1] = z1 * (2. * z1 - 1.);
            values[2] = z2 * (2. * z2 - 1.);
            values[3] = z3 * (2. * z3 - 1.);
            values[4] = 4. * z0 * z1;
            values[5] = 4. * z1 * z2;
            values[6] = 4. * z2 * z0;
            values[7] = 4. * z0 * z3;
            values[8] = 4. * z1 * z3;
            values[9] = 4. * z2 * z3;
        }

        template<> inline
        void LagrangeTetrahedron<2>:: gradient( const VecDim& xi,
                                                GradArray& values ) const
        {
            // Barycentric coordinates
            const double z1 = xi[0];
            const double z2 = xi[1];
            const double z3 = xi[2];
            const double z0 = 1. - z1 - z2 - z3;

            // Derivatives of Barycentric with respect to Cartesian coordinates
            const double dz0dxi = -1.; const double dz0deta = -1.; const double dz0dzeta = -1.;
            const double dz1dxi =  1.; const double dz1deta =  0.; const double dz1dzeta =  0.;
            const double dz2dxi =  0.; const double dz2deta =  1.; const double dz2dzeta =  0.;
            const double dz3dxi =  0.; const double dz3deta =  0.; const double dz3dzeta =  1.;

            // d()/dxi
            values[0][0] = (4. * z0 - 1.) * dz0dxi;
            values[1][0] = (4. * z1 - 1.) * dz1dxi;
            values[2][0] = (4. * z2 - 1.) * dz2dxi;
            values[3][0] = (4. * z3 - 1.) * dz3dxi;
            values[4][0] = 4. * (z0 * dz1dxi + dz0dxi * z1);
            values[5][0] = 4. * (z1 * dz2dxi + dz1dxi * z2);
            values[6][0] = 4. * (z0 * dz2dxi + dz0dxi * z2);
            values[7][0] = 4. * (z0 * dz3dxi + dz0dxi * z3);
            values[8][0] = 4. * (z1 * dz3dxi + dz1dxi * z3);
            values[9][0] = 4. * (z2 * dz3dxi + dz2dxi * z3);

            // d()/deta
            values[0][1] = (4. * z0 - 1.) * dz0deta;
            values[1][1] = (4. * z1 - 1.) * dz1deta;
            values[2][1] = (4. * z2 - 1.) * dz2deta;
            values[3][1] = (4. * z3 - 1.) * dz3deta;
            values[4][1] = 4. * (z0 * dz1deta + dz0deta * z1);
            values[5][1] = 4. * (z1 * dz2deta + dz1deta * z2);
            values[6][1] = 4. * (z0 * dz2deta + dz0deta * z2);
            values[7][1] = 4. * (z0 * dz3deta + dz0deta * z3);
            values[8][1] = 4. * (z1 * dz3deta + dz1deta * z3);
            values[9][1] = 4. * (z2 * dz3deta + dz2deta * z3);

            // d()/dzeta
            values[0][2] = (4. * z0 - 1.) * dz0dzeta;
            values[1][2] = (4. * z1 - 1.) * dz1dzeta;
            values[2][2] = (4. * z2 - 1.) * dz2dzeta;
            values[3][2] = (4. * z3 - 1.) * dz3dzeta;
            values[4][2] = 4. * (z0 * dz1dzeta + dz0dzeta * z1);
            values[5][2] = 4. * (z1 * dz2dzeta + dz1dzeta * z2);
            values[6][2] = 4. * (z0 * dz2dzeta + dz0dzeta * z2);
            values[7][2] = 4. * (z0 * dz3dzeta + dz0dzeta * z3);
            values[8][2] = 4. * (z1 * dz3dzeta + dz1dzeta * z3);
            values[9][2] = 4. * (z2 * dz3dzeta + dz2dzeta * z3);

        }

        template<> inline
        void LagrangeTetrahedron<2>::supportPoints( boost::array<VecDim,
                                                                 numFun> & values )
        {
            values[0][0] = 0.; values[0][1] = 0.; values[0][2] = 0.;
            values[1][0] = 1.; values[1][1] = 0.; values[1][2] = 0.;
            values[2][0] = 0.; values[2][1] = 1.; values[2][2] = 0.;
            values[3][0] = 0.; values[3][1] = 0.; values[3][2] = 1.;

            values[4][0] = 0.5; values[4][1] = 0.0; values[4][2] = 0.0;
            values[5][0] = 0.5; values[5][1] = 0.5; values[5][2] = 0.0;
            values[6][0] = 0.0; values[6][1] = 0.5; values[6][2] = 0.0;

            values[7][0] = 0.0; values[7][1] = 0.0; values[7][2] = 0.5;
            values[8][0] = 0.5; values[8][1] = 0.0; values[8][2] = 0.5;
            values[9][0] = 0.0; values[9][1] = 0.5; values[9][2] = 0.5;
        }
        //@}
        
    }
}
