//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   LagrangeTriangle.ipp
//! @author Thomas Rueberg
//! @date   2012

namespace base{
    namespace sfun{

        //----------------------------------------------------------------------
        //! @name Constant function
        //@{
        template<>
        void LagrangeTriangle<0>::fun( const VecDim& xi,
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
        template<>
        void LagrangeTriangle<1>::fun( const VecDim& xi,
                                            FunArray& values ) const
        {
            values[0] = 1. - xi[0] - xi[1];
            values[1] = xi[0];
            values[2] = xi[1];
        }

        template<>
        void LagrangeTriangle<1>::gradient( const VecDim& xi,
                                            GradArray& values ) const
        {
            values[0][0] = -1.; values[0][1] = -1.;
            
            values[1][0] =  1.; values[1][1] =  0.;
            
            values[2][0] =  0.; values[2][1] =  1.; 
        }

        template<>
        void LagrangeTriangle<1>::supportPoints( boost::array<VecDim,
                                                              numFun> & values )
        {
            values[0][0] = 0.; values[0][1] = 0.;
            
            values[1][0] = 1.; values[1][1] = 0.;
            
            values[2][0] = 0.; values[2][1] = 1.;
        }
        //@}

        //----------------------------------------------------------------------
        //! @name Quadratic functions
        //@{
        template<>
        void LagrangeTriangle<2>::fun( const VecDim& xi,
                                            FunArray& values ) const
        {
            values[ 0 ] = (1. - xi[0] - xi[1]) * (1. - 2.*xi[0] - 2.*xi[1] );
            values[ 1 ] = xi[0] * ( 2. * xi[0] - 1.);
            values[ 2 ] = xi[1] * ( 2. * xi[1] - 1.);

            values[ 3 ] = 4. * xi[0] * (1. - xi[0] - xi[1] );
            values[ 4 ] = 4. * xi[0] * xi[1];
            values[ 5 ] = 4. * xi[1] * (1. - xi[0] - xi[1] );
        }

        template<>
        void LagrangeTriangle<2>:: gradient( const VecDim& xi,
                                             GradArray& values ) const
        {
            values[0][0] = 4. * xi[0] + 4. * xi[1] - 3.;
            values[0][1] = values[0][0];

            values[1][0] = 4. * xi[0] - 1.;
            values[1][1] = 0.;

            values[2][0] = 0.;
            values[2][1] = 4. * xi[1] - 1.;

            values[3][0] = 4. - 8. * xi[0] - 4. * xi[1];
            values[3][1] = - 4. * xi[0];

            values[4][0] = 4. * xi[1];
            values[4][1] = 4. * xi[0];

            values[5][0] = - 4. * xi[1];
            values[5][1] = 4. - 4. * xi[0] - 8. * xi[1];
        }

        template<>
        void LagrangeTriangle<2>::hessian( const VecDim& xi,
                                           HessianArray& values ) const
        {
            values[0](0, 0) = 4.; values[0](0, 1) = values[0](1, 0) = 4.; values[0](1, 1) = 4.;
            values[1](0, 0) = 4.; values[1](0, 1) = values[1](1, 0) = 0.; values[1](1, 1) = 0.;
            values[2](0, 0) = 0.; values[2](0, 1) = values[2](1, 0) = 0.; values[2](1, 1) = 4.;
            values[3](0, 0) =-8.; values[3](0, 1) = values[3](1, 0) =-4.; values[3](1, 1) = 0.;
            values[4](0, 0) = 0.; values[4](0, 1) = values[4](1, 0) = 4.; values[4](1, 1) = 0.;
            values[5](0, 0) = 0.; values[5](0, 1) = values[5](1, 0) =-4.; values[5](1, 1) =-8.;
        }

        template<>
        void LagrangeTriangle<2>::supportPoints( boost::array<VecDim,
                                                              numFun> & values )
        {
            values[0] = base::constantVector<dim>( 0. );
            values[1][0] = 1.;  values[1][1] = 0.;
            values[2][0] = 0.;  values[2][1] = 1.;

            values[1][0] = 0.5; values[1][1] = 0.0;
            values[1][0] = 0.5; values[1][1] = 0.5;
            values[1][0] = 0.0; values[1][1] = 0.5;
        }
        //@}
        
    }
}
