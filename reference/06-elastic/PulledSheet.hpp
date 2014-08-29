//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   PulledSheet
//! @author Thomas Rueberg
//! @date   2014

#ifndef elasticity_pulledsheet_hpp
#define elasticity_pulledsheet_hpp

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace ref06{
    template<unsigned DIM> class PulledSheet;
}

//------------------------------------------------------------------------------
/** Test problem for compressible hyperelasticity.
 *  Block of material, occupying (0,1)^DIM, fix x_1 = 0 and pull at x_1 = 1.
 *  Optionally, at x_1 = 1, a surface traction is applied or a normal
 *  displacement.
 */
template<unsigned DIM>
class ref06::PulledSheet
{
public:
    typedef typename base::Vector<DIM>::Type VecDim;

    // Fix x_0=0 and optionally pull at x_1=1
    template<typename DOF>
    static void dirichletBC( const VecDim& x, DOF* doFPtr,
                             const bool pullRightSide,
                             const double value ) 
    {
        // location at x_1 = 0 or x_1 = 1
        const bool onLeftBdr =  ( std::abs( x[0] -  0. ) < coordTol );
        const bool onRightBdr = ( std::abs( x[0] -  1. ) < coordTol );

        // Fix left boundary at x_0 = 0
        if ( onLeftBdr ) {
            for ( unsigned d = 0; d < DOF::size; d++ ) {
                if ( doFPtr -> isActive(d) )
                    doFPtr -> constrainValue( d, 0.0 );
            }
        }

        // If assked for, apply normal displacement at x_1=1
        if (  onRightBdr and pullRightSide ) {
            if ( doFPtr -> isActive(0) )
                doFPtr -> constrainValue( 0, value );
        }

        return;
    }

    // apply surface traction at x_1 = 1 
    static VecDim neumannBC( const VecDim& x,
                             const VecDim& normal,
                             const double value )
    {
        VecDim result = VecDim::Constant( 0. );

        const double coordTol = 1.e-5;
        
        const bool onTractionBdr =
            ( std::abs( x[0] -  1. ) < coordTol );

        if ( onTractionBdr )
            // result[1] = value;
            result[0] = value * (x[1] - 0.5);
        
        return result;
    }
    
};


#endif
