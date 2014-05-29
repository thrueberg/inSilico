#ifndef boundaryconditions_h
#define boundaryconditions_h

const double coordTol = 1.e-4;

//------------------------------------------------------------------------------
// Enum for switching the BCs at runtime
enum BC
{
    SHEAR,
    SEMISHEAR,
    CAVITY,
    PIPE,
    TANK
};

std::ostream& operator <<(std::ostream &os, const BC& bc )
{
    switch( bc ) {
        case SHEAR:      os << "SHEAR"    ; break;
        case SEMISHEAR:  os << "SEMISHEAR"; break;
        case CAVITY:     os << "CAVITY"   ; break;
        case PIPE:       os << "PIPE"     ; break;
        case TANK:       os << "TANK"     ; break;
        default:         os << "UNDEFINED"; 
    }
    return os; 
} 

std::istream& operator >>(std::istream &is, BC& bc )
{
    std::string buffer;
    is >> buffer;
    if      ( buffer == "SHEAR"     ) bc = SHEAR;
    else if ( buffer == "SEMISHEAR" ) bc = SEMISHEAR;
    else if ( buffer == "CAVITY"    ) bc = CAVITY;
    else if ( buffer == "PIPE"      ) bc = PIPE;
    else if ( buffer == "TANK"      ) bc = TANK;
    else VERIFY_MSG( false, "Unknown BC found" );
    return is; 
} 

//------------------------------------------------------------------------------
/**  Symmetric shear flow applied at top and bottom boundaries, rest open
 *
 *               o-> -> -> ->o- 
 *               ¦           ¦
 *               ¦   Imm.    ¦
 *         open  ¦   solid   ¦ open
 *               ¦   domain  ¦
 *               ¦           ¦
 *               o<- <- <- <-o- 
 */
template<unsigned DIM, typename DOF>
void shear( const typename base::Vector<DIM>::Type& x,
            DOF* doFPtr, const double factor,
            const BoundingBox<DIM>& bbox  ) 
{
    const bool top = bbox.isOnUpperBoundary( x, DIM-1, coordTol );
    const bool bot = bbox.isOnLowerBoundary( x, DIM-1, coordTol );

    if ( (not top) and (not bot) ) return;

    for ( unsigned d = 0; d < DIM; d++ ) {
        if ( doFPtr -> isActive(d) ) {
            const double value = ( d > 0 ? 0. : (top ? 1.0 : -1.0) ) * factor;
            doFPtr -> constrainValue( d, value );
        }
    }
}

//------------------------------------------------------------------------------
/** Apply linear shear inflow at left boundary and constant flow at top,
 *  bottom closed, right boundary open.
 *
 *               o--> --> -->  o- 
 *               |-->          ¦
 *    linear     |->   Imm.    ¦
 *    variation  |->   solid   ¦ open
 *               |>    domain  ¦
 *               |>            ¦
 *               o-------------o- 
 */
template<unsigned DIM, typename DOF>
void semiShear( const typename base::Vector<DIM>::Type& x,
                DOF* doFPtr, const double factor,
                const BoundingBox<DIM>& bbox ) 
{
    const bool left= bbox.isOnLowerBoundary( x, 0, coordTol );
    const bool top = bbox.isOnUpperBoundary( x, DIM-1, coordTol );
    const bool bot = bbox.isOnLowerBoundary( x, DIM-1, coordTol );

    if ( (not top) and (not bot) and (not left) ) return;

    for ( unsigned d = 0; d < DIM; d++ ) {
        if ( doFPtr -> isActive(d) ) {
            // normalised coordinate
            const typename BoundingBox<DIM>::VecDim xi =
                bbox.localCoordinate( x );
            const double value = ( d > 0 ? 0. : xi[DIM-1] ) * factor;
            doFPtr -> constrainValue( d, value );
        }
    }
}

//------------------------------------------------------------------------------
/**  Driven cavity flow on top boundary, rest closed walls
 *
 *               o-> -> -> ->o
 *               |           |
 *               |   Imm.    |
 *        closed |   solid   | closed
 *               |   domain  |
 *               |           |
 *               o-----------o
 *                   closed
 *
 */
template<unsigned DIM, typename DOF>
void cavity( const typename base::Vector<DIM>::Type& x,
             DOF* doFPtr, const double factor,
            const BoundingBox<DIM>& bbox  ) 

{
    // if d-th coordinate has the value 1.0
    bool activeBdry = bbox.isOnUpperBoundary( x, DIM-1, coordTol );
    // remove the corner/edge locations
    const typename BoundingBox<DIM>::VecDim xi = bbox.localCoordinate( x );
    for ( unsigned d = 0; d < DIM-1; d ++ ) {
        if ( std::abs( xi[d] - 0. ) < coordTol ) activeBdry = false;
        if ( std::abs( xi[d] - 1. ) < coordTol ) activeBdry = false;
    }

    // boundary condition is either 0 or the e_1 vector 
    for ( unsigned d = 0; d < DIM; d ++ ) {
        const double value = ( (d==0) and activeBdry ) ? factor : 0.0;
        if ( doFPtr -> isActive(d) ) 
            doFPtr -> constrainValue( d, value );
    }

}

//------------------------------------------------------------------------------
/** Poisseuille flow in a pipe 
 *
 *                        closed           
 *               o-------------------------o
 *               |>                        ¦
 *       parabol.|->     Imm.              ¦
 *               |-->    solid             ¦ open
 *        inflow |->     domain            ¦
 *               |>                        ¦
 *               o-------------------------o
 *                              closed      
 *
 */
template<unsigned DIM, typename DOF>
void pipe( const typename base::Vector<DIM>::Type& x,
           DOF* doFPtr, const double factor,
           const BoundingBox<DIM>& bbox  ) 
{
    VERIFY_MSG( DIM==2, "3D Poiseuille flow not implemented" );

    // Outflow on left-hand side
    const bool isNeumann = bbox.isOnUpperBoundary( x, 0, coordTol );
    if ( isNeumann ) return;
    
    const bool isInflow = bbox.isOnLowerBoundary( x, 0, coordTol );

    // boundary condition is either 0 or the e_1 vector
    const typename BoundingBox<DIM>::VecDim xi = bbox.localCoordinate( x );
    const double inflowU = 1.5 * factor * (1. - xi[1]*xi[1]);
    const double value = (isInflow ? inflowU : 0. );
    if ( doFPtr -> isActive(0) ) 
        doFPtr -> constrainValue( 0, value );
    
    for ( unsigned d = 1; d < DIM; d ++ ) {
        if ( doFPtr -> isActive(d) ) 
            doFPtr -> constrainValue( d, 0. );
    }
}

//------------------------------------------------------------------------------
/**  Tank 
 *
 *               .............  
 *               |           |
 *               |   Imm.    |
 *               |   solid   |     
 *               |   domain  | 
 *               |           |
 *               -------------  
 */
template<unsigned DIM, typename DOF>
void tank( const typename base::Vector<DIM>::Type& x,
           DOF* doFPtr, const double factor,
           const BoundingBox<DIM>& bbox  ) 
{
    const bool topBot =
        ( bbox.isOnUpperBoundary( x, DIM-1, coordTol ) or
          bbox.isOnLowerBoundary( x, DIM-1, coordTol ) );

    for ( unsigned d = 0; d < DIM; d++ ) {
        if ( topBot or (d==0) ) {
            if ( doFPtr -> isActive(d) ) {
                doFPtr -> constrainValue( d, 0. );
            }
        }
    }
}

template<unsigned DIM, typename DOF>
void tankP( const typename base::Vector<DIM>::Type& x,
            DOF* doFPtr, const BoundingBox<DIM>& bbox  ) 
{
    const bool top = bbox.isOnUpperBoundary( x, DIM-1, coordTol );

    if (top and (doFPtr -> isActive(0)) ) doFPtr -> constrainValue( 0, 0. );
}

//------------------------------------------------------------------------------
//! Apply BC in function of the given enum
template<unsigned DIM, typename DOF>
void dirichletBC( const typename base::Vector<DIM>::Type& x,
                  DOF* doFPtr, const double factor, 
                  const BoundingBox<DIM>& bbox, const BC& bc )
{
    if ( bc == SHEAR )
        return shear<DIM,DOF>( x, doFPtr, factor, bbox );
    else if ( bc == SEMISHEAR )
        return semiShear<DIM,DOF>( x, doFPtr, factor, bbox );
    else if ( bc == CAVITY )
        return cavity<DIM,DOF>( x, doFPtr, factor, bbox );
    else if ( bc == PIPE )
        return pipe<DIM,DOF>( x, doFPtr, factor,   bbox  );
    else if ( bc == TANK )
        return tank<DIM,DOF>( x, doFPtr, factor,   bbox  );
    else VERIFY_MSG( false, "Wrong BC identifier" );

    return;
}

//------------------------------------------------------------------------------
//! Set value of active degrees of freedom to zero.
template<unsigned DIM, typename DOF>
void fix( const typename base::Vector<DIM>::Type& x, DOF* doFPtr )
{
    if ( doFPtr -> isConstrained(0) ) {
        doFPtr -> clearConstraints();
    }

    for ( unsigned d = 0; d < DIM; d++ ) {
        if ( doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, 0. );
    }
}

//------------------------------------------------------------------------------
//! Set vertical value of active degrees of freedom to zero.
template<unsigned DIM, typename DOF>
void verticalFix( const typename base::Vector<DIM>::Type& x, DOF* doFPtr )
{
    if ( doFPtr -> isConstrained(DIM-1) ) doFPtr -> clearConstraints();

    if ( doFPtr -> isActive(DIM-1) ) doFPtr -> constrainValue( DIM-1, 0. );
}

//==============================================================================
template<unsigned DIM>
typename base::Vector<DIM>::Type
gravity( const typename base::Vector<DIM>::Type& x,
         const double value )
{
    typename base::Vector<DIM>::Type result = base::constantVector<DIM>( 0. );
    result[DIM-1] = -0.98 * value;
    return result;
}


//==============================================================================

#endif
