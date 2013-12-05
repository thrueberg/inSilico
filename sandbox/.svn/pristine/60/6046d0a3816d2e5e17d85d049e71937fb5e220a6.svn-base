//------------------------------------------------------------------------------
// Tolerance for coordinate comparison
static const double coordTol = 1.e-6;

//------------------------------------------------------------------------------
template<unsigned DIM>
class BoundaryConditions
{
public:
    // spatial dimension
    static const unsigned dim = DIM;

    // current implementation wants this
    STATIC_ASSERT_MSG( dim == 3, "Implementation requires dim=3" );

    // convenience typedefs
    typedef typename base::Vector<DIM>::Type VecDim;
    typedef typename base::Vector<1>::Type   VecDoF;

    // boundary points with x_1 > 0.7 are on the Neumann boundary
    // or top and bottom points;
    static bool isNeumann( const VecDim & x ) 
    {
        const bool minX = ( std::abs( x[0] + 0.5 ) < coordTol );
        const bool maxX =  x[0] > 0.3;
        const bool maxZ =
            ( std::abs( x[2] - 0.    ) < coordTol ) or
            ( std::abs( x[2] + 0.145 ) < coordTol );
	return maxX or ( maxZ and not minX );
    }

    // boundary points with x_1 < 0.8 are on the Dirichlet boundary
    static bool isDirichlet( const VecDim& x )
    {
        return not( isNeumann( x ) );
    }

    // boundary points with x_1 = -0.5 have a non-zero value
    static bool hasAppliedValue( const VecDim& x )
    {
	
	/* All sides but bottom
	const bool boundx1 = std::abs( x[0] + 0.5 ) < coordTol;
	const bool boundx2 = std::abs( x[0] - 0.5 ) < coordTol;
        const bool boundy1 = std::abs( x[1] - 0.5 ) < coordTol;
	const bool boundz1 = std::abs( x[2] + 0.5 ) < coordTol;
	const bool boundz2 = std::abs( x[2] - 0.5 ) < coordTol;
	return boundx1 or boundx2 or boundy1 or boundz1 or boundz2; 
    */
	const bool boundx1 = std::abs( x[0] - 0.5 ) < coordTol;
	const bool boundy1 = std::abs( x[1] + 0.5 ) < coordTol;
	const bool boundz1 = std::abs( x[2] - 0.5 ) < coordTol;
	return (boundx1 or boundy1);
    }

    // Dirichlet boundary condition
    template<typename DOF>
    static void dirichleBC( const VecDim& x, DOF* doFPtr )
    {
        const double appliedValue = 20.0;

        if ( isDirichlet( x ) ) {

            const double value =
                (hasAppliedValue(x) ? appliedValue : 0. );
            
            if ( doFPtr -> isActive(0) )
                doFPtr -> constrainValue( 0, value );
        }
        return;
    }

    // Body force function
    static VecDoF forceFun( const VecDim& x )
    {
        VecDoF result = base::constantVector<1>( 0. );
        return result;
    }

    // Flux boundary condition
    static VecDoF neumannBC( const VecDim& x, const VecDim& normal )
    {
        VecDoF result = base::constantVector<1>( 1. );
        return result;
    }

    // Initial condition
    template<typename DOF>
    static void initialState( const VecDim& x, DOF* doFPtr ) 
    {
        const double value = 0.;
        doFPtr -> setValue( 0, value );
    }
    
};
	    
//-----------------------------------------------------------------------
// Diffusion Coefficients - IMCOMPLETE
template<typename ELEMENT>
double diffusionConstant( const ELEMENT* geomEp,
		          const typename ELEMENT::GeomFun::VecDim& xi,
                          const double D1, const double D2 )
{
    const typename ELEMENT::Node::VecDim x =
		base::Geometry<ELEMENT>()( geomEp, xi );

    const bool domain2 = (x[0] >= (-0.1) and x[0] <= (0.1));

    return (domain2 ? D2 : D1);
}


//--------------------------------------------------------------------------
// output to a VTK file
template<typename MESH, typename FIELD>
void writeVTKFile( const MESH& mesh, const FIELD& field, 
                   const std::string outFileName )
{
    std::ofstream vtk(  outFileName.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );
    base::io::vtk::writePointData( vtkWriter, mesh, field, "temperature" );
    vtk.close();
}
    

