#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/post/ErrorNorm.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <heat/Static.hpp>
#include <mat/thermal/FenicsTest.hpp>
#include <mat/thermal/IsotropicConstant.hpp>

//------------------------------------------------------------------------------
namespace ref04 {

    //--------------------------------------------------------------------------
    //! Reference solution, see mat::thermal::FenicsTest
    template<unsigned DIM>
    base::Vector<1,double>::Type
    referenceSolution( const typename base::Vector<DIM,double>::Type& x,
                       const double m )
    {
        const double fac = std::pow( 2.0, m+1.0 );
        const double mantissa = (fac-1.0) * x[0] + 1.0;
        const double result = std::pow( mantissa, 1./(m+1.0) ) - 1.0;
        return base::Vector<1,double>::Type::Constant( result );
    }

    //--------------------------------------------------------------------------
    //! Set simple boundary conditions u(0) = 0 and u(1) = 1
    template<unsigned DIM, typename DOF>
    void simpleBC( const typename base::Vector<DIM,double>::Type& x,
                   DOF* doFPtr ) 
    {
        const double tol = 1.e-5;
        const bool left  = ( std::abs( x[0] - 0. ) < tol );
        const bool right = ( std::abs( x[0] - 1. ) < tol );

        if ( left and ( doFPtr -> isActive(0) ) )
            doFPtr -> constrainValue( 0, 0. );

        if ( right and ( doFPtr -> isActive(0) ) )
            doFPtr -> constrainValue( 0, 1. );

    }

    int nonlinearHeat( int argc, char * argv[] );
}

//------------------------------------------------------------------------------
/** Solve a simple Sturm-Liouville equation.
 *  Here an essentially one-dimensional boundary value problem with a mild
 *  non-linearity is solved. The BVP is described in mat::thermal::FenicsTest
 *  and tackled here with a Newton iteration.
 *
 *  The main feature of this function is
 *  - nonlinear iteration with a Newton method
 *  - check of convergence
 *
 *  Execution
 *  \code{.txt}
 *  ./nonlinearHeat square_200.smf
 *  \endcode
 *  gives the output of
 *  \include nonlinearHeat.ref.dat
 *
 *  \param[in] argc Number of command line arguments
 *  \param[in] argv Values of command line arguments
 */
int ref04::nonlinearHeat( int argc, char * argv[] )
{
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf \n\n";
        return -1;
    }
        
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    const unsigned maxIter   = 10;
    const double   tolerance = 1.e-8;
    const double   solM      = 3.0;

    //--------------------------------------------------------------------------
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;
    
    const unsigned    doFSize = 1; // temperature

    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned    dim     = Mesh::Node::dim;

    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature 
    const unsigned kernelDegEstimate = 5;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // Finite element basis
    typedef base::fe::Basis<shape,fieldDeg> FEBasis;

    // DOF handling
    typedef base::Field<FEBasis,doFSize>           Field;
    Field field;

    base::dof::generate<FEBasis>( mesh, field );

    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, field, 
                                           boost::bind( &simpleBC<dim,
                                                        Field::DegreeOfFreedom>,
                                                        _1, _2 ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << " Number of dofs " << numDofs << std::endl;

    //typedef mat::thermal::IsotropicConstant Material;
    typedef mat::thermal::FenicsTest Material;
    Material material( solM );

    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type FieldTupleBinder;

    
    typedef heat::Static<Material,FieldTupleBinder::Tuple> StaticHeat;
    StaticHeat staticHeat( material );

    unsigned iter = 0;
    while ( iter < maxIter ) {

        std::cout << iter << ":  " << std::flush;
    
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDofs );

        // Compute system matrix
        base::asmb::stiffnessMatrixComputation<FieldTupleBinder>( quadrature, solver,
                                                                  fieldBinder,
                                                                  staticHeat );

        // compute residual forces
        base::asmb::computeResidualForces<FieldTupleBinder>( quadrature, solver,
                                                             fieldBinder,
                                                             staticHeat);

        const double residualNorm = solver.norm();
        std::cout << "|R| = " << residualNorm << std::flush;
        if ( residualNorm < tolerance ) {
            std::cout << std::endl;
            break;
        }
        
        // Finalise assembly
        solver.finishAssembly();

        // Solve a possibly non-symmetric system
        solver.superLUSolve();
        //solver.luSolve();
        //solver.biCGStabSolve();

        const double incrementNorm = solver.norm();
        if ( incrementNorm < tolerance ) break;
        std::cout << ",  |x| = " << incrementNorm << std::endl;

        // distribute results back to dofs
        base::dof::addToDoFsFromSolver( solver, field );

        iter++;
    }

    // message to user
    if ( iter < maxIter ) {
        std::cout << "Converged to given tolerance in " << iter
                  << " iterations \n";
    }
    else {
        std::cout << "Reached maximal allowed number of iterations\n ";
    }

    // output to a VTK file
    {
        // VTK Legacy
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh );
        base::io::vtk::writePointData( vtkWriter, mesh, field, "temperature" );
        vtk.close();
    }

    // compute L2-error and tell it to the user
    std::cout << "\nL2-error = "
              << base::post::errorComputation<0>(
                  quadrature, mesh, field,
                  boost::bind( &referenceSolution<dim>, _1, solM ) )
              << '\n';
    
    
    return 0;
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    return ref04::nonlinearHeat( argc, argv );
}
