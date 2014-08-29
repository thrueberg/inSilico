#include <iostream> // input- and output streams
#include <string>   // string objects
#include <fstream>  // file streams
#include <boost/lexical_cast.hpp> // lexical cast between objects

//------------------------------------------------------------------------------

#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/Quadrature.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/post/ErrorNorm.hpp>

#include <base/nitsche/Parameters.hpp>
#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>

#include <base/kernel/Laplace.hpp>

#include "Helper.hpp"

STATIC_ASSERT_MSG( (SPACEDIM>1) and (SPACEDIM<4), "Inapt choice of dimension" );

//------------------------------------------------------------------------------
// namespace for the analytic solution of the problem
namespace analytic{

    //--------------------------------------------------------------------------
    // simple 1D functions
    const double alpha = 5.;
    double fun(   const double x ) { return                std::sin( alpha * x ); }
    double funD(  const double x ) { return        alpha * std::cos( alpha * x ); }
    double funDD( const double x ) { return -alpha*alpha * std::sin( alpha * x ); }

    //--------------------------------------------------------------------------
    // Tensor-product realisations of 1D solution function
    template<unsigned DIM>
    base::Vector<1>::Type 
    solution( const typename base::Vector<DIM>::Type& x )
    {
        base::Vector<1>::Type u;
        u[0] = 1.;
        
        for ( unsigned d = 0; d < DIM; d++ ) u[0] *= fun( x[d] );

        return u;
    }

    // gradient of the tensor product function
    template<unsigned DIM>
    typename base::Matrix<DIM,1>::Type
    gradient( const typename base::Vector<DIM>::Type& x )
    {
        typename base::Matrix<DIM,1>::Type gradU;
        
        for ( unsigned d = 0; d < DIM; d++ ) {
            gradU[d] = funD( x[d] );
            
            for ( unsigned d2 = 0; d2 < DIM; d2++ ) {
                if ( d2 != d ) gradU[d] *= fun( x[d2] );
            }
        }
        
        return gradU;
    }

    // boundary fluxes
    template<unsigned DIM>
    base::Vector<1>::Type 
    neumannBC( const typename base::Vector<DIM>::Type& x,
               const typename base::Vector<DIM>::Type& normal )
    {
        const typename base::Vector<DIM>::Type gradU = gradient<DIM>( x );
        return (gradU.transpose() * normal );
    }


    // negative Laplacian gives the force term
    template<unsigned DIM>
    base::Vector<1>::Type
    forceFun( const typename base::Vector<DIM>::Type& x ) 
    {
        base::Vector<1>::Type f;
        f[0] = 0.;

        for ( unsigned d = 0; d < DIM; d++ ) {
            double aux = funDD( x[d] );
            for ( unsigned d2 = 0; d2 < DIM; d2++ ) {
                if ( d2 != d ) aux *= fun( x[d2] );
            }
            f[0] += aux;
        }

        return -f;
    }

} // end namespace analytic


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    const unsigned dim      = SPACEDIM;
    const unsigned geomDeg  = 1;
    const unsigned fieldDeg = 1;

#ifdef STRUCTURED
    typedef apps::nitsche::StructuredHelper<dim,geomDeg,fieldDeg>   Helper;
#else
    typedef apps::nitsche::UnstructuredHelper<dim,geomDeg,fieldDeg> Helper;
#endif
    
    // Check the number of input arguments
    if ( argc != 2 ) { 
        std::cout << "Usage:  " << argv[0] << " file" << Helper::suffix() << "\n";
        return 0;
    }
    
    // convert input argument to a string object
    const std::string fileName = boost::lexical_cast<std::string>( argv[1] );
    // extract the basename of the file
    const std::string baseName =
        fileName.substr( 0, fileName.find( Helper::suffix() ) );

    //--------------------------------------------------------------------------
    const double penaltyFactor = 50.;
    const bool   nitsche = true; 
    const unsigned doFSize  = 1;

    //--------------------------------------------------------------------------
    typedef Helper::Mesh Mesh;
    Mesh mesh;
    {
        std::ifstream smgf( fileName.c_str() );
        Helper::readFromFile( smgf, mesh );
    }

    //--------------------------------------------------------------------------
    // Quadrature 
    const unsigned kernelDegEstimate = 5;
    typedef base::Quadrature<kernelDegEstimate,Helper::shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,Helper::shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // Finite element basis
    typedef Helper::FEBasis FEBasis;
    
    // DOF handling
    typedef base::Field<FEBasis,doFSize>           Field;
    Field field;
    base::dof::generate<FEBasis>( mesh, field );

#ifdef NEUMANN
    // fix a dof
    (*field.doFsBegin()) -> constrainValue( 0, 0. );
#endif

    // Number of DoFs 
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    //std::cout << " Number of dofs " << numDofs << std::endl;

    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type FieldTupleBinder;

    // Body force
    base::asmb::bodyForceComputation<FieldTupleBinder>(
        quadrature, solver, fieldBinder,
        boost::bind( &analytic::forceFun<dim>, _1 ) );

    //--------------------------------------------------------------------------
    // Creates a list of <Element,faceNo> pairs
    base::mesh::MeshBoundary meshBoundary;
    Helper::meshBoundary( meshBoundary, mesh );

    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
    }

    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, field );
    typedef SurfaceFieldBinder::TupleBinder<1,1>::Type STB;

    //--------------------------------------------------------------------------
    // Stiffness matrix
    typedef base::kernel::Laplace<FieldTupleBinder::Tuple> Laplace;
    Laplace laplace( 1. );
    base::asmb::stiffnessMatrixComputation<FieldTupleBinder>( quadrature, solver,
                                                              fieldBinder, laplace );



#ifdef NEUMANN
    // Neumann boundary condition
    base::asmb::neumannForceComputation<STB>( surfaceQuadrature, solver,
                                              surfaceFieldBinder,
                                              boost::bind( &analytic::neumannBC<dim>, _1, _2 ) );
#else
    // Dirichlet BCs a la Nitsche
    base::nitsche::OuterBoundary ob( 1. );
    base::nitsche::penaltyLHS<STB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                    ob, penaltyFactor );
    
    base::nitsche::penaltyRHS<STB>( surfaceQuadrature, solver, surfaceFieldBinder,
                                    boost::bind( &analytic::solution<dim>, _1 ), ob,
                                    penaltyFactor );

    if ( nitsche ) {
    
        base::nitsche::energyLHS<STB>( laplace, surfaceQuadrature, solver,
                                       surfaceFieldBinder, ob );
        base::nitsche::energyRHS<STB>( laplace, surfaceQuadrature, solver,
                                       surfaceFieldBinder,
                                       boost::bind( &analytic::solution<dim>, _1 ), ob );
    }
#endif

    // Finalise assembly
    solver.finishAssembly();

    // Solve
    //solver.choleskySolve();
    const unsigned cgIter = solver.cgSolve();

    // distribute results back to dofs
    base::dof::setDoFsFromSolver( solver, field );

    //--------------------------------------------------------------------------
    // VTK file of the input mesh
    {
        const std::string vtkMeshFileName = baseName + ".vtk";
        std::ofstream vtk( vtkMeshFileName.c_str() );

        Helper::writeVTK( vtk, mesh, field );
        vtk.close();
    }
    
    // compute L2-error and tell it to the user
    std::cout 
        << base::post::errorComputation<0>( quadrature, mesh, field,
                                            boost::bind( &analytic::solution<dim>, _1 ) )
        << "  " 
        << base::post::errorComputation<1>( quadrature, mesh, field,
                                            boost::bind( &analytic::gradient<dim>, _1 ) )
        << "  "
        << cgIter
        << '\n';

    
    return 0;
}
