#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>

#include <base/mesh/MeshBoundary.hpp>
#include <base/Quadrature.hpp>

#include <base/io/Format.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/fe/Basis.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>

#include <base/post/ErrorNorm.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/asmb/FieldIterator.hpp>
#include <base/asmb/StiffnessMatrix.hpp>

#include <base/asmb/FieldBinder.hpp>

#include <base/aux/FundamentalSolution.hpp>

#include <heat/Laplace.hpp>

//------------------------------------------------------------------------------
// Function for the point-wise constraint of the Boundary
//[dirichlet]{
template<typename FUN, typename DOF>
void dirichletBC( const typename FUN::arg1_type x, DOF* doFPtr, FUN fun )
{
    const typename FUN::result_type value = fun( x );
    
    if ( doFPtr -> isActive(0) ) {
        doFPtr -> constrainValue( 0, value[0] );
    }

}
//[dirichlet]}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << " file.smf \n\n";
        return -1;
    }
        
    const std::string smfFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string baseName = base::io::baseName( smfFile, ".smf" );

    //--------------------------------------------------------------------------
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 2;
    const base::Shape shape    = base::QUAD;
    
    const unsigned    doFSize = 1; // temperature

    //--------------------------------------------------------------------------
    typedef base::Unstructured<shape,geomDeg>       Mesh;
    const unsigned dim = Mesh::Node::dim;

    Mesh mesh;
    {
        std::ifstream smf( smfFile.c_str() );
        VERIFY_MSG( smf.is_open(), "Cannot open mesh file" );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Quadrature 
    const unsigned kernelDegEstimate = 3;
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

    // Use a fundamental solution
    const base::aux::FundSolLaplace<dim>::VecDim sourcePoint
        = base::constantVector<dim>( -.5 );

    typedef base::aux::FundSolLaplace<dim> FSol;
    FSol fSol;

    typedef boost::function< FSol::VecDoF( const FSol::VecDim& ) > FSolFun;
    FSolFun fSolFun = boost::bind( &FSol::fun, &fSol, _1, sourcePoint );
    
    // Object to constrain the boundary 
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, field, 
                                           boost::bind( &dirichletBC<FSolFun,
                                                                     Field::DegreeOfFreedom>,
                                                        _1, _2, fSolFun ) );

    // Number of DoFs after constraint application!
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << " Number of dofs " << numDofs << std::endl;
    
    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );

    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // Object providing the element stiffness matrix
    const double conductivity = 1.0;
    typedef heat::Laplace<FTB::Tuple> Laplace;
    Laplace laplace( conductivity );

    base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                 fieldBinder, laplace );

    // Finalise assembly
    solver.finishAssembly();

    // Solve
    solver.choleskySolve();

    // distribute results back to dofs
    base::dof::setDoFsFromSolver( solver, field );

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
    std::cout << "L2-error = "
              << base::post::errorComputation<0>( quadrature, mesh, field, fSolFun )
              << '\n';

    return 0;
}
