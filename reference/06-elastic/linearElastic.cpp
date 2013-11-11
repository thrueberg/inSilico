#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>

#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>
#include <base/fe/Basis.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/FieldBinder.hpp>

#include <base/solver/Eigen3.hpp>

#include <mat/hypel/StVenant.hpp>
#include <mat/Lame.hpp>
#include <solid/HyperElastic.hpp>
#include <solid/Stress.hpp>

#include <base/post/ErrorNorm.hpp>
#include <base/aux/FundamentalSolution.hpp>

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF, typename FUN>
void dirichletU( const typename base::Vector<DIM>::Type& x,
                 DOF* doFPtr, const FUN& fun )
{
    const typename base::Vector<DIM>::Type U = fun(x); 
    for ( unsigned d = 0; d < DIM; d++ )
        doFPtr -> constrainValue( d, U[d] );
}

//------------------------------------------------------------------------------
template<typename MESH, typename DISP, typename MATERIAL>
void writeVTKFile( const std::string& baseName,
                   const MESH&        mesh,
                   const DISP&        disp,
                   const MATERIAL&    material )
{
    // create file name with step number
    const std::string vtkFile = baseName + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );

    base::io::vtk::writePointData( vtkWriter, mesh, disp, "disp" );
    base::io::vtk::writeCellData( vtkWriter, mesh, disp, 
                                  boost::bind( solid::cauchy<typename MESH::Element,
                                                             typename DISP::Element,
                                                             MATERIAL>,
                                               _1, _2, material ), "sigma" );
            
    vtk.close();
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // material behaviour
    const double E  = 1000.0;
    const double nu = 0.25;

    
    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 1;
    const base::Shape shape    = base::QUAD;

    typedef mat::hypel::StVenant Material;

    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  mesh.smf \n";
        return 0;
    }

    // read name of input file
    const std::string meshFile = boost::lexical_cast<std::string>( argv[1] );

    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
    // define a mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;
    const unsigned dim = Mesh::Node::dim;

    // create a mesh and read from input
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // Use elastostatic fundmental solution as reference solution
    typedef base::aux::FundSolElastoStatic<dim> FSES;
    FSES fses( mat::Lame::lambda( E, nu ), mat::Lame::mu( E, nu ) );
    const FSES::VecDim y   = base::constantVector<dim>( -0.1 );
    FSES::VecDim dir = base::constantVector<dim>( 0.   );
    dir[0] = 0.; dir[1] = 1.;

    typedef boost::function< FSES::VecDoF( const FSES::VecDim& ) > SolFun;
    SolFun solFun = boost::bind( &FSES::fun, &fses, _1, y, dir );


    // quadrature objects for volume and surface
    const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    // Create a field
    const unsigned    doFSize = dim;
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize>           Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field field;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, field );

    // Creates a list of <Element,faceNo> pairs along the boundary
    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    // Create a boundary mesh from this list
    typedef base::mesh::BoundaryMeshBinder<Mesh::Element>::Type BoundaryMesh;
    BoundaryMesh boundaryMesh;
    {
        // Create a real mesh object from this list
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh, boundaryMesh );
    }
    
    // constrain the boundary
    base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                           meshBoundary.end(),
                                           mesh, field, 
                                           boost::bind( &dirichletU<dim,DoF,SolFun>,
                                                        _1, _2, solFun ) );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // material object
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );

    // matrix kernel
    typedef solid::HyperElastic<Material,FTB::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    // Number the degrees of freedom
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs << std::endl;


    // Create a solver object
    typedef base::solver::Eigen3           Solver;
    Solver solver( numDofs );

    // Compute element stiffness matrices and assemble them
    base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                 fieldBinder,
                                                 hyperElastic );

    // Finalise assembly
    solver.finishAssembly();


    // Solve
    solver.choleskySolve();
            
    // distribute results back to dofs
    base::dof::setDoFsFromSolver( solver, field );

    // write a vtk file
    writeVTKFile( baseName, mesh, field, material );
        
    // Finished load steps

    //--------------------------------------------------------------------------
    // compute L2-error and tell it to the user
    std::cout << "L2-error =  "
              << base::post::errorComputation<0>( quadrature, mesh, field, solFun )
              << "\n";
    
    return 0;
}
