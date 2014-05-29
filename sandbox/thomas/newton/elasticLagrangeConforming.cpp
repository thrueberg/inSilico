#include <string>
#include <boost/lexical_cast.hpp>

#include <base/verify.hpp>
#include <base/Unstructured.hpp>
#include <base/Field.hpp>

#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/mesh/MeshBoundary.hpp>

#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>
#include <base/io/smf/Reader.hpp>

#include <base/Quadrature.hpp>
#include <base/fe/Basis.hpp>

#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/location.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <solid/HyperElastic.hpp>
#include <mat/hypel/NeoHookeanCompressible.hpp>
#include <mat/Lame.hpp>

#include <base/post/findLocation.hpp>
#include <base/post/Monitor.hpp>

const double coordTol = 1.e-6;

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
void writeVTKFile( const std::string& baseName,
                   const unsigned step,
                   const MESH&    mesh,
                   const FIELD&   displacement )
{
    const std::string vtkFile = baseName + "." + 
        base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );

    vtkWriter.writeUnstructuredGrid( mesh );
    base::io::vtk::writePointData( vtkWriter, mesh, displacement, "disp" );
    vtk.close();
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename base::Vector<DIM>::Type
tractionFun( const typename base::Vector<DIM>::Type& x,
             const typename base::Vector<DIM>::Type& normal,
             const double value )
{
    const double phi = std::atan2( x[1], x[0] );
    
    
    typename base::Vector<DIM>::Type f = value * normal
        * (std::cos(2. * phi));
    return f;
}


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    // spatial dimension
    const unsigned    dim = 2; 
    
    if ( argc != 3 ) {
        std::cerr << "Usage: " << argv[0] << " file.smf  input.dat\n"
                  << "(Compiled for dim=" << dim << ")\n\n";
        return -1;
    }

    // read name of input file
    const std::string   meshFile  = boost::lexical_cast<std::string>( argv[1] );
    const std::string   inputFile = boost::lexical_cast<std::string>( argv[2] );

    // read from input file
    double E, nu, f, tolerance;
    unsigned numLoadSteps, maxIter;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "E",            E );
        prop.registerPropertiesVar( "nu",           nu );
        prop.registerPropertiesVar( "f",            f );
        prop.registerPropertiesVar( "numLoadSteps", numLoadSteps );
        prop.registerPropertiesVar( "maxIter",      maxIter );
        prop.registerPropertiesVar( "tolerance",    tolerance );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        VERIFY_MSG( prop.readValuesAndCheck( inp ), "Input error" );
        inp.close( );
    }

    // basic attributes of the computation
    const unsigned             geomDeg  = 1;
    const unsigned             fieldDeg = 1;
    const base::Shape             shape = base::HyperCubeShape<dim>::value;
    const unsigned    kernelDegEstimate = 3;
    const unsigned              doFSize = dim;
    const unsigned              nHist   = 1;

    typedef base::Unstructured<shape,geomDeg>  Mesh;
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
    }
    typedef Mesh::Node::VecDim VecDim;

    base::mesh::MeshBoundary meshBoundary;
    meshBoundary.create( mesh.elementsBegin(), mesh.elementsEnd() );

    typedef base::mesh::BoundaryMeshBinder<Mesh>::Type SurfaceMesh;
    SurfaceMesh boundaryMesh;
    {
        base::mesh::generateBoundaryMesh( meshBoundary.begin(),
                                          meshBoundary.end(),
                                          mesh,
                                          boundaryMesh );
                                         
    }

    //--------------------------------------------------------------------------
    // FE
    typedef base::fe::Basis<shape,fieldDeg>               FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>            Field;
    Field displacement;
    base::dof::generate<FEBasis>( mesh, displacement );

    // for domain field
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, displacement );
    typedef FieldBinder::TupleBinder<1,1>::Type UU;

    // for surface field
    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Field> SurfaceFieldBinder;
    SurfaceFieldBinder surfaceFieldBinder( boundaryMesh, displacement );
    typedef SurfaceFieldBinder::TupleBinder<1>::Type SU;

    //--------------------------------------------------------------------------
    // Quadratures
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;
    SurfaceQuadrature surfaceQuadrature;

    writeVTKFile( "conform", 0, mesh, displacement );

    // fix some points
    {
        std::vector<std::pair<std::size_t,VecDim> > doFLocation;
        base::dof::associateLocation( displacement, doFLocation );
        
        VecDim point = base::constantVector<dim>( 0. );
        const std::size_t index = base::dof::findDoFWithLocation( doFLocation, mesh,
                                                                  point, coordTol );
        
        Field::DegreeOfFreedom* dPtr = displacement.doFPtr( index );
        for ( unsigned d = 0; d < dim; d++ ) dPtr -> constrainValue( d, 0. );

        if ( dim > 1 ) {
            point[0] = 1.0;
            const std::size_t index = base::dof::findDoFWithLocation( doFLocation, mesh,
                                                                      point, coordTol );
            
            Field::DegreeOfFreedom* dPtr = displacement.doFPtr( index );
            for ( unsigned d = 1; d < dim; d++ ) dPtr -> constrainValue( d, 0. );
        }

    }
    
    // number DoFs
    const std::size_t activeDoFsU = 
        base::dof::numberDoFsConsecutively( displacement.doFsBegin(), displacement.doFsEnd() );

    typedef mat::hypel::NeoHookeanCompressible Material;
    Material material( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) );
    typedef solid::HyperElastic<Material,UU::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < numLoadSteps; step++ ) {

        const double fFun = f *
            static_cast<double>(step+1) / static_cast<double>( numLoadSteps );

        for ( unsigned iter = 0; iter < maxIter; iter++ ) {

            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( activeDoFsU );

            base::asmb::computeResidualForces<UU>( quadrature, solver,
                                                   fieldBinder, hyperElastic );

            base::asmb::stiffnessMatrixComputation<UU>( quadrature, solver,
                                                        fieldBinder, hyperElastic,
                                                        iter > 0 );

            base::asmb::neumannForceComputation<SU>( surfaceQuadrature,
                                                     solver, surfaceFieldBinder,
                                                     boost::bind( &tractionFun<dim>,
                                                                  _1, _2, fFun ) );


            // Finalise assembly
            solver.finishAssembly();

            // norm of residual
            const double conv1 = solver.norm() / E;
            std::cout << "* " << iter << " " << conv1 << " ";

            if ( conv1 < tolerance ) {
                std::cout << std::endl;
                break;
            }

            // Solve
            solver.superLUSolve();
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, displacement, iter > 0 );

            const double conv2 = solver.norm();
            std::cout << conv2 << std::endl;

            if ( conv2 < tolerance ) break;
        }

        // write a vtk file
        writeVTKFile( "conform", step+1, mesh, displacement );


        std::cout << step << "  " << fFun << "  ";
        std::cout << std::endl;

    }

    return 0;
}
