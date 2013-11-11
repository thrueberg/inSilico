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
#include <base/io/PropertiesParser.hpp>
#include <base/io/vtk/LegacyWriter.hpp>
#include <base/io/Format.hpp>

#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>
#include <base/fe/Basis.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <solid/HyperElastic.hpp>

#include <surf/HyperElasticMembrane.hpp>
#include <surf/Skalak.hpp>
#include <surf/NeoHookean.hpp>


//------------------------------------------------------------------------------
template<typename ELEMENT>
typename base::Vector<ELEMENT::Node::dim>::Type
internalPressure( const ELEMENT* gep, 
                  const typename ELEMENT::GeomFun::VecDim& xi,
                  const double value )
{
    typename base::Vector<ELEMENT::Node::dim>::Type normal;
    base::SurfaceNormal<ELEMENT>()( gep, xi, normal );
    return value * normal;
}

//------------------------------------------------------------------------------
// Statically determined positioning
template<typename FIELD>
void constrainSupports( FIELD& field )
{
    VERIFY_MSG( (FIELD::DegreeOfFreedom::size==2),
                "This method only works for 2D" );
    
    typename FIELD::DoFPtrIter dIter = field.doFsBegin();

    const std::size_t numDoFs = std::distance( dIter, field.doFsEnd() );

    VERIFY_MSG( (numDoFs % 4) == 0, "Trying to fix the quarter points" );

    (*dIter) -> constrainValue( 1, 0. );

    std::advance( dIter, numDoFs/4 );
    (*dIter) -> constrainValue( 0, 0. );
        
    std::advance( dIter, numDoFs/4 );
    (*dIter) -> constrainValue( 1, 0. );

    std::advance( dIter, numDoFs/4 );
    (*dIter) -> constrainValue( 0, 0. );

}

//------------------------------------------------------------------------------
// Assume the circle centre to be at (0,0)
template<typename MESH, typename FIELD>
double giveRadius( const MESH& mesh, const FIELD& field )
{
    typename MESH::NodePtrConstIter nIter = mesh.nodesBegin();
    
    typename FIELD::DoFPtrConstIter dIter = field.doFsBegin();

    return
        ((*nIter) -> getX())[0] +
        ((*dIter) -> getValue(0));
}

//------------------------------------------------------------------------------
template<typename MESH, typename DISP, typename MATERIAL>
void writeVTKFile( const std::string& baseName,
                   const unsigned     step,
                   const MESH&        mesh,
                   const DISP&        disp,
                   const MATERIAL&    material )
{
    // create file name with step number
    const std::string vtkFile =
        baseName + "." + base::io::leadingZeros( step ) + ".vtk";
    std::ofstream vtk( vtkFile.c_str() );
    base::io::vtk::LegacyWriter vtkWriter( vtk );
    vtkWriter.writeUnstructuredGrid( mesh );

    base::io::vtk::writePointData( vtkWriter, mesh, disp, "disp" );
    vtk.close();
}


//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    // basic attributes of the computation
    const unsigned    geomDeg  = 1;
    const unsigned    fieldDeg = 1;
    const unsigned    dim      = 2;
    const base::Shape shape    = base::SimplexShape<dim-1>::value;

    // usage message
    if ( argc != 2 ) {
        std::cout << "Usage:  " << argv[0] << "  input.dat \n";
        return 0;
    }

    // read name of input file
    const std::string inputFile = boost::lexical_cast<std::string>( argv[1] );

    // read from input file
    std::string meshFile;
    double A, B, tolerance, pressure;
    unsigned maxIter, loadSteps;
    {    
        //Feed properties parser with the variables to be read
        base::io::PropertiesParser prop;
        prop.registerPropertiesVar( "meshFile",         meshFile );
        prop.registerPropertiesVar( "A",                A );
        prop.registerPropertiesVar( "B",                B );
        prop.registerPropertiesVar( "maxIter",          maxIter );
        prop.registerPropertiesVar( "loadSteps",        loadSteps );
        prop.registerPropertiesVar( "tolerance",        tolerance );
        prop.registerPropertiesVar( "pressure",         pressure );

        // Read variables from the input file
        std::ifstream inp( inputFile.c_str()  );
        VERIFY_MSG( inp.is_open(), "Cannot open input file" );
        prop.readValues( inp );
        inp.close( );

        // Make sure all variables have been found
        if ( not prop.isEverythingRead() ) {
            prop.writeUnread( std::cerr );
            VERIFY_MSG( false, "Could not find above variables" );
        }
    }

    // find base name from mesh file
    const std::string baseName = base::io::baseName( meshFile, ".smf" );

    //--------------------------------------------------------------------------
    // define a mesh
    typedef base::Unstructured<shape,geomDeg,dim>    Mesh;

    // create a mesh and read from input
    Mesh mesh;
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh );
        smf.close();
    }

    // quadrature objects for volume and surface
    const unsigned kernelDegEstimate = 1; //3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;
    
    // Create a field
    const unsigned    doFSize = dim;
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize>           Field;
    typedef Field::DegreeOfFreedom                 DoF;
    Field field;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, field );


    // constrain a dof
    constrainSupports( field );

    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    FieldBinder fieldBinder( mesh, field );
    typedef FieldBinder::TupleBinder<1,1>::Type FTB;

    // material object
    typedef surf::Skalak Energy;
    Energy   energy( A, B );

    //typedef surf::NeoHookean Energy;
    //Energy   energy( A  );

    typedef surf::HyperElasticMembrane<Energy> Material;
    Material material( energy );

    // matrix kernel
    typedef solid::HyperElastic<Material,FTB::Tuple> HyperElastic;
    HyperElastic hyperElastic( material );

    // Number the degrees of freedom
    const std::size_t numDofs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );
    std::cout << "# Number of dofs " << numDofs << std::endl;

    // write a vtk file
    writeVTKFile( baseName, 0, mesh, field, material );

    // load increment
    const double deltaP = pressure / static_cast<double>( loadSteps );

    std::cout << "# pressure   radius \n"
              << 0.0  << "   "  << giveRadius( mesh, field ) << "\n";
    
    //--------------------------------------------------------------------------
    // Loop over load steps
    //--------------------------------------------------------------------------
    for ( unsigned step = 0; step < loadSteps; step++ ) {

        // current load
        const double p = (step+1) * deltaP;

        //----------------------------------------------------------------------
        // Nonlinear iterations
        //----------------------------------------------------------------------
        unsigned iter = 0;
        while ( iter < maxIter ) {

            //std::cout << "Iter " << iter;

            // Create a solver object
            typedef base::solver::Eigen3           Solver;
            Solver solver( numDofs );

            
            base::asmb::computeResidualForces<FTB>( quadrature, solver,
                                                    fieldBinder, hyperElastic );
            
            // Compute element stiffness matrices and assemble them
            base::asmb::stiffnessMatrixComputation<FTB>( quadrature, solver,
                                                         fieldBinder, hyperElastic,
                                                         iter > 0 );

            // internal pressure load
            base::asmb::bodyForceComputation2<FTB>( quadrature, solver, fieldBinder,
                                                    boost::bind( &internalPressure<Mesh::Element>,
                                                                 _1, _2, p ) );


            // Finalise assembly
            solver.finishAssembly();

            // norm of residual 
            const double conv1 = solver.norm();

            //std::cout << " |F| = " << conv1;

            // convergence via residual norm
            if ( conv1 < tolerance * A ) { // note the tolerance multiplier
                //std::cout << "\n";
                break;
            }

            // Solve
            solver.choleskySolve();
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, field, iter > 0 );

            // norm of displacement increment
            const double conv2 = solver.norm();
            //std::cout << " | DU | = " << conv2 << "\n";
            
            // convergence via increment
            if ( conv2 < tolerance ) break;

            iter++;

        }
        // Finished non-linear iterations
        //----------------------------------------------------------------------

        // warning
        if ( iter == maxIter ) {
            std::cout << "# (WW) Step " << step << " has not converged within "
                      << maxIter << " iterations \n";
        }

        // write a vtk file
        writeVTKFile( baseName, step+1, mesh, field, material );

        std::cout << p << " " << giveRadius( mesh, field ) << std::endl;

    }
    // Finished load steps
    //--------------------------------------------------------------------------

    
    return 0;
}
