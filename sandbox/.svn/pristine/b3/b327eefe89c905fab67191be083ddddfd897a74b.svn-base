#include <base/Unstructured.hpp>
#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/generate.hpp>
#include <base/fe/Basis.hpp>

#include <base/shape.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/io/vtk/LegacyWriter.hpp>

//#include <mat/hypel/StVenant.hpp>
//#include <mat/hypel/CompNeoHookean.hpp>
#include "NeoHookeanCompressible.hpp"
#include <mat/Lame.hpp>

//#include <solid/Stress.hpp>

#include "HyperElasticMaterial.hpp"
#include "HyperElasticSpatial.hpp"

#include "InitialiseGeometry.hpp"


//------------------------------------------------------------------------------
template<typename MATERIAL, unsigned DIM, unsigned FDEG, unsigned GDEG=1>
class Solid
{
public:
    typedef MATERIAL Material;
    
    static const unsigned dim = DIM;
    static const unsigned fieldDeg = FDEG;
    static const unsigned geomDeg  = GDEG;

    static const unsigned nHist    = 1;

    static const base::Shape shape    = base::HyperCubeShape<dim>::value;

    // define a mesh
    typedef base::Unstructured<shape,geomDeg>    Mesh;


    // quadrature objects for volume and surface
    static const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;

    // define a field
    static const unsigned    doFSize = dim;
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>     Field;
    typedef typename Field::DegreeOfFreedom        DoF;
    
    // Bind the fields together
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    typedef typename FieldBinder::template TupleBinder<1,1>::Type FTB;

    // matrix kernel
    typedef solid::HyperElasticMaterial<Material,typename FTB::Tuple> HyperElasticMaterial;
    typedef solid::HyperElasticSpatial< Material,typename FTB::Tuple> HyperElasticSpatial;

    // create a mesh and read from input
    Solid( const std::string meshFile,
           const double E, const double nu )
        : E_( E ),
          material_( mat::Lame::lambda( E, nu), mat::Lame::mu( E, nu ) ),
          hyperElasticSpatial_(  material_ ),
          hyperElasticMaterial_( material_ ),
          fieldBinder_( mesh_, current_ )
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh_ );
        smf.close();

        // generate DoFs from mesh
        base::dof::generate<FEBasis>( mesh_, current_ );

        solid::initialiseGeometry( mesh_, current_ );
    }

    template<typename DIRIFUN>
    void constrainBoundary( const DIRIFUN& diriFun )
    {
        // Creates a list of <Element,faceNo> pairs along the boundary
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh_.elementsBegin(), mesh_.elementsEnd() );
        
        base::dof::constrainBoundary<FEBasis>( meshBoundary.begin(),
                                               meshBoundary.end(),
                                               mesh_, current_,
                                               diriFun );
    }

    std::size_t numberDoFs()
    {
        // Number the degrees of freedom
        numDoFs_ =
            base::dof::numberDoFsConsecutively( current_.doFsBegin(),
                                                current_.doFsEnd() );
        //std::cout << "# Number of dofs " << numDofs << std::endl;

        return numDoFs_;
    }


    //------------------------------------------------------------------------------
    void writeVTKFile( const std::string& baseName,
                       const unsigned     step ) const
    {
        // create file name with step number
        const std::string vtkFile =
            baseName + "." + base::io::leadingZeros( step ) + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh_ );
        
        base::io::vtk::writePointData( vtkWriter, mesh_, current_, "current" );

        //{
        //    std::vector<typename Mesh::Node::VecDim> oldDisp;
        //    base::post::evaluateHistoryAtNodes<1>( mesh_, current_, oldDisp );
        //    vtkWriter.writePointData( oldDisp.begin(), oldDisp.end(), "disp2" );
        //
        //}

#if 0
        typedef mat::hypel::SpatialWrapper<Material> AuxMat;
        AuxMat auxMat( material_ );
        
        base::io::vtk::writeCellData(  vtkWriter, mesh_, current_, 
                                       boost::bind( ::cauchy<typename Mesh::Element,
                                                    typename Field::Element,
                                                    AuxMat>, //Material>,
                                                    _1, _2,
                                                    auxMat ), "sigma" );
        //material_ ), "sigma" );
#endif       
        vtk.close();
    }


    template<typename BFUN>
    std::pair<double,double> iterate( const unsigned iter,
                                      const bool updated,
                                      const BFUN& bodyForceFun )
    {
        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDoFs_ );

        const bool zeroConstraints = iter > 0;

        if ( updated ) {
            
            base::asmb::computeResidualForces<FTB>( quadrature_, solver,
                                                    fieldBinder_,
                                                    hyperElasticSpatial_ );
            
            // Compute element stiffness matrices and assemble them
            base::asmb::stiffnessMatrixComputation<FTB>( quadrature_, solver,
                                                         fieldBinder_,
                                                         hyperElasticSpatial_,
                                                         zeroConstraints );

        }
        else{
            
            base::asmb::computeResidualForces<FTB>( quadrature_, solver,
                                                    fieldBinder_,
                                                    hyperElasticMaterial_ );
            
            // Compute element stiffness matrices and assemble them
            base::asmb::stiffnessMatrixComputation<FTB>( quadrature_, solver,
                                                         fieldBinder_,
                                                         hyperElasticMaterial_,
                                                         zeroConstraints );
        }
        
        base::asmb::bodyForceComputation<FTB>( quadrature_, solver,
                                               fieldBinder_, bodyForceFun );

        // Finalise assembly
        solver.finishAssembly();

        // norm of residual 
        const double conv1 = solver.norm();// / E_;


        {
            std::string solvName = "solver." + base::io::leadingZeros( iter );
            std::ofstream slv( solvName.c_str() );
            solver.debugLHS( slv );
            solver.debugRHS( slv );

        }

        // Solve
        solver.choleskySolve();
            
        // distribute results back to dofs
        //if ( iter > 0 ) 
        base::dof::addToDoFsFromSolver( solver, current_, zeroConstraints );
            //else 
            
        // update?

        // norm of displacement increment
        const double conv2 = solver.norm();

        return std::make_pair( conv1, conv2 );
    }

    void updateGeometry()
    {
#if 0
        typedef typename Mesh::Node::VecDim VecDim;
        
        // get nodal displacements
        std::vector<VecDim> nodalDisp;
        base::post::evaluateAtNodes( mesh_, current_, nodalDisp );

        // add displacements to reference configuration
        typename Mesh::NodePtrIter nIter = mesh_.nodesBegin();
        typename Mesh::NodePtrIter nEnd  = mesh_.nodesEnd();
        typename std::vector<VecDim>::iterator dispIter = nodalDisp.begin();
        for ( ; nIter != nEnd; ++nIter, ++dispIter )
        {
            const VecDim u = *dispIter;

            VecDim X;
            (*nIter) -> getX( &(X[0]) );

            const VecDim x = X + u;
            (*nIter) -> setX( &(x[0]) );
        }

        //base::dof::clearDoFs( current_ );

        typename Field::DoFPtrIter dIter = current_.doFsBegin();
        typename Field::DoFPtrIter dEnd  = current_.doFsEnd();
        for ( ; dIter != dEnd; ++dIter ) {

            for ( unsigned d = 0; d < dim; d++ ) {
                const double dU = (*dIter) -> getValue( d );
                const double  U = (*dIter) ->template getHistoryValue<1>( d );
                (*dIter) ->template setHistoryValue<1>( d, U + dU );
                (*dIter) -> setValue( d, 0. );
            }

        }
#endif   
        return;
    }

    //--------------------------------------------------------------------------
    Mesh&  accessMesh() { return mesh_; }
    Field& accessDisplacement() { return current_; }
    Quadrature& accessQuadrature() { return quadrature_; }
    
private:
    const double E_;
    Mesh         mesh_;
    Quadrature   quadrature_;
    Field        current_;
    Material     material_;
    
    HyperElasticSpatial  hyperElasticSpatial_;
    HyperElasticMaterial hyperElasticMaterial_;
    FieldBinder  fieldBinder_;


    std::size_t numDoFs_;
};
