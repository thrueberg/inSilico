#ifndef chemical_hpp
#define chemical_hpp

#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/Quadrature.hpp>
#include <base/io/Format.hpp>

#include <base/fe/Basis.hpp>
#include <base/cut/ScaledField.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/setField.hpp>

#include <base/dof/Distribute.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/BodyForce.hpp>

#include <base/io/vtk/LegacyWriter.hpp>

#include <heat/Laplace.hpp>
#include <heat/Convection.hpp>

#include <base/solver/Eigen3.hpp>

#include <base/cut/LevelSet.hpp>
#include <base/cut/bruteForce.hpp>
#include <base/cut/Cell.hpp>
#include <base/cut/Quadrature.hpp>
#include <base/cut/generateCutCells.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/cut/ComputeSupport.hpp>
#include <base/cut/TransferSurfaceDatum.hpp>
#include <base/cut/ComputeSurfaceForces.hpp>

#include <base/nitsche/Penalty.hpp>
#include <base/nitsche/Energy.hpp>
#include <base/nitsche/Parameters.hpp>


//------------------------------------------------------------------------------
template<typename MESH, unsigned DEG=1>
class Chemical
{
public:
    //! @name Template parameter
    //@{
    typedef MESH Mesh;
    static const unsigned fieldDeg = DEG;
    //@}

    //! Shape
    static const base::Shape shape = Mesh::Element::shape;

    typedef typename Mesh::Node::VecDim VecDim;

    //! Type of quadrature
    static const unsigned kernelDegEstimate = 3;
    typedef base::cut::Quadrature<kernelDegEstimate,shape>   Quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;

    
    //! @name Chemical Field
    //@{
    static const unsigned    doFSize = 1;
    typedef base::fe::Basis<shape,fieldDeg>          FEBasis;
    typedef base::cut::ScaledField<FEBasis,doFSize>  Concentration;
    //@}

    // Solver type
    typedef base::solver::Eigen3           Solver;

    //
    typedef typename base::cut::LevelSet<dim> LevelSet;
    typedef typename base::cut::Cell<shape>   Cell;

    //! Constructor
    FluidBox( const Mesh& mesh, 
              const double density, const double diffusion,
              const bool inside = true )
        : mesh_( mesh ),
          isInside_( inside ), 
          quadrature_( cells_, isInside_ ),
          density_( density ), diffusion_( diffusion )
    { 
        // generate DoFs from mesh
        base::dof::generate<FEBasis>( mesh_, concentration_ );
        cells_.resize( std::distance( mesh_.elementsBegin(), mesh_.elementsEnd() ) );
    }

    //! Constrain the box-boundary values of the velocity
    template<typename FUN>
    void constrainConcentration( const FUN& fun )
    {
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh_.elementsBegin(), mesh_.elementsEnd() );

        base::dof::constrainBoundary<FEBasis>( meshBoundary.boundaryBegin(),
                                               meshBoundary.boundaryEnd(),
                                               mesh_, concentration_, fun );
    }

    //! Immerse a surface
    template<typename SMESH>
    void immerseSurface( const SMESH& surfaceMesh, const double suppTolerance = 1.e-8 )
    {
        concentration_.activateAll();
        
        // compute signed-distance function
        const bool isSigned = true;
        levelSet_.resize( 0 );
        base::cut::bruteForce( mesh_, surfaceMesh, isSigned, levelSet_ );

        // generate structure of cut-cells
        base::cut::generateCutCells( mesh_, levelSet_, cells_ );
        
        {
            typedef base::asmb::FieldBinder<Mesh,Concentration>   TmpField;
            typedef typename TmpField::template TupleBinder<1,1>::Type TFTB;

            TmpField field_( mesh_, concentration_ );
            
            const std::size_t numDoFs =
                std::distance( concentration_.doFsBegin(), concentration_.doFsEnd() );
            supportSizes_.resize( numDoFs );
    
            base::cut::supportComputation<TFTB>( field_, quadrature_,
                                                 supportSizes_, isInside_ );
            velocity_.scaleAndTagBasis( supportSizes_, suppTolerance );
        }
        
    }

    
    //! Generate dofs and assign dof numbers
    std::size_t numberDoFs( )
    {
        numDoFs_ =
            base::dof::numberDoFsConsecutively( concentration_.doFsBegin(),
                                                concentration_.doFsEnd() );
        
        return numDoFs_;
    }

    template<typename SMESH, typename VELOCITY>
    void solveProblem( const SMESH&    surfaceMesh,
                       VELOCITY& velocity_,
                       const double   deltaT,
                       const unsigned numStep, 
                       const double   penaltyFactor, 
                       const unsigned maxIter, 
                       const double   tolerance,
                       const bool     dirichlet )
    {
        base::dof::clearDoFs( concentration_ );
        
        // generate mesh from immersed boundary
        typedef typename base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
        SurfaceMesh immersedMesh;
        base::cut::generateSurfaceMesh<Mesh,Cell>( mesh_, cells_, immersedMesh );

        // bind concentration to velocity
        typedef base::asmb::FieldBinder<Mesh,Concentration,VELOCITY>     Field;
        typedef typename Field::template TupleBinder<1,1,2>::Type FTB;
        Field field( mesh_, concentration_, velocity_ );

        // bind a field around it
        typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Concentration> SurfaceFieldBinder;
        SurfaceFieldBinder immersedFieldBinder( immersedMesh, concentration_)

        typedef typename SurfaceFieldBinder::template TupleBinder<1,1>::Type STBU;
        
        //typename base::cut::TransferSurfaceDatum<SMESH,SVELOC,
        //                                         typename Mesh::Element>
        //    s2d( surfaceMesh, surfaceVelocity, levelSet_ );

        
        // Create a solver object
        Solver solver( numDoFsU_ + numDoFsP_ );
    
        // Compute system matrix
        heat::Laplace<typename FTB::Tuple> laplace( diffusion_ );
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature_, solver,
                                                     field, laplace );
        

        
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature_, solver,
                                                     field, convection );

        // compute inertia terms, d/dt, due to time integration
        base::time::computeInertiaTerms<FTB,MSM>( quadrature_, solver,
                                                  field, deltaT, numStep,
                                                  density_ );

        heat::Convection<typename FTB::Tuple> convection( density );
        base::asmb::stiffnessMatrixComputation<FTB>( quadrature_, solver,
                                                     field, convection );

        
        base::nitsche::ImmersedBoundary<Cell> ib( 1.0, cells_ );

        if ( dirichlet ) {
            
            base::nitsche::penaltyLHS<STBUU>( surfaceQuadrature_, solver,
                                              immersedFieldBinder, ib, penaltyFactor );
            
            base::nitsche::energyLHS<STBUU>( laplace, surfaceQuadrature_, solver,
                                             immersedFieldBinder, ib, isInside_ );
            //base::nitsche::penaltyRHS2<STBUU>( surfaceQuadrature_, solver,
            //                                   immersedFieldBinder, s2d, ib, penaltyFactor );
            //base::nitsche::energyRHS2<STBUU>( stressDivergence_, surfaceQuadrature_, solver,
            //                                  immersedFieldBinder, s2d, ib, isInside_ );

        }

        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.superLUSolve();

        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, concentration_ );

    }


    //--------------------------------------------------------------------------
    void writeVTKFile( const std::string& baseName,
                       const unsigned step )
    {
        // VTK Legacy
        const std::string vtkFile = baseName
            + "." + base::io::leadingZeros( step ) + ".vtk";
        
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh_ );

        base::io::vtk::writePointData( vtkWriter, mesh_, concentration_, "C" );
        
        {
            std::vector<double> distances;
            std::transform( levelSet_.begin(), levelSet_.end(),
                            std::back_inserter( distances ),
                            boost::bind( &LevelSet::getSignedDistance, _1 ) );
            vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
        }

        vtk.close();
    }


    //--------------------------------------------------------------------------
    Mesh&          accessMesh()          { return mesh_; }
    Concentration& accessConcentration() { return velocity_; }

    std::vector<double>& supportSizes() { return supportSizes_; }
    
private:
    const double isInside_;
    
    Mesh              mesh_;
    Quadrature        quadrature_;
    SurfaceQuadrature surfaceQuadrature_;
    Concentration     concentration_;

    std::size_t numDoFs_;

    std::vector<LevelSet> levelSet_;
    std::vector<Cell>     cells_;

    std::vector<double>   supportSizes_;

    const double density_;
    const double diffusion_;
};


#endif
