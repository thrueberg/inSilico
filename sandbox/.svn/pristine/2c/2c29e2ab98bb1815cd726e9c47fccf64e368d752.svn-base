#ifndef fluidbox_hpp
#define fluidbox_hpp

#include <fstream>
#include <string>
#include <utility>

#include <base/shape.hpp>
#include <base/Unstructured.hpp>

#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/Quadrature.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/smf/Writer.hpp>
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

#include <base/io/Format.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <fluid/Stokes.hpp>
#include <fluid/Convection.hpp>

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

#include "Robin.hpp"

//------------------------------------------------------------------------------
template<unsigned DIM, typename DOF>
void fixDoF( const typename base::Vector<DIM,double>::Type& x,
             DOF* doFPtr,
             const typename base::Vector<DIM,double>::Type& where,
             const double tol,
             const std::vector<double>& supportSizes,
             const double value, 
             bool& jobDone )
{
    if ( jobDone ) return; // one is enough

    const double supportSize = supportSizes[ doFPtr -> getID() ];
    const double aux = value * supportSize;
    
    if ( (x-where).norm() < tol ) {

        jobDone = true;
        
        for ( unsigned d = 0; d < DOF::size; d++ )
            doFPtr -> constrainValue( d, aux );
    }
}

//------------------------------------------------------------------------------

template<unsigned DIM, unsigned GDEG=1, unsigned UDEG=2, unsigned PDEG=1>
class FluidBox
{
public:
    //! @name Template parameter
    //@{
    static const unsigned dim       = DIM;
    static const unsigned geomDeg   = GDEG;
    static const unsigned fieldDegU = UDEG;
    static const unsigned fieldDegP = PDEG;
    //@}

    //! Shape
    static const base::Shape shape = base::SimplexShape<dim>::value;

    //! Type of mesh
    typedef base::Unstructured<shape,geomDeg,dim>    Mesh;

    typedef typename Mesh::Node::VecDim VecDim;

    //! Type of quadrature
    static const unsigned kernelDegEstimate = 3;
    typedef base::cut::Quadrature<kernelDegEstimate,shape>   Quadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape> SurfaceQuadrature;

    
    //! @name Displacement Field
    //@{
    static const unsigned    doFSizeU = dim;
    static const unsigned    doFSizeP = 1;
    typedef base::fe::Basis<shape,fieldDegU>          FEBasisU;
    typedef base::fe::Basis<shape,fieldDegP>          FEBasisP;
    typedef base::cut::ScaledField<FEBasisU,doFSizeU> Velocity;
    typedef base::cut::ScaledField<FEBasisP,doFSizeP> Pressure;
    //@}

    //! @name Field binding
    //@{
    typedef base::asmb::FieldBinder<Mesh,Velocity,Pressure>   Field;
    typedef typename Field::template TupleBinder<1,1,1>::Type TopLeft;
    typedef typename Field::template TupleBinder<1,2>::Type   TopRight;
    typedef typename Field::template TupleBinder<2,1>::Type   BottomLeft;
    typedef typename Field::template TupleBinder<2,2>::Type   BottomRight;
    //@}

    // kernels
    typedef fluid::StressDivergence<  typename    TopLeft::Tuple> StressDivergence;
    typedef fluid::Convection<        typename    TopLeft::Tuple> Convection;
    typedef fluid::PressureGradient<  typename   TopRight::Tuple> GradP;
    typedef fluid::VelocityDivergence<typename BottomLeft::Tuple> DivU;

    // Solver type
    typedef base::solver::Eigen3           Solver;

    //
    typedef typename base::cut::LevelSet<dim> LevelSet;
    typedef typename base::cut::Cell<shape>   Cell;

    //! Constructor
    FluidBox( const double density, const double viscosity,
              const bool inside = true )
        : isInside_( inside ), 
          quadrature_( cells_, isInside_ ),
          field_( mesh_, velocity_, pressure_ ),
          stressDivergence_( viscosity ),
          convection_(       density ),
          gradP_(), divU_()
    { }

    //! Read mesh from file
    void readMesh( const std::string& meshFile )
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh_ );
        smf.close();
        
        // generate DoFs from mesh
        base::dof::generate<FEBasisU>( mesh_, velocity_ );
        base::dof::generate<FEBasisP>( mesh_, pressure_ );

        cells_.resize( std::distance( mesh_.elementsBegin(), mesh_.elementsEnd() ) );

    }

    //! Constrain the box-boundary values of the velocity
    template<typename FUN>
    void constrainVelocity( const FUN& fun )
    {
        base::mesh::MeshBoundary meshBoundary;
        meshBoundary.create( mesh_.elementsBegin(), mesh_.elementsEnd() );

        base::dof::constrainBoundary<FEBasisU>( meshBoundary.begin(),
                                                meshBoundary.end(),
                                                mesh_, velocity_, fun );
    }

    //! Fix the first pressure dof for solvability
    void fixPressureDoF( const VecDim& where, const double tol, const double value )
    {
        // 
        bool jobDone = false;
        
        base::dof::setField( mesh_, pressure_,
                             boost::bind( &fixDoF<dim,
                                          typename Pressure::DegreeOfFreedom>, _1, _2,
                                          where, tol,
                                          boost::ref( supportSizesP_ ),
                                          value, 
                                          boost::ref(jobDone) ) );

        VERIFY_MSG( jobDone, "Could not find a pressure node" );

    }

    //! Immerse a surface
    template<typename SMESH>
    void immerseSurface( const SMESH& surfaceMesh, const double suppTolerance = 1.e-8 )
    {
        velocity_.activateAll();
        pressure_.activateAll();
        
        // compute signed-distance function
        const bool isSigned = true;
        levelSet_.resize( 0 );
        base::cut::bruteForce( mesh_, surfaceMesh, isSigned, levelSet_ );

        // generate structure of cut-cells
        base::cut::generateCutCells( mesh_, levelSet_, cells_ );

        {
            const std::size_t numDoFs =
                std::distance( velocity_.doFsBegin(), velocity_.doFsEnd() );
            supportSizesU_.resize( numDoFs );
    
            base::cut::supportComputation<TopLeft>( field_, quadrature_,
                                                     supportSizesU_, isInside_ );
            velocity_.scaleAndTagBasis( supportSizesU_, suppTolerance );
        }

        {
            const std::size_t numDoFs =
                std::distance( pressure_.doFsBegin(), pressure_.doFsEnd() );
            supportSizesP_.resize( numDoFs );
    
            base::cut::supportComputation<BottomRight>( field_, quadrature_,
                                                        supportSizesP_, isInside_ );
            pressure_.scaleAndTagBasis( supportSizesP_, suppTolerance );
        }
    }

    
    //! Generate dofs and assign dof numbers
    std::pair<std::size_t,std::size_t> numberDoFs( )
    {
        numDoFsU_ =
            base::dof::numberDoFsConsecutively( velocity_.doFsBegin(),
                                                velocity_.doFsEnd() );
        
        numDoFsP_ =
            base::dof::numberDoFsConsecutively( pressure_.doFsBegin(),
                                                pressure_.doFsEnd(),
                                                numDoFsU_ );

        return std::make_pair( numDoFsU_, numDoFsP_ );
    }

    //template<typename SMESH, typename SVELOC>
    unsigned solveProblem( //const SMESH&   surfaceMesh,
                           //const SVELOC&  surfaceVelocity,
                           const double   penaltyFactor, 
                           const unsigned maxIter, 
                           const double   tolerance,
                           const bool     withConvection )
    {
        base::dof::clearDoFs( velocity_ );
        base::dof::clearDoFs( pressure_ );
        
        // generate mesh from immersed boundary
        typedef typename base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
        SurfaceMesh immersedMesh;
        base::cut::generateSurfaceMesh<Mesh,Cell>( mesh_, cells_, immersedMesh );

        // bind a field around it
        typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Velocity,Pressure> SurfaceFieldBinder;
        SurfaceFieldBinder immersedFieldBinder( immersedMesh, velocity_, pressure_ );

        typedef typename SurfaceFieldBinder::template TupleBinder<1,1,1>::Type STBUU;
        typedef typename SurfaceFieldBinder::template TupleBinder<1,2>::Type   STBUP;


        
        //typename base::cut::TransferSurfaceDatum<SMESH,SVELOC,
        //                                         typename Mesh::Element>
        //    s2d( surfaceMesh, surfaceVelocity, levelSet_ );

        
        // Nonlinear Picard iterations
        unsigned iter = 0;
        double prevNorm = 0.;
        bool   converged = false;
        while( (iter < maxIter) and not converged ) {

            // Create a solver object
            Solver solver( numDoFsU_ + numDoFsP_ );
    
            // Compute system matrix
            base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature_, solver,
                                                             field_, stressDivergence_ );

            if ( withConvection )
                base::asmb::stiffnessMatrixComputation<TopLeft>( quadrature_, solver,
                                                                 field_, convection_ );
            else converged = true;

            base::asmb::stiffnessMatrixComputation<TopRight>( quadrature_, solver,
                                                              field_, gradP_ );

            base::asmb::stiffnessMatrixComputation<BottomLeft>( quadrature_, solver,
                                                                field_, divU_ );

            base::nitsche::ImmersedBoundary<Cell> ib( 1.0, cells_ );
            
            base::nitsche::penaltyLHS<STBUU>( surfaceQuadrature_, solver,
                                              immersedFieldBinder, ib, penaltyFactor );
            //base::nitsche::penaltyRHS2<STBUU>( surfaceQuadrature_, solver,
            //                                   immersedFieldBinder, s2d, ib, penaltyFactor );
            
            base::nitsche::energyLHS<STBUU>( stressDivergence_, surfaceQuadrature_, solver,
                                             immersedFieldBinder, ib, isInside_ );
            //base::nitsche::energyRHS2<STBUU>( stressDivergence_, surfaceQuadrature_, solver,
            //                                  immersedFieldBinder, s2d, ib, isInside_ );
            
            base::nitsche::energyLHS<STBUP>( gradP_, surfaceQuadrature_, solver,
                                             immersedFieldBinder, ib, isInside_ );
            //base::nitsche::energyRHS2<STBUP>( gradP_, surfaceQuadrature_, solver,
            //                                  immersedFieldBinder, s2d, ib, isInside_ );

            // Finalise assembly
            solver.finishAssembly();

            // Solve
            solver.superLUSolve();

            // distribute results back to dofs
            base::dof::setDoFsFromSolver( solver, velocity_ );
            base::dof::setDoFsFromSolver( solver, pressure_ );

            // check convergence via solver norms
            const double newNorm = solver.norm();
            const double convCrit = (std::abs( prevNorm - newNorm ) / std::abs( newNorm ) );
            if ( (iter > 2) and (convCrit < tolerance) ) converged = true;
        
            // update storage of previous norm
            prevNorm = newNorm;
        
            iter++;
        }

        return iter;
    }


    //--------------------------------------------------------------------------
    void writeVTKFile( const std::string& baseName,
                       const unsigned step )
    {
        // VTK Legacy
        const std::string vtkFile = baseName
            + "F." + base::io::leadingZeros( step ) + ".vtk";
        
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh_ );

        base::io::vtk::writePointData( vtkWriter, mesh_, velocity_, "U" );
        base::io::vtk::writePointData( vtkWriter, mesh_, pressure_, "P" );
        
        {
            std::vector<double> distances;
            std::transform( levelSet_.begin(), levelSet_.end(),
                            std::back_inserter( distances ),
                            boost::bind( &LevelSet::getSignedDistance, _1 ) );
            vtkWriter.writePointData( distances.begin(), distances.end(), "distances" );
            
            //vtkWriter.writePointData( supportSizes_.begin(), supportSizes_.end(), "supports");
            // cannot write these, as there are more dofs than nodes :(
        }



        vtk.close();
    }


    //--------------------------------------------------------------------------
    template<typename SMESH, typename FFIELD>
    VecDim computeImmersedSurfaceForce( SMESH& surfaceMesh,
                                        FFIELD& surfForces, 
                                        const double viscosity,
                                        const double velocityWeight )
    {
        // generate mesh from immersed boundary
        typedef typename base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;
        SurfaceMesh immersedMesh;
        base::cut::generateSurfaceMesh<Mesh,Cell>( mesh_, cells_, immersedMesh );

        // bind a field around it
        typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,Velocity,Pressure> SurfaceFieldBinder;
        SurfaceFieldBinder immersedFieldBinder( immersedMesh, velocity_, pressure_ );

        typedef typename SurfaceFieldBinder::template TupleBinder<1,2>::Type   STBUP;
                
        VecDim sumOfForces = base::constantVector<dim>( 0. );
        
        typedef typename Field::template TupleBinder<1,2>::Type UP;
        typedef typename ::Robin<typename UP::Tuple> Traction;
        Traction traction( viscosity, velocityWeight );
        
        base::cut::ComputeSurfaceForces<SMESH,FFIELD,SurfaceQuadrature,
                                        typename STBUP::Tuple,Traction>
            computeSurfaceForces( surfaceMesh, surfForces, surfaceQuadrature_,
                                  levelSet_, traction, not isInside_ );

        typename SurfaceFieldBinder::FieldIterator first = immersedFieldBinder.elementsBegin();
        typename SurfaceFieldBinder::FieldIterator  last = immersedFieldBinder.elementsEnd();
        for ( ; first != last; ++first ) {

            sumOfForces +=
                computeSurfaceForces( STBUP::makeTuple( *first ) );
                
        }

        return sumOfForces;
    }



    //--------------------------------------------------------------------------
    Mesh&     accessMesh()     { return mesh_; }
    Velocity& accessVelocity() { return velocity_; }
    Pressure& accessPressure() { return pressure_; }

    std::vector<double>& supportSizes() { return supportSizesU_; }
    
private:
    const double isInside_;
    
    Mesh              mesh_;
    Quadrature        quadrature_;
    SurfaceQuadrature surfaceQuadrature_;
    Velocity          velocity_;
    Pressure          pressure_;
    Field             field_;

    const StressDivergence stressDivergence_;
    const Convection       convection_;
    const GradP            gradP_;
    const DivU             divU_;

    std::size_t numDoFsU_;
    std::size_t numDoFsP_;

    std::vector<LevelSet> levelSet_;
    std::vector<Cell>     cells_;

    std::vector<double>   supportSizesU_;
    std::vector<double>   supportSizesP_;
};
#endif
