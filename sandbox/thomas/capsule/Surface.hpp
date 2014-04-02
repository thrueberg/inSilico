#ifndef surface_hpp
#define surface_hpp

#include <fstream>
#include <string>

#include <base/io/smf/Reader.hpp>
#include <base/shape.hpp>
#include <base/Unstructured.hpp>
#include <base/Quadrature.hpp>
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
#include <solid/Stress.hpp>

#include <surf/HyperElasticMembrane.hpp>
#include <surf/Skalak.hpp>
#include <surf/NeoHookean.hpp>
#include <surf/Moments.hpp>

#include <base/time/BDF.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>


//------------------------------------------------------------------------------
template<unsigned DIM, unsigned GDEG=1, unsigned FDEG=GDEG>
class Surface
{
public:
    //! @name Template parameter
    //@{
    static const unsigned dim      = DIM;
    static const unsigned geomDeg  = GDEG;
    static const unsigned fieldDeg = FDEG;
    //@}

    //! Time integration hard-wired to Euler-backward
    static const unsigned tiOrder = 1;
    typedef base::time::BDF<tiOrder> MSM;
    static const unsigned nHist = MSM::numSteps;

    //! Use always DIM-1-simplex elements
    static const base::Shape shape = base::SimplexShape<dim-1>::value;

    //! Type of mesh
    typedef base::Unstructured<shape,geomDeg,dim>    Mesh;

    typedef typename Mesh::Node::VecDim  VecDim;

    //! Type of quadrature
    static const unsigned kernelDegEstimate = 3;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    
    //! @name Displacement Field
    //@{
    static const unsigned    doFSize = dim;
    typedef base::fe::Basis<shape,fieldDeg>        FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>     Field;
    //@}

    //! @name Field binding
    //@{
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    typedef typename FieldBinder::template TupleBinder<1,1>::Type FTB;
    //@}

    // material object
    typedef surf::Skalak Energy;
    typedef surf::HyperElasticMembrane<Energy> Material;
    typedef solid::HyperElastic<Material,typename FTB::Tuple> HyperElastic;

    typedef base::kernel::Mass<typename FTB::Tuple> Mass;

    // Solver type
    typedef base::solver::Eigen3           Solver;

    //! Constructor
    Surface( const double A, const double B, const double gamma = 0. )
        : fieldBinder_( mesh_, displacement_ ),
          energy_( A, B ),
          material_(     energy_ ),
          hyperElastic_( material_ ),
          mass_(         gamma )
    { }

    //! Read mesh from file
    void readMesh( const std::string& meshFile )
    {
        std::ifstream smf( meshFile.c_str() );
        base::io::smf::readMesh( smf, mesh_ );
        smf.close();
        
        // generate DoFs from mesh
        base::dof::generate<FEBasis>( mesh_, displacement_ );
        base::dof::generate<FEBasis>( mesh_, forces_ );
        base::dof::generate<FEBasis>( mesh_, velocity_ );
    }

    //! Generate dofs and assign dof numbers
    std::size_t numberDoFs( )
    {
        // Number the degrees of freedom
        numDoFs_ =
            base::dof::numberDoFsConsecutively( displacement_.doFsBegin(),
                                                displacement_.doFsEnd() );

        base::dof::numberDoFsConsecutively( forces_.doFsBegin(),
                                            forces_.doFsEnd() );
        base::dof::numberDoFsConsecutively( velocity_.doFsBegin(),
                                            velocity_.doFsEnd() );
        return numDoFs_;
    }

    //!
    void assembleForcesFromField( Solver& solver )
    {
        static const unsigned doFSize = Field::DegreeOfFreedom::size;
        
        typename Field::DoFPtrConstIter dIter = forces_.doFsBegin();
        typename Field::DoFPtrConstIter dEnd  = forces_.doFsEnd();

        for ( ; dIter != dEnd; ++dIter ) {

            // get dof status
            std::vector<base::dof::DoFStatus> status;
            (*dIter) -> getStatus(  std::back_inserter( status ) );

            // get dof IDs
            std::vector<std::size_t>          ids;
            (*dIter) -> getIndices( std::back_inserter( ids ) );

            // get dof constraints
            std::vector< std::pair<unsigned,
                                   std::vector<std::pair<double,
                                                         std::size_t> >
                                   > > constraints;
            
            for ( unsigned d = 0; d < doFSize; d++ ) {
                if ( status[d] == base::dof::CONSTRAINED ) {

                    // get pairs of weight and global DoF ID
                    std::vector<std::pair<double,std::size_t> > weightedDoFIDs;
                    (*dIter) -> getConstraint(d) -> getWeightedDoFIDs( weightedDoFIDs );

                    // store pair of local ID and weighted dof IDs
                    constraints.push_back( std::make_pair( d, weightedDoFIDs ) );
                }
            }

            // get values from dof
            base::VectorD fVec = base::VectorD::Zero( doFSize );
            for ( unsigned d = 0; d < doFSize; d++ )
                fVec[d] = ( (*dIter) -> getValue(d) );

            // assemble to solver
            base::asmb::assembleForces( fVec, status, ids, constraints, solver );
            
        }
        
        return;
    }

    void computeSurfaceVelocity( const double dt ) 
    {
        typename Field::DoFPtrConstIter dIter = displacement_.doFsBegin();
        typename Field::DoFPtrConstIter dEnd  = displacement_.doFsEnd();
        typename Field::DoFPtrIter vIter = velocity_.doFsBegin();
        for ( ; dIter != dEnd; ++dIter, ++vIter ) {

            for ( unsigned d = 0; d < dim; d++ ) {
                const double un2 = (*dIter) ->template getHistoryValue<0>( d );
                const double un1 = (*dIter) ->template getHistoryValue<1>( d );
                (*vIter) -> setValue( d, (un2 - un1) / dt );
            }
            
        }
        return;
    }

    //!
    unsigned findEquilibrium( const double dt,
                              const unsigned step, 
                              const unsigned maxIter,
                              const double tolerance )
    {
        
        unsigned iter = 0;
        while ( iter < maxIter ) {

            // Create a solver object

            Solver solver( numDoFs_ );

            
            base::asmb::computeResidualForces<FTB>( quadrature_, solver,
                                                    fieldBinder_, hyperElastic_ );
            
            // Compute element stiffness matrices and assemble them
            base::asmb::stiffnessMatrixComputation<FTB>( quadrature_, solver,
                                                         fieldBinder_, hyperElastic_,
                                                         iter > 0 );

            // assemble forces from field
            this -> assembleForcesFromField( solver );

            // compute inertia terms, d/dt, due to time integration
            base::time::computeReactionTerms<FTB,MSM>( mass_, quadrature_, solver,
                                                       fieldBinder_, dt, step,
                                                       iter > 0 );


            // Finalise assembly
            solver.finishAssembly();

            // norm of residual 
            const double conv1 = solver.norm();


            // convergence via residual norm
            const double matFactor = energy_.energy( 1., 1. );
            if ( conv1 < tolerance * matFactor ) { // note the tolerance multiplier
                break;
            }

            // Solve
            solver.choleskySolve();
            
            // distribute results back to dofs
            base::dof::addToDoFsFromSolver( solver, displacement_, iter > 0 );

            // norm of displacement increment
            const double conv2 = solver.norm();
            
            // convergence via increment
            if ( conv2 < tolerance ) break;

            iter++;

        }
        // Finished non-linear iterations
        //----------------------------------------------------------------------

        return iter;
    }


    //------------------------------------------------------------------------------
    void writeVTKFile( const std::string& baseName,
                       const unsigned     step )
    {
        // create file name with step number
        const std::string vtkFile =
            baseName + "." + base::io::leadingZeros( step ) + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );
        vtkWriter.writeUnstructuredGrid( mesh_ );

        base::io::vtk::writePointData( vtkWriter, mesh_, displacement_, "U" );
        base::io::vtk::writePointData( vtkWriter, mesh_, forces_, "F" );
        base::io::vtk::writePointData( vtkWriter, mesh_, velocity_, "V" );

        base::io::vtk::writeCellData( vtkWriter, mesh_, displacement_, 
                                      boost::bind( solid::cauchy<typename Mesh::Element,
                                                   typename Field::Element,
                                                   Material>,
                                                   _1, _2, material_ ), "sigma" );

        vtk.close();
    }

    //------------------------------------------------------------------------------
    Mesh&  accessMesh() {          return mesh_; }
    Field& accessDisplacements() { return displacement_; }
    Field& accessVelocity()      { return velocity_; }
    Field& accessForces()        { return forces_; }
    const Quadrature& accessQuadrature() const { return quadrature_; }

    //------------------------------------------------------------------------------
    void currentConfiguration( Mesh& currentConfig ) const
    {
        VERIFY_MSG( std::distance( currentConfig.nodesBegin(), currentConfig.nodesEnd() )
                    == 0, "Mesh node container has to be empty!" );
        VERIFY_MSG( std::distance( currentConfig.elementsBegin(), currentConfig.elementsEnd() )
                    == 0, "Mesh element container has to be empty!" );

        // copy reference configuration into the mesh
        currentConfig.append( mesh_ );

        // get nodal displacements
        std::vector<VecDim> nodalDisp;
        base::post::evaluateAtNodes( mesh_, displacement_, nodalDisp );

        // add displacements to reference configuration
        typename Mesh::NodePtrIter newNIter = currentConfig.nodesBegin();
        typename Mesh::NodePtrIter newNEnd  = currentConfig.nodesEnd();
        typename std::vector<VecDim>::iterator dIter = nodalDisp.begin();
        for ( ; newNIter != newNEnd; ++newNIter, ++dIter )
        {
            const VecDim u = *dIter;

            VecDim X;
            (*newNIter) -> getX( &(X[0]) );

            const VecDim x = X + u;
            (*newNIter) -> setX( &(x[0]) );
        }
        return;
    }

    //------------------------------------------------------------------------------
    void fieldCopy( Field& copy ) const
    {
        // generate DoFs from mesh
        base::dof::generate<FEBasis>( mesh_, copy );
    }

    //------------------------------------------------------------------------------
    double enclosedVolume() const
    {
        Mesh currentConfig;
        this -> currentConfiguration( currentConfig );

        Field dummy;
        this -> fieldCopy( dummy );
        FieldBinder tmp( currentConfig, dummy );

        return surf::enclosedVolume<FTB>( quadrature_, tmp );
    }

    //--------------------------------------------------------------------------
    VecDim moment() const
    {
        Mesh currentConfig;
        this -> currentConfiguration( currentConfig );

        Field dummy;
        this -> fieldCopy( dummy );
        FieldBinder tmp( currentConfig, dummy );

        return surf::enclosedVolumeMoment<FTB>( quadrature_, tmp );
    }


    //------------------------------------------------------------------------------
    void updateHistory()
    {
        // push history
        std::for_each( displacement_.doFsBegin(),
                       displacement_.doFsEnd(),
                       boost::bind( &Surface::Field::DegreeOfFreedom::pushHistory, _1 ) );
    }


private:
    Mesh               mesh_;
    const Quadrature   quadrature_;
    Field              displacement_;
    Field              velocity_;
    Field              forces_;
    FieldBinder        fieldBinder_;
    const Energy       energy_;
    const Material     material_;
    const HyperElastic hyperElastic_;
    const Mass         mass_;

    std::size_t numDoFs_;
};



#endif
