//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   base/BoundaryValueProblem.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_boundaryvalueproblem_hpp
#define base_boundaryvalueproblem_hpp

//------------------------------------------------------------------------------
#include <fstream>
#include <boost/function.hpp>

#include <base/verify.hpp>
#include <base/shape.hpp>
#include <base/linearAlgebra.hpp>
#include <base/Field.hpp>

#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>

#include <base/fe/Basis.hpp>

#include <base/time/DummyMethod.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>

#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>
#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/asmb/NeumannForce.hpp>

#include <base/dof/generate.hpp>
#include <base/dof/constrainBoundary.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/Distribute.hpp>
#include <base/dof/setField.hpp>

#include <base/io/vtk/LegacyWriter.hpp>
//------------------------------------------------------------------------------
namespace base{

    //! \ingroup driver
    template<typename MESH,
             unsigned DOFSIZE, 
             unsigned FDEG     = MESH::Element::GeomFun::degree,
             typename TIMEINT  = base::time::NoTimeIntegration>
    class BoundaryValueProblem;

}

//------------------------------------------------------------------------------
/** Handle a (initial) boundary value problem.
 *  The prototype for this equation is
 *  \f[
 *       \rho \dot{u} - L u = f \quad x \in \Omega
 *  \f]
 *  with possible Dirichlet boundary conditions
 *  \f[
 *        u = \bar{u} \quad x \in \Gamma_D
 *  \f]
 *  and Neumann boundary conditions
 *  \f[
 *       t(u)  = \bar{t} x \in \Gamma_N
 *  \f]
 *  and initial condition
 *  \f[
 *       u = u_0  \quad x \in \Omega, t = 0
 *  \f]
 *
 *  \tparam MESH     Type of mesh for description of \f$ \Omega \f$
 *  \tparam DOFSIZE  Number of components per degree of freedom
 *  \tparam FDGE     Polynomial degree for the FE basis functions
 *  \tparam TIMEINT  Type of time integration
 */
template<typename MESH, unsigned DOFSIZE, unsigned FDEG, typename TIMEINT>
class base::BoundaryValueProblem
{
public:
    //! @name Template parameter
    //@{
    typedef MESH          Mesh;
    typedef TIMEINT       TimeIntegration;
    static const unsigned doFSize  = DOFSIZE;
    static const unsigned fieldDeg = FDEG;
    //@}

    //! @name Convenience parameters
    //@{
    static const unsigned    dim     = Mesh::Node::dim;
    static const base::Shape shape   = Mesh::Element::shape;
    static const unsigned    nHist   = TimeIntegration::numSteps;
    //@}

    //! @name FE Field types
    //@{
    typedef base::fe::Basis<shape,fieldDeg>                       FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>                    Field;
    typedef typename Field::DegreeOfFreedom                       DoF;
    typedef base::asmb::FieldBinder<Mesh,Field>                   FieldBinder;
    typedef typename FieldBinder::template TupleBinder<1,1>::Type UU;
    //@}

    //! Coordinate type
    typedef typename base::Vector<dim,double>::Type  VecDim;

    //! Field value type
    typedef typename base::Vector<doFSize>::Type     VecDoF;

    //! List of faces on the boundary of the mesh
    typedef base::mesh::MeshBoundary MeshBoundary;
    
    //! Surface mesh of the volume's mesh boundary
    typedef typename base::mesh::BoundaryMeshBinder<Mesh>::Type   BoundaryMesh;

    //! @name Surface FE field
    //@{
    typedef base::asmb::SurfaceFieldBinder<BoundaryMesh,Field> SurfaceFieldBinder;
    typedef typename SurfaceFieldBinder::template TupleBinder<1,1>::Type      SUU;
    //@}

    //! @name BC and body force functions
    //@{
    typedef boost::function<void  (const VecDim&, DoF*)>          DirichletFun;
    typedef boost::function<VecDoF(const VecDim&)>                BodyForceFun;
    typedef boost::function<VecDoF(const VecDim&, const VecDim&)> SurfaceForceFun;
    //@}

    //--------------------------------------------------------------------------
    /** Pass mesh reference and initialise the local data.
     *  The geometry triangulation, which forms the foundation of the FE
     *  approximation, is represented by the mesh object whose reference is
     *  passed here. The mesh object has to be set before calling this
     *  constructor!
     *
     *  In this constructor the degrees of freedom which represent the FE field
     *  are constructed and a face-list of the elements adjacent to the domain
     *  boundary is generated. The state variables are given their initial
     *  values.
     *
     *  \param[in] mesh          FE triangulation of the computational domain
     *  \param[in] density       Material density, needed for time integration
     */
    BoundaryValueProblem( Mesh& mesh, double density = 1. )
        : mesh_(     mesh ),
          density_(  density )
    {
        // Check if mesh is non-empty
        VERIFY_MSG( std::distance( mesh_.elementsBegin(), mesh_.elementsEnd() ),
                    "Mesh is empty - a proper mesh is needed here " );
        
        // generate degrees of freedom first
        base::dof::generate<FEBasis>( mesh_, field_ );
        
        // extract boundary faces from mesh, necessary for Dirichlet constraints
        meshBoundary_.create( mesh_.elementsBegin(), mesh_.elementsEnd() );

        // set state indicators to their initial values
        this -> initiateStates_();
    }

    /** Make a FE triangulation out of the boundary of the domain mesh.
     *  Such a mesh is needed for the application of surface forces.
     */
    void generateBoundaryMesh()
    {
        if ( not haveBoundaryMesh_ ) {
            // only necessary for Neumann boundary conditions
            base::mesh::generateBoundaryMesh( meshBoundary_.begin(),
                                              meshBoundary_.end(),
                                              mesh_, boundaryMesh_ );
            haveBoundaryMesh_ = true;
        }
    }

    /** Set the initial values by means of a given function.
     *  For time-dependent problems, often an initial state \f$ u_0 \f$ is
     *  given from which the problem develops. Such state is defined by
     *  function given by the caller.
     *  \param[in] initialValueFun Function for the initial state
     */
    void setInitialValue( DirichletFun initialValueFun )
    {
        base::dof::setField( mesh_, field_, initialValueFun );
        base::dof::pushHistory( field_ );
    }
 
    /** Set the boundary values by means of a given function.
     *  For Dirichlet boundary conditions the datum \f$ \bar{u} \f$ is
     *  given. The function is given by the caller.
     *  \note This function needs to be called \e before numbering the
     *        degrees of freedom
     *
     *  \param[in] diriFun Function for the Dirichlet datum
     */
    void applyDirichletConstraints( DirichletFun diriFun )
    {
        // Important: first apply constraints, then number doFs
        VERIFY_MSG(
            (doFsAreNumbered_ == false),
            "Constraints have to be applied before numbering the DoFs" );
        
        base::dof::constrainBoundary<FEBasis>(
            meshBoundary_.begin(), meshBoundary_.end(),
            mesh_, field_, diriFun );
    }

    /** Assign a unique number to each degree of freedom component.
     *  Starting from a given number (defaults to zero), every unconstrained
     *  Degree of freedom receives a unique number.
     *  \note This function must not be called before application of the
     *        Dirichlet constraints in applyDirichletConstraints()
     *
     *  \param[in] init Initial number for the DoFs
     *  \return         The number of DoFs numbered by this function
     */
    std::size_t numberDoFs( const std::size_t init = 0 )
    {
        const std::size_t totalNumber =
            base::dof::numberDoFsConsecutively( field_.doFsBegin(),
                                                field_.doFsEnd(), init );
        doFsAreNumbered_ = true;
        
        return totalNumber - init;
    }

    /** In case of multi-threaded assembly, the sparsity pattern of the
     *  matrix is registered in the solver.
     *  \tparam SOLVER Type of solver
     *  \param[in,out] solver The solver object
     */
    template<typename SOLVER>
    void registerInSolver( SOLVER& solver )
    {
        VERIFY_MSG( doFsAreNumbered_, "DoFs need to be numbered first" );
        FieldBinder fb( mesh_, field_ );
        solver.template registerFields<UU>( fb );
    }

    /** Assemble the element stiffness matrix of the Laplace operator.
     *  The kernel of this p.d.e. is in fact the bilinear form stemming
     *  from the linearisation of the term
     *  \f[
     *       - \nabla \cdot ( D(u) \nabla u )
     *  \f]
     *  and is not necessarily the same as the one from the Laplace operator,
     *  but in most applications this is the case. In all cases where
     *  \f$ D \f$ is \e not a function of \f$ u \f$, the Laplace situation
     *  is recovered.
     *  \tparam KERNEL Type of integration kernel defining the bilinear form
     *  \tparam QUAD   Type of quadrature
     *  \tparam SOLVER Type of equation solver
     *  \param[in] kernel      The integrand
     *  \param[in] quadrature  Numerical integration object
     *  \param[in] solver      Equation solver object
     *  \param[in] incremental Type of analysis
     */
    template<typename KERNEL, typename QUAD, typename SOLVER>
    void assembleBilinearForm( const KERNEL& kernel,
                               const QUAD& quadrature,
                               SOLVER& solver,
                               const bool incremental = false )
    {
        VERIFY_MSG( doFsAreNumbered_, "DoFs need to be numbered first" );
        FieldBinder fb( mesh_, field_ );
        base::asmb::stiffnessMatrixComputation<UU>( quadrature, solver, fb,
                                                    kernel, incremental );

        // residual forces
        base::asmb::computeResidualForces<UU>( quadrature, solver,
                                               fb, kernel );
    }

    /** Assemble nodal forces due to a domain force term.
     *  \note Here the function has to have the format
     *  \code{.cpp}
     *  VecDoF forceFun( const VecDim& x ) { };
     *  \endcode
     *       Other variants have to be handled outside of this class.
     *
     *  \tparam QUAD   Type of quadrature
     *  \tparam SOLVER Type of equation solver
     *  \param[in] quadrature Numerical integration object
     *  \param[in] solver     Equation solver object
     *  \param[in] forceFun   Function representing the domain force.
     */
    template<typename QUAD, typename SOLVER>
    void applyBodyForce( const QUAD&   quadrature,
                         SOLVER&       solver,
                         BodyForceFun  forceFun )
    {
        VERIFY_MSG( doFsAreNumbered_, "DoFs need to be numbered first" );
        FieldBinder fb( mesh_, field_ );
        base::asmb::bodyForceComputation<UU>( quadrature, solver, fb,
                                              forceFun );
    }

    /** Assemble nodal forces due to a surface force term.
     *  \note Here the function has to have the format
     *  \code{.cpp}
     *  VecDoF forceFun( const VecDim& x, const VecDim& normal ) { };
     *  \endcode
     *       Other variants have to be handled outside of this class.
     *
     *  \note First the function generateBoundaryMesh() has to be called
     *        in order to have a proper surface mesh description.
     *
     *  \tparam SQUAD  Type of surface quadrature
     *  \tparam SOLVER Type of equation solver
     *  \param[in] surfaceQuadrature Numerical integration object
     *  \param[in] solver            Equation solver object
     *  \param[in] neumannForceFun   Function representing the boundary force
     */
    template<typename SQUAD, typename SOLVER>
    void applyBoundaryForces( const SQUAD&    surfaceQuadrature,
                              SOLVER&         solver,
                              SurfaceForceFun neumannForceFun )
    {
        VERIFY_MSG( doFsAreNumbered_, "DoFs need to be numbered first" );
        VERIFY_MSG( haveBoundaryMesh_,
                    "Cannot apply Neumann BC without boundary mesh" );

        SurfaceFieldBinder sfb( boundaryMesh_, field_ );
        base::asmb::neumannForceComputation<SUU>( surfaceQuadrature,
                                                  solver, sfb,
                                                  neumannForceFun );
    }

    /** Assemble the terms stemming from a time integration method.
     */
    template<typename KERNEL, typename QUAD, typename SOLVER>
    void applyTimeIntegrator( const unsigned stepNum,
                              const double   stepSize,
                              const KERNEL&  kernel, 
                              const QUAD&    quadrature,
                              SOLVER&        solver,
                              const unsigned iter = 0 )
    {
        VERIFY_MSG( doFsAreNumbered_, "DoFs need to be numbered first" );
        FieldBinder fb( mesh_, field_ );

        // compute inertia terms, d/dt, due to time integration
        base::time::computeInertiaTerms<UU,TimeIntegration>(
            quadrature, solver, fb, stepSize, stepNum, density_, iter > 0 );

        // compute history of residual forces due to time integration
        base::time::computeResidualForceHistory<UU,TimeIntegration>(
            kernel, quadrature, solver, fb, stepNum );
    }

    /** After solving, pass back the values to the Degrees of freedom.
     */
    template<typename SOLVER>
    void setDoFsFromSolver( const SOLVER& solver)
    {
        base::dof::setDoFsFromSolver( solver, field_ );
        base::dof::pushHistory( field_ );
    }

    /** After solving, pass back the values to the Degrees of freedom.
     */
    template<typename SOLVER>
    void addToDoFsFromSolver( const SOLVER& solver)
    {
        base::dof::addToDoFsFromSolver( solver, field_ );
    }

    //--------------------------------------------------------------------------
    //! @name Accessors
    //@{
    Mesh&         getMesh()         { return mesh_;  }
    Field&        getField()        { return field_; }
    BoundaryMesh& getBoundaryMesh() { return boundaryMesh_; }
    //@}

    //--------------------------------------------------------------------------
    //! Write field state to a VTK file
    void writeVTKFile( const std::string& baseName,
                       const std::string& fieldName = "u",
                       const std::string& gradName  = "gradU" ) const
    {
        // open file
        const std::string vtkFile = baseName + ".vtk";
        std::ofstream vtk( vtkFile.c_str() );
        base::io::vtk::LegacyWriter vtkWriter( vtk );

        // write mesh
        vtkWriter.writeUnstructuredGrid( mesh_ );
        base::io::vtk::writePointData( vtkWriter, mesh_, field_, fieldName );
            
        // write solution gradient
        const typename base::Vector<dim>::Type xi =
            base::ShapeCentroid<Mesh::Element::shape>::apply();
        base::io::vtk::writeCellData( vtkWriter, mesh_, field_,
                                      boost::bind(
                                          base::post::template
                                          evaluateFieldGradient<
                                          typename Mesh::Element,
                                          typename Field::Element>, _1, _2, xi ),
                                      gradName );
        vtk.close();
    }

private:
    //! Initial value of the flags
    void initiateStates_()
    {
        haveBoundaryMesh_ = false;
        doFsAreNumbered_  = false;
    }
    
private:
    //! Reference to the domain triangulation
    Mesh&        mesh_;
    //! Field representing the FE approximation
    Field        field_;
    //! The boundary of the mesh is a list of element faces
    MeshBoundary meshBoundary_;
    //! A mesh made from the mesh-boundary 
    BoundaryMesh boundaryMesh_;

    //! Material density (for dynamic problems only)
    const double density_;

    //! @name State indicators
    //@{
    
    /** For the application of Neumann boundary condition, one has to
        integrated over the boundary of the domain and, moreover, the
        outward normal vector is needed. This flag ensures that the
        user first calls generateBoundaryMesh() before using surface-
        related functions.
     */
    bool haveBoundaryMesh_;

    /** This flag serves two purposes
     *  1. Dirichlet constraints must not be set after numbering of the
     *     degrees of freedom
     *  2. All assembly-related functionality requires numbered DoFs
     */
    bool doFsAreNumbered_;
    //@}
};

#endif
