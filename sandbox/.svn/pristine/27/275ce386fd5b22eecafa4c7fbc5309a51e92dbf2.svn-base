#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <base/Unstructured.hpp>
#include <base/Quadrature.hpp>

#include <base/io/Format.hpp>
#include <base/io/smf/Reader.hpp>
#include <base/io/vtk/LegacyWriter.hpp>

#include <base/Field.hpp>
#include <base/dof/numbering.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/setField.hpp>
#include <base/fe/Basis.hpp>
#include <base/dof/Distribute.hpp>

#include <base/asmb/StiffnessMatrix.hpp>
#include <base/asmb/ForceIntegrator.hpp>
#include <base/asmb/NeumannForce.hpp>
#include <base/asmb/FieldBinder.hpp>
#include <base/asmb/BodyForce.hpp>
#include <base/solver/Eigen3.hpp>
#include <base/time/BDF.hpp>
#include <base/time/AdamsMoulton.hpp>
#include <base/time/ReactionTerms.hpp>
#include <base/time/ResidualForceHistory.hpp>
#include <base/time/derivative.hpp>

#include <base/post/evaluateField.hpp>
//#include <base/post/ErrorNorm.hpp>


//------------------------------------------------------------------------------
template<typename FIELDTUPLE>
class Kernel
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    typedef typename FieldTuple::TestElement  TestElement;
    typedef typename FieldTuple::TrialElement TrialElement;
    //@}

    //! Local coordinate
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim  LocalVecDim;

    Kernel( const double a ) : a_( a ) { }

    //--------------------------------------------------------------------------
    void tangentStiffness( const FieldTuple&  fieldTuple, 
                           const LocalVecDim& xi,
                           const double       weight,
                           base::MatrixD&     matrix ) const
    {
        matrix( 0, 0 ) += a_ * weight;
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST>
    void residualForceHistory( const FieldTuple&   fieldTuple,
                               const LocalVecDim&  xi,
                               const double        weight,
                               base::VectorD&      vector ) const
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        //const TestElement*  testEp  = fieldTuple.testElementPtr();
        const TrialElement* trialEp = fieldTuple.trialElementPtr();
        
        // Evalute temperature gradient
        const typename base::Vector<1>::Type y
            = base::post::evaluateFieldHistory<HIST>( geomEp, trialEp, xi );
        
        vector[0] += a_ * y[0] * weight;
    }

private:
    const double a_;
};

//------------------------------------------------------------------------------
template<unsigned DIM>
base::Vector<1,double>::Type
forceFun( const typename base::Vector<DIM,double>::Type& x,
          const double time )
{
    base::Vector<1,double>::Type result;
    result[0] = -std::sin( time ) + 5. * std::cos( 5. * time );
    return result;
}


//------------------------------------------------------------------------------
// initial condition
template<unsigned DIM, typename DOF>
void setField( const typename base::Vector<DIM>::Type& x, DOF* doFPtr )
{
    doFPtr -> setValue( 0, 1.0 ); 
    doFPtr -> pushHistory();
}

// analytic solution: exp( -a * t )
double analytic( const double time, const double a )
{
    return std::exp( -a * time );
}

// derivative: -a * exp( -a * t )
double analyticDerivative( const double time, const double a )
{
    return -a * analytic( time, a );
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    if ( argc != 4 ) {

        std::cout << "Usage: " << argv[0] << "  numSteps  deltaT  a  \n"
                  << '\n'
                  << "Solves  dot{y} + a y = 0, y(0) = 1 \n";
        return 0;
    }

    const unsigned tiOrder  = 3;
    typedef base::time::BDF<tiOrder> MSM;
    //typedef base::time::AdamsMoulton<tiOrder> MSM;
    
    // user input
    const unsigned numSteps = boost::lexical_cast<unsigned>( argv[1] );
    const double   deltaT   = boost::lexical_cast<double>(   argv[2] );
    const double   a        = boost::lexical_cast<double>(   argv[3] );

    // FE Stuff
    const unsigned geomDeg  = 1;
    const unsigned fieldDeg = 0;
    const base::Shape shape = base::LINE;
    const unsigned nHist    = MSM::numSteps; 
    const unsigned dim      = base::ShapeDim<shape>::value;

    // Mesh
    typedef base::Unstructured<shape,geomDeg>     Mesh;

    Mesh mesh;
    {
        std::stringstream buffer;
        buffer << "! elementShape     line \n"
               << "! elementNumPoints 2    \n"
               << "2  1 \n"
               << "0. 0. 0. \n "
               << "1. 0. 0. \n "
               << "0  1 \n";

        base::io::smf::Reader<Mesh> smfReader;
        smfReader( mesh, buffer ); 
    }

    // quadrature objects for volume and surface
    const unsigned kernelDegEstimate = 1;
    typedef base::Quadrature<kernelDegEstimate,shape> Quadrature;
    Quadrature quadrature;

    // Create a field
    const unsigned    doFSize = 1;
    typedef base::fe::Basis<shape,fieldDeg>           FEBasis;
    typedef base::Field<FEBasis,doFSize,nHist>        Field;
    Field field;

    // generate DoFs from mesh
    base::dof::generate<FEBasis>( mesh, field );

    // Number the degrees of freedom
    const std::size_t numDoFs =
        base::dof::numberDoFsConsecutively( field.doFsBegin(), field.doFsEnd() );

    // Definition of the field combination
    typedef base::asmb::FieldBinder<Mesh,Field> FieldBinder;
    typedef FieldBinder::TupleBinder<1,1>::Type TupleBinder;
    FieldBinder fieldBinder( mesh, field );

    typedef Kernel<TupleBinder::Tuple> Kernel;
    Kernel kernel( a );

    // set initial conditions
    base::dof::setField( mesh, field,
                         boost::bind( &setField<dim,Field::DegreeOfFreedom>, _1, _2 ) );


    Mesh::Element*  geomEPtr  = mesh.elementPtr( 0 );
    Field::Element* fieldEPtr = field.elementPtr( 0 );
    Field::Element::FEFun::VecDim xi =
        base::ShapeCentroid<shape>::apply();

    double error0 = 0.;
    double error1 = 0.;


    for ( unsigned n = 0; n < numSteps; n ++ ) {

        const double time = (n+1) * deltaT;

        // Create a solver object
        typedef base::solver::Eigen3           Solver;
        Solver solver( numDoFs );

        //------------------------------------------------------------------
        base::asmb::stiffnessMatrixComputation<TupleBinder>( quadrature, solver, 
                                                             fieldBinder, kernel );

        base::time::computeInertiaTerms<TupleBinder,MSM>( quadrature, solver,
                                                          fieldBinder, deltaT, n, 1.0 );


        base::time::computeResidualForceHistory<TupleBinder,MSM>( kernel, quadrature,
                                                                  solver, fieldBinder, n );

        // Body force
        base::asmb::bodyForceComputation<TupleBinder>( quadrature, solver, fieldBinder,
                                                       boost::bind( &forceFun<dim>, _1, time ) );


        // Finalise assembly
        solver.finishAssembly();

        // Solve
        solver.choleskySolve();
            
        // distribute results back to dofs
        base::dof::setDoFsFromSolver( solver, field );

        const double aux1= std::abs( std::cos(time) + std::sin(5. * time) - 
                                     ( base::post::evaluateField( geomEPtr, fieldEPtr, xi ) )[0] );
        error0 += (aux1 * aux1) * deltaT;

        const double aux2 = std::abs( -std::sin(time) + 5. * std::cos(5.*time) - 
                                           ( base::time::evaluateTimeDerivative( geomEPtr, fieldEPtr, xi,
                                                                                 deltaT, n ) )[0] );
        error1 += (aux2 * aux2) * deltaT;

#if 0
        std::cout << time << " "
                  << //std::abs( analytic( time, a ) -
            (//std::abs( std::cos(time) + std::sin(5. * time) - 
                      ( base::post::evaluateField( geomEPtr, fieldEPtr, xi ) )[0] )
                  << "  "
                  << //std::abs( analyticDerivative( time, a ) -
            (//std::abs( -std::sin(time) + 5. * std::cos(5. * time) - 
                      ( base::time::evaluateTimeDerivative( geomEPtr, fieldEPtr, xi,
                                                            deltaT, n ) )[0] )
                  << '\n';
#endif
        
        // push history
        std::for_each( field.doFsBegin(), field.doFsEnd(),
                       boost::bind( &Field::DegreeOfFreedom::pushHistory, _1 ) );
        

        
        // Finished time steps
        //--------------------------------------------------------------------------
    }

    std::cout << deltaT << "  " << error0 << "  " << error1 << '\n';
    
}
