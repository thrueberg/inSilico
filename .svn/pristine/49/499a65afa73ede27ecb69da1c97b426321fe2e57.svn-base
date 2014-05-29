//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ForceIntegrator.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_forceintegrator_hpp
#define base_asmb_forceintegrator_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
// boost includes
#include <boost/bind.hpp>
#include <boost/function.hpp>
// base includes
#include <base/linearAlgebra.hpp>
// base/asmb includes
#include <base/asmb/collectFromDoFs.hpp>
#include <base/asmb/assembleForces.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{
        
        template<typename QUAD, typename SOLVER, typename FIELDTUPLE>
        class ForceIntegrator;

        //--------------------------------------------------------------------------
        // Convenience function to compute the residual forces
        template<typename FIELDTUPLEBINDER,
                 typename QUADRATURE, typename SOLVER,
                 typename FIELDBINDER, typename KERNEL>
        void computeResidualForces( const QUADRATURE& quadrature,
                                    SOLVER&           solver,
                                    const FIELDBINDER& fieldBinder,
                                    const KERNEL& kernelObj )
        {
            typedef ForceIntegrator<QUADRATURE,SOLVER,
                                    typename FIELDTUPLEBINDER::Tuple>
                ForceIntegrator;

            typename ForceIntegrator::ForceKernel
                residualForce = boost::bind( &KERNEL::residualForce, 
                                             &kernelObj,
                                             _1, _2, _3, _4 );

            // Note the -1 for moving the forces to the RHS of the system
            ForceIntegrator forceInt( residualForce, quadrature, solver, -1.0 );

            // apply to all elements
            typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
            typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
            for ( ; iter != end; ++iter ) {
                forceInt( FIELDTUPLEBINDER::makeTuple( *iter ) );
            }

            //const std::size_t numElements = std::distance( iter, end );
            //#pragma omp parallel for
            //for ( std::size_t e = 0; e < numElements; e++ ) {
            //    forceInt( FIELDTUPLEBINDER::makeTuple( fieldBinder.elementPtr( e ) ) );
            //}

            return;
        }

        
    } // namespace asmb    
} // namespace base

//------------------------------------------------------------------------------
/** Computation of the 'nodal' forces due to a given element  force function.
 *  The contribution of a force term to the system RHS comes from
 *  the linear form 
 *  \f[
 *        L( \phi ) =  \int_{\Omega} f(\phi) d x, 
 *  \f]
 *  where \f$ f \f$ is the force function. In the implementation, the evaluation
 *  of \f$ f \f$ will be represented by means of an element pointer and local
 *  quadrature point with weight. Finally, this term will be approximated with a
 *  given quadrature rule.
 *
 *  \tparam QUAD        Type of quadrature 
 *  \tparam SOLVER      Type of system solver with RHS storage
 *  \tparam FIELDTUPLE  Type of tuple holding field element pointers
 */
template<typename QUAD, typename SOLVER, typename FIELDTUPLE>
class base::asmb::ForceIntegrator
    : public boost::function<void( const FIELDTUPLE& )>
{
public:
    //! @name Template parameter
    //@{
    typedef QUAD        Quadrature;
    typedef SOLVER      Solver;
    typedef FIELDTUPLE  FieldTuple;
    //@}

    typedef typename FieldTuple::TestElementPtr   TestElementPtr;
    
    //! Type of forcing function
    typedef boost::function< void( const FieldTuple&,
                                   const typename Quadrature::VecDim&,
                                   const double,
                                   base::VectorD& )>  ForceKernel;

    //--------------------------------------------------------------------------
    //! Constructor setting all references
    ForceIntegrator( ForceKernel&         forceKernel,
                     const Quadrature&    quadrature,
                     Solver&              solver,
                     const double factor = 1.0 )
        : forceKernel_( forceKernel ),
          quadrature_(  quadrature ),
          solver_(      solver ),
          factor_(      factor )
    { }

    //--------------------------------------------------------------------------
    void operator()( const FieldTuple& fieldTuple )
    {
        // extract test and trial elements from tuple
        TestElementPtr  testEp  = fieldTuple.testElementPtr();
        
        // dof activities and IDs
        std::vector<base::dof::DoFStatus> doFStatus;
        std::vector<std::size_t> doFIDs;
        std::vector<base::number> prescribedValues; // placeholder
        std::vector<
            std::pair<unsigned,std::vector<std::pair<base::number,std::size_t> > >
            > constraints;

        //  get all dof data
        const bool doSomething = 
            base::asmb::collectFromDoFs( testEp, doFStatus, doFIDs,
                                         prescribedValues, constraints, false );

        // if no dof is ACTIVE or CONSTRAINED, just return
        if ( not doSomething ) return;
        
        // Compute the element contribution to all its dofs
        base::VectorD forceVec = base::VectorD::Zero( doFIDs.size() );

        // apply quadrature
        quadrature_.apply( forceKernel_, fieldTuple, forceVec );
        
        // apply factor
        forceVec *= factor_;
        
        // Assemble to solver
        base::asmb::assembleForces( forceVec, doFStatus, doFIDs, constraints,  solver_ );
        
        return;
    }

private:
    ForceKernel&         forceKernel_; //!< Function producing the force term
    const Quadrature&    quadrature_;  //!< Quadrature
    Solver&              solver_;      //!< Solver with RHS storage
    const double         factor_;      //!< Scalar multiplier
};


#endif

