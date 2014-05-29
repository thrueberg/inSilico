//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   StiffnessMatrix.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_asmb_stiffnessmatrix_hpp
#define base_asmb_stiffnessmatrix_hpp

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
#include <base/auxi/EqualPointers.hpp>
#include <base/auxi/parallel.hpp>
// base/asmb includes
#include <base/asmb/collectFromDoFs.hpp>
#include <base/asmb/assembleMatrix.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        template<typename QUAD, typename SOLVER, typename FIELDTUPLE>
        class StiffnessMatrix;

        //----------------------------------------------------------------------
        /** Convenience function for the computation of the system's stiffness
         *  matrix.
         *  This function combines the steps of (i) creating a kernel function
         *  object, (ii) creating an object of StiffnessMatrix and (iii) applying
         *  the latter to the range of elements in the given mesh and fields.
         *  \tparam FIELDTUPLEBINDER Binding of the right field tuple
         *  \tparam QUADRATURE   Type of quadrature
         *  \tparam SOLVER       Type of solver
         *  \tparam FIELDBINDER  Type of field compound
         *  \tparam KERNEL       Type of object with kernel function implementation
         */
        template<typename FIELDTUPLEBINDER,
                 typename QUADRATURE, typename SOLVER, typename FIELDBINDER,
                 typename KERNEL>
        void stiffnessMatrixComputation( const QUADRATURE& quadrature,
                                         SOLVER& solver,
                                         const FIELDBINDER& fieldBinder,
                                         const KERNEL&      kernelObj,
                                         const bool         incremental = true )
        {
            typedef typename FIELDTUPLEBINDER::Tuple ElementPtrTuple;
            
            // type of stiffness matrix assembly object
            typedef StiffnessMatrix<QUADRATURE,SOLVER,ElementPtrTuple> StiffMat;

            // create a kernel function
            typename StiffMat::Kernel kernel =
                boost::bind( &KERNEL::tangentStiffness, 
                             &kernelObj, _1, _2, _3, _4 );

            // Object of the stiff matrix assembler
            StiffMat stiffness( kernel, quadrature, solver, incremental );

            // Apply to all elements
            base::auxi::applyToAllFieldTuple<FIELDTUPLEBINDER>( fieldBinder, stiffness );

            //typename FIELDBINDER::FieldIterator iter = fieldBinder.elementsBegin();
            //typename FIELDBINDER::FieldIterator end  = fieldBinder.elementsEnd();
            //for ( ; iter != end; ++iter ) {
            //    stiffness( FIELDTUPLEBINDER::makeTuple( *iter ) );
            //}

//            const std::size_t numElements = std::distance( iter, end );
//#pragma omp parallel for
//            for ( std::size_t e = 0; e < numElements; e++ ) {
//                stiffness( FIELDTUPLEBINDER::makeTuple( fieldBinder.elementPtr( e ) ) );
//            }

            
        }


    } // namespace asmb
} // namespace base

//------------------------------------------------------------------------------
/** Computation (via Quadrature) of an element matrix and direct assembly.
 *
 *  Given a matrix valued kernel function \f$k\f$, the element matrix reads
 *  \f[
 *       K = \int_{\Omega_e} k dx \approx \sum_g k(\xi_g) w_g
 *  \f]
 *  where the approximation of the integral has already been carried out.
 *  Once, \f$ K \f$ is obtained it has to be inserted into some global storage
 *  according to the degree of freedom numbering associated with the considered
 *  element. Moreover, some of these degrees of freedom might be constrained
 *  to fixed value. Assume that \f$i\f$ is the row counter, \f$j\f$ the column
 *  counter, \f$ I = doF(i) \f$ is local-to-global numbering map and
 *  \f$ i = c \f$ is constrained degree of freedom with value \f$v_c\f$.
 *  Then we get the assignments for the global system \f$ A x = b \f$
 *  \f[
 *      A[ doF(i), doF(j) ] += K[i,j] \quad i,j \neq c
 *  \f]
 *  and
 *  \f[
 *      b[ doF(i) ]  -= v_c K[i, c]   \quad j=c
 *  \f]
 *  Note that there is a symmetry flag which is activated if rows and columns
 *  represent the same finite element space. In such a case, a few function
 *  calls can be omitted.
 *
 *  \tparam QUAD         Quadrature
 *  \tparam SOLVER       Solver
 *  \tparam FIELDTUPLE   Tuple of field element pointers
 */
template<typename QUAD, typename SOLVER, typename FIELDTUPLE>
class base::asmb::StiffnessMatrix
    : public boost::function<void( const FIELDTUPLE& )>
{
public:
    //! @name Template parameter
    //@{
    typedef QUAD          Quadrature;
    typedef SOLVER        Solver;
    typedef FIELDTUPLE    FieldTuple;
    //@}

    //! @name Access types of the tuple
    //@{
    typedef typename FieldTuple::TestElement   TestElement;
    typedef typename FieldTuple::TrialElement  TrialElement;
    //@}

    //! Kernel function
    typedef boost::function< void( const FieldTuple&, 
                                   const typename Quadrature::VecDim&,
                                   const double,
                                   base::MatrixD& ) >  Kernel;

    //! Constructor with kernel function, quadrature and solver
    StiffnessMatrix( Kernel&              kernel,
                     const Quadrature&    quadrature,
                     Solver&              solver,
                     const bool           incremental )
        : kernel_(          kernel ),
          quadrature_(      quadrature ),
          solver_(          solver ),
          incremental_(     incremental )
    { }

    //--------------------------------------------------------------------------
    void operator()( const FieldTuple& fieldTuple )
    {
        // extract test and trial elements from tuple
        TestElement*  testEp  = fieldTuple.testElementPtr();
        TrialElement* trialEp = fieldTuple.trialElementPtr();

        // if pointers are identical, Galerkin-Bubnov scheme
        const bool isBubnov =
            base::auxi::EqualPointers<TestElement,TrialElement>::apply( testEp,
                                                                       trialEp );
        
        // dof activities
        std::vector<base::dof::DoFStatus> rowDoFStatus, colDoFStatus;

        // dof IDs
        std::vector<std::size_t> rowDoFIDs, colDoFIDs;

        // dof values (for constraints, rowDoFValues is just a placeholder)
        std::vector<base::number> rowDoFValues, colDoFValues;

        // dof constraints
        typedef std::pair<unsigned, std::vector< std::pair<base::number,std::size_t> > >
            WeightedDoFIDs;
        std::vector<WeightedDoFIDs> rowConstraints, colConstraints;

        // Collect dof entities from element
        bool doSomething = 
            base::asmb::collectFromDoFs( testEp, rowDoFStatus,
                                         rowDoFIDs, rowDoFValues,
                                         rowConstraints,
                                         incremental_ );

        // if no row dof is ACTIVE or CONSTRAINED, just return
        if ( not doSomething ) return;

        // In case of identical test and trial spaces, just copy
        if ( isBubnov ) {
            colDoFStatus   = rowDoFStatus;
            colDoFIDs      = rowDoFIDs;
            colDoFValues   = rowDoFValues;
            colConstraints = rowConstraints;
        }
        else // otherwise, collect for trial space
            doSomething =
                base::asmb::collectFromDoFs( trialEp, colDoFStatus,
                                             colDoFIDs, colDoFValues,
                                             colConstraints,
                                             incremental_ );

        // if no col dof is ACTIVE or CONSTRAINED, just return
        if ( not doSomething ) return;


        // Compute the element matrix contribution
        base::MatrixD elemMatrix = base::MatrixD::Zero( rowDoFIDs.size(),
                                                        colDoFIDs.size() );
        quadrature_.apply( kernel_, fieldTuple, elemMatrix );

        // assemble element matrix to global system
        base::asmb::assembleMatrix( elemMatrix,
                                    rowDoFStatus, colDoFStatus,
                                    rowDoFIDs, colDoFIDs,
                                    colDoFValues,
                                    rowConstraints, colConstraints,
                                    solver_, isBubnov );
        return;
    }
    
private:
    Kernel&              kernel_;        //!< Kernel function 
    const Quadrature&    quadrature_;    //!< Quadrature object
    Solver&              solver_;        //!< Solver object

    //! Use difference between prescribed and current value
    const bool incremental_;
};


#endif

