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
#include <base/aux/algorithms.hpp>

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
         *  \tparam QUADRATURE   Type of quadrature
         *  \tparam SOLVER       Type of solver
         *  \tparam BOUNDFIELD   Type of field compound
         *  \tparam KERNEL       Type of object with kernel function implementation
         */
        template<typename QUADRATURE, typename SOLVER, typename BOUNDFIELD,
                 typename KERNEL>
        void stiffnessMatrixComputation( const QUADRATURE& quadrature,
                                         SOLVER& solver,
                                         const BOUNDFIELD& boundField,
                                         const KERNEL&     kernelObj,
                                         const bool zeroConstraints = false )
        {
            // type of stiffness matrix assembly object
            typedef StiffnessMatrix<QUADRATURE,SOLVER,
                                    typename BOUNDFIELD::ElementPtrTuple> StiffMat;

            // create a kernel function
            typename StiffMat::Kernel kernel =
                boost::bind( &KERNEL::tangentStiffness, 
                             &kernelObj, _1, _2, _3, _4 );

            // Object of the stiff matrix assembler
            StiffMat stiffness( kernel, quadrature, solver, zeroConstraints );

            // Apply to all elements
            std::for_each( boundField.elementsBegin(),
                           boundField.elementsEnd(), stiffness );
            
        }


        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            //! Helper to collect relevant data from the dof objects
            template<typename FELEMENT>
            void collectFromDoFs( const FELEMENT*           fElement, 
                                  std::vector<bool>&        activity,
                                  std::vector<std::size_t>& ids,
                                  std::vector<number>&      values ) 
            {
                // Get the element's DoF objects
                std::vector<typename FELEMENT::DegreeOfFreedom*> doFs;
                std::copy( fElement -> doFsBegin(),
                           fElement -> doFsEnd(),
                           std::back_inserter( doFs ) );

                for ( unsigned d = 0; d < doFs.size(); d ++ ) {
                    doFs[d] -> getActivity(    std::back_inserter( activity ) ); 
                    doFs[d] -> getIndices(     std::back_inserter( ids ) );
                    doFs[d] -> getConstraints( std::back_inserter( values ) );
                }
            }

            //------------------------------------------------------------------
            // Helper function for assembling the system matrix
            template<typename SOLVER>
            void assembleMatrix( const base::MatrixD&            elemMatrix,
                                 const std::vector<bool>&        rowDoFActivity,
                                 const std::vector<bool>&        colDoFActivity,
                                 const std::vector<std::size_t>& rowDoFIDs,
                                 const std::vector<std::size_t>& colDoFIDs,
                                 const std::vector<number>&      colDoFValues,
                                 SOLVER& solver,
                                 const bool isBubnov, 
                                 const bool zeroConstraints )
            {
                // number of active dofs in the row space
                const unsigned numActiveRowDoFs =
                    std::count_if( rowDoFActivity.begin(), rowDoFActivity.end(),
                                   boost::bind( std::equal_to<bool>(), _1, true ) );
            
                // number of active dofs in the column space
                const unsigned numActiveColDoFs =
                    std::count_if( colDoFActivity.begin(), colDoFActivity.end(),
                                   boost::bind( std::equal_to<bool>(), _1, true ) );

                // Result container
                MatrixD sysMatrix( numActiveRowDoFs, numActiveColDoFs );
                VectorD sysVector = VectorD::Zero( numActiveRowDoFs );

                // Insert to result container
                unsigned rowCtr = 0;
                unsigned colCtr = 0;

                // Go through all rows
                for ( unsigned r = 0; r < rowDoFIDs.size(); r ++ ) {

                    if ( rowDoFActivity[r] ) {

                        colCtr = 0;
                        // Go through all columns
                        for ( unsigned c = 0; c < colDoFIDs.size(); c ++ ) {

                            if ( colDoFActivity[c] ) { // If active insert to matrix container

                                sysMatrix( rowCtr, colCtr ) = elemMatrix( r, c );
                                colCtr++;
                            }
                            else { // If constrained, pass on to RHS

                                sysVector[ rowCtr ] -= colDoFValues[ c ] * elemMatrix( r, c );

                            }
                        }
                        rowCtr++;
                    }
                }

                // Collect active ID numbers
                std::vector<std::size_t> activeRowDoFIDs( numActiveRowDoFs );
                rowCtr = 0;
                for ( unsigned r = 0; r < rowDoFIDs.size(); r ++ ) {
                    if ( rowDoFActivity[r] ) activeRowDoFIDs[ rowCtr++ ] = rowDoFIDs[r];
                }

                std::vector<std::size_t> activeColDoFIDs( numActiveColDoFs );
                if ( isBubnov ) activeColDoFIDs = activeRowDoFIDs;
                else {
                    colCtr = 0;
                    for ( unsigned c = 0; c < colDoFIDs.size(); c ++ )
                        if ( colDoFActivity[c] ) activeColDoFIDs[ colCtr++ ] = colDoFIDs[c];
                }

                // Pass on to system matrix
                solver.insertToLHS( sysMatrix, activeRowDoFIDs, activeColDoFIDs );
                if ( not zeroConstraints )
                    solver.insertToRHS( sysVector, activeRowDoFIDs );
            
                return;
            }

        } // namespace detail_
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
    typedef typename FieldTuple::TestElementPtr   TestElementPtr;
    typedef typename FieldTuple::TrialElementPtr  TrialElementPtr;
    //@}

    //! Bubnov-Galerkin method
    static const bool isBubnov =
        boost::is_same<TestElementPtr,TrialElementPtr>::value;

    //! Kernel function
    typedef boost::function< void( const FieldTuple&, 
                                   const typename Quadrature::VecDim&,
                                   const double,
                                   base::MatrixD& ) >  Kernel;

    //! Constructor with kernel function, quadrature and solver
    StiffnessMatrix( Kernel&              kernel,
                     const Quadrature&    quadrature,
                     Solver&              solver,
                     const bool zeroConstraints = false )
        : kernel_(          kernel ),
          quadrature_(      quadrature ),
          solver_(          solver ),
          zeroConstraints_( zeroConstraints )
    { }

    //--------------------------------------------------------------------------
    void operator()( const FieldTuple& fieldTuple )
    {
        // extract test and trial elements from tuple
        TestElementPtr  testEp  = fieldTuple.testElementPtr();
        TrialElementPtr trialEp = fieldTuple.trialElementPtr();
        
        // dof activities
        std::vector<bool> rowDoFActivity, colDoFActivity;

        // dof IDs
        std::vector<std::size_t> rowDoFIDs, colDoFIDs;

        // dof values (for constraints)
        std::vector<number> rowDoFValues, colDoFValues;

        // Collect dof entities from element
        detail_::collectFromDoFs( testEp, rowDoFActivity,
                                  rowDoFIDs, rowDoFValues );

        if ( isBubnov ) {
            colDoFActivity = rowDoFActivity;
            colDoFIDs      = rowDoFIDs;
            colDoFValues   = rowDoFValues;
        }
        else
            detail_::collectFromDoFs( trialEp, colDoFActivity,
                                      colDoFIDs, colDoFValues );

        // Compute the element matrix contribution
        base::MatrixD elemMatrix = base::MatrixD::Zero( rowDoFIDs.size(),
                                                        colDoFIDs.size() );
        {
            // do the quadrature loop
            typename Quadrature::Iter qIter = quadrature_.begin();
            typename Quadrature::Iter qEnd  = quadrature_.end();
            for ( ; qIter != qEnd; ++qIter ) {

                // Call kernel function for quadrature point
                kernel_( fieldTuple,
                         qIter -> second, qIter -> first,
                         elemMatrix );
                
            }
        }

        // assemble element matrix to global system
        detail_::assembleMatrix( elemMatrix,
                                 rowDoFActivity, colDoFActivity,
                                 rowDoFIDs, colDoFIDs,
                                 colDoFValues, solver_, isBubnov,
                                 zeroConstraints_ );
        return;
    }
    
private:
    Kernel&              kernel_;        //!< Kernel function 
    const Quadrature&    quadrature_;    //!< Quadrature object
    Solver&              solver_;        //!< Solver object

    //! If set to true the constraint DoF values are passed on as zeros
    const bool zeroConstraints_;

};


#endif

