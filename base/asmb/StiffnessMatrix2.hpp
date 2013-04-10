//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   StiffnessMatrix2.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_asmb_stiffnessmatrix2_hpp
#define base_asmb_stiffnessmatrix2_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
// boost includes
#include <boost/function.hpp>
// base includes
#include <base/linearAlgebra.hpp>
// base/aux includes
#include <base/aux/algorithms.hpp>
// base/asmb includes
#include <base/asmb/StiffnessMatrix.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace asmb{

        template<typename QUAD, typename SOLVER,
                 typename GEOMELEMENT,
                 typename PRIMALELEMENT, typename AUXELEMENT>
        class StiffnessMatrix2;
        
    
        //----------------------------------------------------------------------
        //! Convenience function 
        template<typename QUADRATURE,  typename SOLVER,   typename MESH,
                 typename PRIMALFIELD, typename AUXFIELD, typename KERNEL>
        void stiffnessMatrixComputation2(
            const QUADRATURE&  quadrature,
            SOLVER&            solver,
            const MESH&        mesh,
            const PRIMALFIELD& primalField,
            const AUXFIELD&    auxField,
            const KERNEL&      kernelObj,
            void (KERNEL::*kernelFun)(const typename MESH::Element*,
                                      const typename PRIMALFIELD::Element*,
                                      const typename AUXFIELD::Element*,
                                      const typename QUADRATURE::VecDim&,
                                      const double,
                                      base::MatrixD& ) const, 
            const bool zeroConstraints = false )
        {
            typedef StiffnessMatrix2<QUADRATURE,SOLVER,
                                     typename MESH::Element,
                                     typename PRIMALFIELD::Element,
                                     typename AUXFIELD::Element> StiffMat;
        
            typename StiffMat::Kernel kernel =
                boost::bind( kernelFun, &kernelObj, _1, _2, _3, _4, _5, _6 );

            StiffMat stiffness( kernel, quadrature, solver, zeroConstraints );
            base::aux::forEach3( mesh.elementsBegin(), mesh.elementsEnd(),
                                 primalField.elementsBegin(),
                                 auxField.elementsBegin(), stiffness );
        }

        //----------------------------------------------------------------------
        //! Specialised version for kernel functors
        template<typename QUADRATURE,  typename SOLVER,   typename MESH,
                 typename PRIMALFIELD, typename AUXFIELD, typename KERNEL>
        void stiffnessMatrixComputation2(
            const QUADRATURE&  quadrature,
            SOLVER&            solver,
            const MESH&        mesh,
            const PRIMALFIELD& primalField,
            const AUXFIELD&    auxField,
            const KERNEL&      kernelObj,
            const bool zeroConstraints = false )
        {
            stiffnessMatrixComputation2( quadrature, solver, mesh,
                                         primalField, auxField, kernelObj,
                                         &KERNEL::operator() );
        }

    } // namespace asmb    
} // namespace base


//------------------------------------------------------------------------------
/** Computation and assembly of element stiffness matrices.
 *  Other than in base::StiffnessMatrix, here the matrix is hard-wired to the
 *  (logically) symmetric case, i.e. row and column DoFs stem from the same
 *  space. The extra field provided by AUXELEMENT is passed on to the kernel
 *  function where it will be needed as an additional field. Currently, this
 *  is used for solid::Incompressible, where the main diagonal blocks are
 *  discretised in the symmetric sense but in addition another field has to be
 *  evaluated. In detail, the displacement block contains the pressure variable
 *  and the pressure block contains the volume energy which depends on the
 *  Jacobian which, in turn, depends on the current displacement field.
 *  \tparam QUAD          Quadrature
 *  \tparam SOLVER        Solver
 *  \tparam GEOMELEMENT   Element for the geometry
 *  \tparam PRIMALELEMENT Element for test AND trial space of the discretisation
 *  \tparam AUXELEMENT    Element used for evaluation of an auxiliary field
 */
template<typename QUAD, typename SOLVER,
         typename GEOMELEMENT, typename PRIMALELEMENT, typename AUXELEMENT>
class base::asmb::StiffnessMatrix2
    : public boost::function<void( const GEOMELEMENT*,
                                   const PRIMALELEMENT*,
                                   const AUXELEMENT*) >
{
public:
    //! @name Template parameter
    //@{
    typedef QUAD          Quadrature;
    typedef SOLVER        Solver;
    typedef GEOMELEMENT   GeomElement;
    typedef PRIMALELEMENT PrimalElement;
    typedef AUXELEMENT    AuxElement;
    //@}

    //! Kernel function
    typedef boost::function< void( const GeomElement*,
                                   const PrimalElement*,
                                   const AuxElement*,
                                   const typename Quadrature::VecDim&,
                                   const double,
                                   base::MatrixD& ) >  Kernel;

    //! Constructor with kernel function, quadrature and solver
    StiffnessMatrix2( Kernel&              kernel,
                      const Quadrature&    quadrature,
                      Solver&              solver,
                      const bool zeroConstraints = false )
        : kernel_(          kernel ),
          quadrature_(      quadrature ),
          solver_(          solver ),
          zeroConstraints_( zeroConstraints )
    { }

    //--------------------------------------------------------------------------
    //
    void operator()( const GeomElement*    geomEp,
                     const PrimalElement*  primalEp,
                     const AuxElement*     auxEp )
    {
        // dof activities
        std::vector<bool> rowDoFActivity, colDoFActivity;

        // dof IDs
        std::vector<std::size_t> rowDoFIDs, colDoFIDs;

        // dof values (for constraints)
        std::vector<base::number> rowDoFValues, colDoFValues;

        // Collect dof entities from element
        detail_::collectFromDoFs( primalEp, rowDoFActivity,
                                  rowDoFIDs, rowDoFValues );

        // Hard-wired: test and trial spaces are equal
        colDoFActivity = rowDoFActivity;
        colDoFIDs      = rowDoFIDs;
        colDoFValues   = rowDoFValues;

        // Compute the element matrix contribution
        base::MatrixD elemMatrix = base::MatrixD::Zero( rowDoFIDs.size(),
                                                        colDoFIDs.size() );
        {
            // do the quadrature loop
            typename Quadrature::Iter qIter = quadrature_.begin();
            typename Quadrature::Iter qEnd  = quadrature_.end();
            for ( ; qIter != qEnd; ++qIter ) {

                // Call kernel function for quadrature point
                kernel_( geomEp, primalEp, auxEp,
                         qIter -> second, qIter -> first,
                         elemMatrix );
                
            }
        }

        // assemble element matrix to global system
        detail_::assembleMatrix( elemMatrix,
                                 rowDoFActivity, colDoFActivity,
                                 rowDoFIDs, colDoFIDs,
                                 colDoFValues, solver_, true,
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
