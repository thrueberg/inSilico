//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ComputeSupport.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_computesupport_hpp
#define base_cut_computesupport_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <algorithm>
// base  includes
#include <base/shape.hpp>
#include <base/linearAlgebra.hpp>
#include <base/types.hpp>
// base/kernel includes
#include <base/kernel/Measure.hpp>
// base/cut includes
#include <base/cut/ParametricMeasure.hpp>
// base/asmb includes
#include <base/asmb/FieldBinder.hpp>


//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<typename FIELDTUPLE, typename QUADRATURE, typename MEASURE>
        class ComputeSupport;

        //----------------------------------------------------------------------
        /** Convenience function for the computation of the support area.
         *  Preforms the computation of the degree of freedom support area for
         *  a given field.
         *  \tparam MESH       Type of mesh
         *  \tparam FIELD      Type of field for which support is of interest
         *  \tparam QUADRATURE Type of numerical integration procedure
         *  \param[in]  mesh       Access to the FE mesh
         *  \param[in]  field      Access to the field of interest
         *  \param[in]  quadrature The numerical integration
         *  \param[out] supports   Storage of the support area for all DoFs
         */
        template<typename MESH, typename FIELD, typename QUADRATURE>
        void supportComputation( const MESH&       mesh, 
                                 const FIELD&      field,
                                 const QUADRATURE& quadrature, 
                                 std::vector<double>& supports )
        {
            // number of dofs and storage size
            const std::size_t numDoFs = std::distance( field.doFsBegin(),
                                                       field.doFsEnd() );
            supports.resize( numDoFs );
            
            // bind the field to the mesh
            typedef base::asmb::FieldBinder<const MESH,const FIELD> FieldBinder;
            FieldBinder fieldBinder( mesh, field );
            typedef typename FieldBinder::template TupleBinder<1>::Type TB;

            // support computation object
            base::cut::ComputeSupport<typename TB::Tuple,QUADRATURE,
                                      //base::kernel::Measure<typename TB::Tuple>
                                      base::cut::ParametricMeasure<typename TB::Tuple>
                                      >
                supportComputer( quadrature, supports );

            typename FieldBinder::FieldIterator fIter = fieldBinder.elementsBegin();
            typename FieldBinder::FieldIterator fEnd  = fieldBinder.elementsEnd();
            for ( ; fIter != fEnd; ++fIter )
                supportComputer( TB::makeTuple( *fIter ) );
            
        return;            
        }

    }
}

//------------------------------------------------------------------------------
/** Computation of the support of degrees of freedom.
 *  The support \f$ supp(\phi_K) \f$ of a shape function \f$ \phi_K \f$
 *  associated with the degree of freedom \f$ K \f$ is defined as
 *  \f[
 *       supp(\phi_K) = \{ x \in \Omega: \phi_K(x) \neq 0 \}
 *  \f]
 *  (More precisely, the support is the \e closure of the above defined set).
 *  This object computes the support area element-wise by integration over
 *  all elements and addition of the area contribution to a global storage.
 *  The most common application is the use in immersed methods in which the
 *  the support size is used to turn degrees of freedom on and off or use it
 *  for preconditioning of the basis.
 *  \tparam FIELDTUPLE  Tuple of a geometry and a field element
 *  \tparam QUADRATURE  Area integration (most likely a base::cut::Quadrature)
 *  \tparam MEASURE     Type of measure to be used for the computation
 */
template<typename FIELDTUPLE, typename QUADRATURE, typename MEASURE>
class base::cut::ComputeSupport
{
public:
    //! Template parameters
    //@{
    typedef FIELDTUPLE FieldTuple;
    typedef QUADRATURE Quadrature;
    typedef MEASURE    Measure;
    //@}

    //! The first field element of the tuple is considered
    typedef typename FieldTuple::template Binder<1>::Type FieldElementPtr;

    //! Deduce the basic type of object
    typedef typename base::TypeReduction<FieldElementPtr>::Type  FieldElement;


    //--------------------------------------------------------------------------
    /** Constructor with quadrature and support storage.
     *  Initializes the storage with zero area
     *  \param[in] quadrature  Quadrature to perform numerical area integration
     *  \param[in] doFSupports Storage of the degree of freedom supports
     */
    ComputeSupport( const Quadrature& quadrature, 
                    std::vector<double>& doFSupports )
        : quadrature_( quadrature ),
          doFSupports_( doFSupports )
    {
        std::fill( doFSupports_.begin(), doFSupports_.end(), 0. );
    }

    /** Overloaded function call operator to perform integration.
     *  First, the area of the element is computed with a numerical integration
     *  rule. Note that the used quadrature holds the information of cut-cells
     *  and only acts on the 'physical' part of the element in case of an
     *  immersed method. Second, the computed area is added to the storage
     *  which holds the support area for all degrees of freedom.
     *  \param[in] fieldTuple Tuple of mesh and field element
     */
    void operator()( const FieldTuple& fieldTuple )
    {
        // compute the 'measure' of the element with a quadrature
        Measure measure;
        double area = 0.;
        quadrature_.apply( measure, fieldTuple, area );

        // get field element pointer from tuple
        FieldElementPtr fElement = fieldTuple.template get<1>();
        
        // add to supports
        typename FieldElement::DoFPtrConstIter doFIter = fElement -> doFsBegin();
        typename FieldElement::DoFPtrConstIter doFEnd  = fElement -> doFsEnd();
        
        // add computed area to storage places of element's degrees of freedom
        for ( unsigned d = 0; doFIter != doFEnd; ++doFIter, d++ )
            doFSupports_[ (*doFIter) -> getID() ] += area;

        return;
    }

private:
    const Quadrature&    quadrature_;  //!< Numerical integration object
    std::vector<double>& doFSupports_; //!< Storage of all DoF support areas
};


#endif
