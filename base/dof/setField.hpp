//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   setField.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_setField_hpp
#define base_dof_setField_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/shape.hpp>
#include <base/geometry.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename MESH, typename FIELD, typename FIELDFUN>
        void setField( const MESH& mesh, FIELD& field, FIELDFUN fieldFun );
    }
}

//------------------------------------------------------------------------------
/** Set degrees of freedom to a value as given by a function of the coordinate.
 *
 *  \tparam MESH     Type of domain mesh
 *  \tparam FIELD    Type of FE field
 *  \tparam FIELDFUN Type of function of x giving the field values
 */
template<typename MESH, typename FIELD, typename FIELDFUN>
void base::dof::setField( const MESH&  mesh, FIELD& field, FIELDFUN fieldFun ) 
{
    // element types
    typedef typename FIELD::Element FieldElement;
    typedef typename MESH::Element  GeomElement;
        
    // get support points of the FE function of the field
    typedef typename FieldElement::FEFun::VecDim VecDim;
    boost::array<VecDim,FieldElement::FEFun::numFun> supportPoints;
    FieldElement::FEFun::supportPoints( supportPoints );

    // Go through all elements
    typename MESH::ElementPtrConstIter geomEIter  = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter geomLast   = mesh.elementsEnd();
    typename FIELD::ElementPtrIter     fieldEIter = field.elementsBegin();
    for ( ; geomEIter != geomLast; ++geomEIter, ++fieldEIter ){ 

        // get mesh and field elements
        const GeomElement*   geomElementPtr  = *geomEIter;
        const FieldElement*  fieldElementPtr = *fieldEIter;

        // Get all dofs of current domain element
        std::vector<typename FIELD::DegreeOfFreedom*> elementDoFs;
        for ( typename FieldElement::DoFPtrConstIter doFIter =
                  fieldElementPtr -> doFsBegin();
              doFIter != fieldElementPtr -> doFsEnd();
              ++doFIter ) {
            elementDoFs.push_back( *doFIter );
        }
            
        // Go through dofs 
        for ( unsigned n = 0; n < elementDoFs.size(); n++ ) {
                
            // Local coordinate corresponding to it
            const VecDim xi = supportPoints[ n ];

            // Compute geometry
            const typename GeomElement::Node::VecDim x =
                base::Geometry<GeomElement>()( geomElementPtr, xi );
                
            // Get access to dof
            typename FIELD::DegreeOfFreedom* doFPtr = elementDoFs[n];

            // Apply Dirichlet function to DoF
            fieldFun( x, doFPtr );
                
        } // end loop over face's dofs
            
            
    } // end loop over boundary faces

    return;
}


#endif
