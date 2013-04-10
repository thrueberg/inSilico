//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   constrainBoundary.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_constrainboundary_hpp
#define base_dof_constrainboundary_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/shape.hpp>
#include <base/geometry.hpp>
#include <base/LagrangeShapeFun.hpp>
// base/fe includes
#include <base/fe/LagrangeElement.hpp>
#include <base/fe/Policies.hpp>
// base/mesh includes

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename FEBASIS, typename BITER,
                 typename MESH, typename FIELD, typename DIRIFUN>
        void constrainBoundary( BITER first, BITER last,
                                const MESH&  mesh, FIELD& field, 
                                DIRIFUN diriFun );
        
    }
}

//------------------------------------------------------------------------------
/** Constrain degrees of freedom along the boundary using a given functor.
 *
 *  \tparam FEBASIS Finite element basis
 *  \tparam BITER   Iterator over <elementNo,faceNo> pairs
 *  \tparam MESH    Type of domain mesh
 *  \tparam FIELD   Type of FE field
 *  \tparam DIRIFUN Type of functor for the Dirichlet constraints
 */
template<typename FEBASIS, typename BITER,
         typename MESH, typename FIELD, typename DIRIFUN>
void base::dof::constrainBoundary( BITER first, BITER last,
                                   const MESH&  mesh, FIELD& field, 
                                   DIRIFUN diriFun ) 
{
    // element types
    typedef typename FIELD::Element FieldElement;
    typedef typename MESH::Element  GeomElement;
        
    // Type of face for extraction
    static const base::NFace surface =
        base::ShapeSurface<GeomElement::shape>::value;

    // Finite Element of the field
    typedef typename FEBASIS::FiniteElement FiniteElement;

    // DoF face extraction object
    typedef base::fe::FaceExtraction<FiniteElement,surface> FaceExtraction;

    // get support points of the FE function of the field
    typedef typename FieldElement::FEFun::VecDim VecDim;
    boost::array<VecDim,FieldElement::FEFun::numFun> supportPoints;
    FieldElement::FEFun::supportPoints( supportPoints );

    // Go through all boundary elements and new surface elements
    for ( BITER bIter = first; bIter != last; ++bIter ){ 

        // Get relevant data from iterator
        const std::size_t elemID   = bIter -> first;
        const unsigned  faceNumber = bIter -> second;

        // get mesh and field elements
        const GeomElement*   geomElementPtr  = mesh.elementPtr(  elemID );
        const FieldElement*  fieldElementPtr = field.elementPtr( elemID );

        // Get indices from face extraction
        std::vector<unsigned> faceIndices;
        FaceExtraction::apply( faceNumber, faceIndices );
            
        // Get all dofs of current domain element
        std::vector<typename FIELD::DegreeOfFreedom*> elementDoFs;
        for ( typename FieldElement::DoFPtrConstIter doFIter =
                  fieldElementPtr -> doFsBegin();
              doFIter != fieldElementPtr -> doFsEnd();
              ++doFIter ) {
            elementDoFs.push_back( *doFIter );
        }
            
        // Go through geometry and dof of this face
        for ( unsigned n = 0; n < faceIndices.size(); n ++ ) {
                
            // Local number of DoF
            const unsigned localDofNum = faceIndices[n];

            // Local coordinate corresponding to it
            const VecDim xi = supportPoints[ localDofNum ];

            // Compute geometry
            const typename GeomElement::Node::VecDim x =
                base::Geometry<GeomElement>()( geomElementPtr, xi );
                
            // Get access to dof
            typename FIELD::DegreeOfFreedom* doFPtr = elementDoFs[localDofNum];

            // Apply Dirichlet function to DoF
            diriFun( x, doFPtr );
                
        } // end loop over face's dofs
            
            
    } // end loop over boundary faces

    return;
}


#endif
