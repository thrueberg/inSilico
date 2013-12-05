//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   CreateBoundaryMesh.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_createboundarymesh_hpp
#define base_mesh_createboundarymesh_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/shape.hpp>
// base/fe includes
#include <base/fe/LagrangeElement.hpp>
#include <base/fe/Policies.hpp>
// base/mesh includes
#include <base/mesh/SurfaceElement.hpp>
#include <base/mesh/Unstructured.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{
        template<typename ELEMENT>
        class CreateBoundaryMesh;
    }
}

//------------------------------------------------------------------------------
/** Create the surface elements and their connectivity for a given boundary
 *  range.
 */
template<typename DOMELEM>
class base::mesh::CreateBoundaryMesh
{
public:
    //! Template parameter: domain element
    typedef DOMELEM                                   DomainElement;

    //! Type of the surface element
    typedef base::mesh::SurfaceElement<DomainElement> SurfaceElement;

    //! Type of mesh container
    static const bool sharedNodes = true;
    typedef base::mesh::Unstructured<SurfaceElement,sharedNodes>  BoundaryMesh;


    //--------------------------------------------------------------------------
    /**
     */
    template<typename BITER, typename DOMAINMESH>
    static void apply( BITER first, BITER last,
                       const DOMAINMESH& domainMesh,
                       BoundaryMesh& boundaryMesh )
    {
        //! Allocate space for node and element storage
        const std::size_t numSurfElements = std::distance( first, last );
        const std::size_t numNodes = std::distance( domainMesh.nodesBegin(),
                                                    domainMesh.nodesEnd() );
        boundaryMesh.allocate( numNodes, numSurfElements );

        //! Copy nodes
        std::copy( domainMesh.nodesBegin(), domainMesh.nodesEnd(),
                   boundaryMesh.nodesBegin() );

        //! Type of geometry function of the domain element
        typedef typename DomainElement::GeomFun  DomainGeomFun;

        //! Geometry element is a LagrangeElement
        typedef base::fe::LagrangeElement<DomainElement::shape,
                                          DomainGeomFun::degree> GeomElement;

        //! Type of face for extraction
        static const base::NFace surface =
            base::ShapeSurface<DomainElement::shape>::value;

        //! Face extraction object
        typedef base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;

        //! Local coordinate type
        typedef typename DomainGeomFun::VecDim  VecDim;
        
        //! Support points of the geometry shape function
        boost::array<VecDim,DomainGeomFun::numFun> supportPoints;
        DomainGeomFun::supportPoints( supportPoints );

        //! Go through all boundary elements and new surface elements
        typename BoundaryMesh::ElementPtrIter surfElemIter =
            boundaryMesh.elementsBegin();
        for ( BITER bIter = first; bIter != last; ++bIter, ++surfElemIter ){ 

            //! Get relevant data from iterator
            const DomainElement* domainElement = domainMesh.elementPtr( bIter -> first );
            const unsigned       faceNumber    = bIter -> second;

            //! Pass pointer to domain element to surface element
            (*surfElemIter) -> setDomainElementPointer( domainElement );

            //! Get indicdes from face extraction
            std::vector<unsigned> faceIndices;
            FaceExtraction::apply( faceNumber, faceIndices );

            //! Access to surface elements geometry nodes
            typename SurfaceElement::NodePtrIter nodePtrIter =
                (*surfElemIter) -> nodesBegin();

            //! Go through the parameter points of the surface element
            typename SurfaceElement::ParamIter paramIter =
                (*surfElemIter) -> parametricBegin();
            typename SurfaceElement::ParamIter paramEnd  =
                (*surfElemIter) -> parametricEnd();

            for ( unsigned p = 0;
                  paramIter != paramEnd;
                  ++paramIter, ++nodePtrIter, p++ ) {

                //! Get local coordinate 
                const VecDim xi = supportPoints[ faceIndices[p] ];

                //! Assign to parameter of surface element
                *paramIter = xi;

                // pass node pointer from old to new mesh
                *nodePtrIter = domainElement -> nodePtr( faceIndices[p] );
                                                
            } // end loop over surface element's nodes

            
        } // end loop over boundary elements

        return;
    } // end function
};



#endif
