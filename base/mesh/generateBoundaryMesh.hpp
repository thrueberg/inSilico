//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   generateBoundaryMesh.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_generateboundarymesh_hpp
#define base_mesh_generateboundarymesh_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <iterator>
#include <functional>
// boost includes
#include <boost/array.hpp>
#include <boost/bind.hpp>
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

        //! Create the type of a boundary mesh in function of the domain mesh
        template<typename DOMAINMESH,bool SIMPLEX=false>
        struct BoundaryMeshBinder
        {
            typedef typename DOMAINMESH::Element DomainElement;
            static const base::Shape shape =
                (SIMPLEX ?
                 base::SimplexShape<DomainElement::dim-1>::value :
                 base::FaceShape<DomainElement::shape>::value );
            typedef base::mesh::SurfaceElement<DomainElement,shape> BoundaryElement;
            typedef base::mesh::Unstructured<BoundaryElement>       Type;
        };

        //----------------------------------------------------------------------
        template<typename BITER, typename DOMAINMESH, typename BOUNDARYMESH,
                 typename FILTER>
        void generateBoundaryMesh( BITER first, BITER last,
                                   const DOMAINMESH& domainMesh,
                                   BOUNDARYMESH& boundaryMesh, FILTER filter );
        
        //----------------------------------------------------------------------
        namespace detail_{
            
            struct AlwaysPositive
            {
                template<typename ARG>
                bool operator()( ARG& arg ) const { return true; }
            };
            
        }
        
        template<typename BITER, typename DOMAINMESH, typename BOUNDARYMESH>
        void generateBoundaryMesh( BITER first, BITER last,
                                   const DOMAINMESH& domainMesh,
                                   BOUNDARYMESH& boundaryMesh )
        {
            generateBoundaryMesh( first, last, domainMesh, boundaryMesh,
                                  detail_::AlwaysPositive() );
        }

        namespace detail_{

            //! Dummy
            template<unsigned DEGREE,bool DOTRIANGULATE>
            struct NewIndicesForTriangles
            {
                static unsigned nodeNum( const unsigned n, const unsigned t )
                {
                    return n;
                }
            };

            //------------------------------------------------------------------
            /**  Generate indices for a surface triangulation.
             *   In case of the domain discretisation with hexahedra and the
             *   desire to have triangles for the boundary mesh instead of
             *   quadrilaterals, the nodes need to be reassigned to the two
             *   sub-triangles vertices, edges and faces. This object provides
             *   the map
             *   \f[
             *          (i,t)_p \to  n
             *   \f]
             *   where \f$ i \f$ is the index of the node of the quadrilateral
             *   of degree \f$ p \f$, \f$ t=0\f$ is the first and \f$ t=1\f$
             *   the second sub-triangle, and \f$ n \f$ is the index of the
             *   corresponding node of the triangle. This map is straighforward
             *   for \f$ p=1 \f$ (only vertex nodes) and rather simple for
             *   \f$ p=2 \f$ (additional edge nodes). It is rather a nightmare
             *   for higher degrees due to the numbering of the face-indices.
             *   \tparam DEGREE        Polynomial degree of the elements
             */
            template<unsigned DEGREE>
            struct NewIndicesForTriangles<DEGREE,true>
            {
                //! easier to work with q = p - 1
                static const int q = DEGREE-1;
                static const int numVertices = 3;
                static const int numEdgeNodes = 3 * q;

                //! Give the vertex numbers
                static unsigned vertexNum( const unsigned v, const unsigned t )
                {
                    return (v + 2*t) % 4;
                }

                //! Give the edge node numbers
                static unsigned edgeNodeNum( const unsigned e, const unsigned t )
                {
                    const boost::array<int,6> ES =
                        {{4, 4+q, (q+2)*(q+2)-1, 4+2*q, 4+3*q, 4+4*q }};
                    const boost::array<int,6> EI =
                        {{1, 1, -(q+1), 1, 1, +(q+1) }};

                    // find number of edge and local index on that edge
                    const unsigned numEdge = e / q;
                    const unsigned index   = e - q * numEdge;

                    return ES[3*t + numEdge] + index * EI[3*t + numEdge];
                }

                //! Give the face node number
                static unsigned faceNodeNum( const unsigned f, const unsigned t )
                {
                    static const unsigned N = q * (q-1) / 2;
                    boost::array<unsigned,N> faceNodeNumbers;

                    const int sign = 1 - 2*t;
                    const unsigned fs = (t== 0? 4+4*q : (q+2)*(q+2)-1);

                    unsigned k = 0;
                    for (  int j = 0; j < q-1; j++ )
                        for ( int i = j+1; i < q; i++ )
                            faceNodeNumbers[k++] = fs + sign * (j*q+i);
                    
                    return faceNodeNumbers[f];
                }

                //! Call above functions with shifted indices
                static unsigned nodeNum( const unsigned n, const unsigned t )
                {
                    int i = static_cast<int>(n);
                    if ( i < numVertices ) return vertexNum( i, t );
                    i -= numVertices;
                    if ( i < numEdgeNodes ) return edgeNodeNum( i, t );
                    i -= numEdgeNodes;
                    return faceNodeNum( i, t );
                }
            };

            /** Specialisation for p=2, quadratic elements
             *  <pre>
             *
             *      3--6--2          3--6--2       2
             *      |     |          |    /      / | 
             *      7  8  5   -->    7  8      8   5
             *      |     |          | /      /    |
             *      0--4--1          0       0--4--1
             *
             *  </pre>
             */
            template<>
            struct NewIndicesForTriangles<2,true>
            {
                static unsigned nodeNum( const unsigned n, const unsigned t )
                {
                    // give vertex numbers
                    if ( n < 4 ) return ( n + 2*t ) % 4;
                    // give edge numbers
                    return ( t==0 ? 11 - n : 4 + (n-3)*(n-3) );
                }
            };

            /** Specialisation for p=1, linear elements
             *  <pre>
             *
             *      3---2          3---2       2
             *      |   |          |  /      / | 
             *      |   |   -->    | /      /  |
             *      0---1          0       0---1
             *
             *  </pre>
             */
            template<>
            struct NewIndicesForTriangles<1,true>
            {
                static unsigned nodeNum( const unsigned n, const unsigned t )
                {
                    return ( n + 2*t ) % 4;
                }
            };

            //------------------------------------------------------------------
            /** Do (if needed) a boundary triangulation.
             *  In the only case of spatial dimension equals to three, the
             *  boundary logically two-dimensional and the desire of having
             *  triangles instead of quadrilaterals for the boundary, some
             *  node index resorting needs to be done.
             *  This object actives the NewIndicesForTriangles for that case
             *  and provides some useful attributes.
             *  \tparam DOMSHAPE   Shape of the domain element
             *  \tparam BDRYSHAPE  Shape of the boundary element
             *  \tparam DEGREE     Polynomial degree of the geometry functions
             */
            template<base::Shape DOMSHAPE, base::Shape BDRYSHAPE,
                     unsigned DEGREE>
            struct Tringulate
            {
                //! Inherited shape of the domain (HEX -> QUAD ?)
                static const base::Shape faceShape =
                    base::FaceShape<DOMSHAPE>::value;

                //! If the shapes do not conform, triangulation is needed
                static const bool needToTriangulate =
                    (faceShape != BDRYSHAPE);

                //! Dimension of the surface
                static const unsigned surfaceDim =
                    base::ShapeDim<BDRYSHAPE>::value;

                //! This can only happen in the case of a two-dimensional surface
                STATIC_ASSERT_MSG( (surfaceDim == 2) or not needToTriangulate, 
                                   "Triangulation only possible for 2D boundaries" );

                //! If trianguated, 2 triangles per quadrilateral
                static const unsigned numSubElements =
                    (needToTriangulate ? 2 : 1);

                //! The geometry function used for the inherited shape
                typedef base::LagrangeShapeFun<DEGREE,faceShape> LSF;
                //! Number of nodes per surface segment (macro-element)
                static const unsigned numNodes = LSF::numFun;

                //! Object which does the indes mapping
                typedef NewIndicesForTriangles<DEGREE,
                                               needToTriangulate> NIFT;

                //! Act
                static unsigned apply( const unsigned numSub,
                                       const unsigned index )
                {
                    return NIFT::nodeNum( index, numSub );
                }
            };

        }
    }
}

//------------------------------------------------------------------------------
/** Create the surface elements and their connectivity for a given
 *  boundary range.
 *  Given a range of (element-number,face-number) pairs, this object
 *  applies a given geometry filter, creates a surface element from
 *  each such pair and passes a pointer to the adjacent domain element
 *  to it. For each surface element, the geometry of the domain
 *  element is evaluated at the support points of a corresponding
 *  Lagrange element. The resulting coordinates are used to create
 *  nodes which are stored in the mesh container.
 *  \tparam BITER        Type of iterator over the pairs
 *  \tparam DOMAINMESH   Type of domain mesh
 *  \tparam BOUNDARYMESH Type of boundary mesh
 *  \tparam FILTER       Type of geometry filter
 *  \param[in]  first, last  Begin and end iterators
 *  \param[in]  domainMesh   Mesh of the domain
 *  \param[out] boundaryMesh Mesh of the boundary
 *  \param[in]  filter       A geometry filter (bool(x))
 */
template<typename BITER, typename DOMAINMESH, typename BOUNDARYMESH,
         typename FILTER>
void base::mesh::generateBoundaryMesh( BITER first, BITER last,
                                       const DOMAINMESH& domainMesh,
                                       BOUNDARYMESH& boundaryMesh,
                                       FILTER filter)
{
    // Type of element in the domain
    typedef typename DOMAINMESH::Element DomainElement;

    // Type of element on the boundary
    typedef typename BOUNDARYMESH::Element BoundaryElement;

    // Type of geometry function of the domain element
    typedef typename DomainElement::GeomFun  DomainGeomFun;

    // Geometry element is a LagrangeElement
    typedef base::fe::LagrangeElement<DomainElement::shape,
                                      DomainGeomFun::degree> GeomElement;

    // Tringulation object (for TRI instead of QUAD boundary elements)
    typedef
        base::mesh::detail_::Tringulate<DomainElement::shape,
                                        BoundaryElement::shape,
                                        DomainGeomFun::degree> Triangulate;

    // Use its geometry function
    typedef typename GeomElement::ShapeFun GeomFun;

    // Type of face for extraction
    static const base::NFace surface =
        base::ShapeSurface<GeomElement::shape>::value;

    // Face extraction object
    typedef base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;

    // Local coordinate type
    typedef typename GeomFun::VecDim  VecDim;
        
    // Support points of the geometry shape function
    boost::array<VecDim,GeomFun::numFun> supportPoints;
    GeomFun::supportPoints( supportPoints );

    //--------------------------------------------------------------------------
    // Step 1) collect all coordinates and filter them
    std::vector< std::vector<VecDim> > parameterCoordinates, physicalCoordinates;
    std::vector<bool>   filterValues;
    
    for ( BITER bIter = first; bIter != last; ++bIter ){ 

        // Get relevant data from iterator
        DomainElement* domainElement = domainMesh.elementPtr( bIter -> first );
        const unsigned faceNumber    = bIter -> second;

        // Get indicdes from face extraction
        std::vector<unsigned> faceIndices;
        FaceExtraction::apply( faceNumber, faceIndices );

        bool passesFilter = true;

        std::vector<VecDim> xiTmp, xTmp;
        
        // go through all nodes of the face
        for ( unsigned f = 0; f < faceIndices.size(); f++ ) {

            // parameter coordinate inherited from the domain element
            const VecDim xi = supportPoints[ faceIndices[f] ];
            xiTmp.push_back( xi );
            // physical coordinate via geometry evaluation 
            const VecDim x  =
                base::Geometry<DomainElement>()( domainElement, xi );
            xTmp.push_back( x );
            // check filter
            if ( not filter(x) ) passesFilter = false; 
        }

        physicalCoordinates.push_back(  xTmp  );
        parameterCoordinates.push_back( xiTmp );

        // store if this element passes the filter        
        filterValues.push_back( passesFilter );
        
    }

    //--------------------------------------------------------------------------
    // Step 2) Count and Allocate

    // number of valid elements
    const std::size_t numValidElements =
        std::count_if( filterValues.begin(), filterValues.end(),
                       boost::bind( std::equal_to<bool>(), _1, true ) );

    // Allocate space for node and element storage
    const unsigned    numSubElements   = Triangulate::numSubElements;
    const std::size_t numSurfElements  = numValidElements * numSubElements;
    const std::size_t numNodes         = numValidElements * Triangulate::numNodes;
    boundaryMesh.allocate( numNodes, numSurfElements );

    //--------------------------------------------------------------------------
    // Step 3) Construct boundary mesh from valid faces

    // ID of the nodes of the surface mesh
    std::size_t nodeId = 0;
    std::size_t eCtr   = 0;

    // Go through all boundary elements and new boundary elements
    typename BOUNDARYMESH::ElementPtrIter bdryElemIter =
        boundaryMesh.elementsBegin();
    for ( BITER bIter = first; bIter != last; ++bIter, eCtr++ ) {

        // if element has not passed the filter, just continue
        if ( not filterValues[eCtr] ) continue;

        // pointer to domain element
        DomainElement* domainElement = domainMesh.elementPtr( bIter -> first );

        // extract the parameter and physical coordinates
        const std::vector<VecDim> xiTmp = parameterCoordinates[ eCtr ];
        const std::vector<VecDim>  xTmp = physicalCoordinates[  eCtr ];

        // pass data to node pointers
        std::vector<typename DOMAINMESH::Node*> nodePointers;
        for ( unsigned p = 0; p < xiTmp.size(); p++ ) {
            
            // get node pointer and equip with data
            typename DOMAINMESH::Node* np = boundaryMesh.nodePtr( nodeId );
            np -> setID( nodeId );
            const VecDim x = xTmp[p];
            np -> setX( &(x[0]) );

            // store for later use
            nodePointers.push_back( np );

            // increment node ID
            nodeId++;
        }
        
        // go through the sub-elements of the boundary element
        for ( unsigned s = 0; s < numSubElements; s++ ) {
            
            // Pass pointer to domain element to boundary element
            (*bdryElemIter) -> setDomainElementPointer( domainElement );

            // Access to boundary elements geometry nodes
            typename BoundaryElement::NodePtrIter nodePtrIter =
                (*bdryElemIter) -> nodesBegin();

            // Go through the parameter points of the boundary element
            typename BoundaryElement::ParamIter paramIter =
                (*bdryElemIter) -> parametricBegin();
            typename BoundaryElement::ParamIter paramEnd  =
                (*bdryElemIter) -> parametricEnd();

            for ( unsigned p = 0;
                  paramIter != paramEnd;
                  ++paramIter, ++nodePtrIter, p++ ) {

                // Re-assign node index (needed for TRI instead of QUAD)
                const unsigned n = Triangulate::apply( s, p );

                // Assign to pointers
                *paramIter   = xiTmp[n];
                *nodePtrIter = nodePointers[n];
                                                
            } // end loop over boundary element's nodes

            // increment iterator over boundary elements
            ++bdryElemIter;
            
        } // end loop over boundary sub-elements
        
    } // end loop over boundary elements

    return;
} // end function

#endif
