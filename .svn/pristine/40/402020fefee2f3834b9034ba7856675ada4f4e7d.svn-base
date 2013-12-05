//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   GenerateBoundaryTriangulationFromBox.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_mesh_generateboundarytriangulationfrombox_hpp
#define base_mesh_generateboundarytriangulationfrombox_hpp

//------------------------------------------------------------------------------
// std includes
#include <iterator>
// boost includes
#include <boost/array.hpp>
// base includes
#include <base/shape.hpp>
#include <base/LagrangeShapeFun.hpp>
// base/mesh includes
#include <base/mesh/SurfaceElement.hpp>
#include <base/mesh/Unstructured.hpp>
#include <base/mesh/ParameterSurface.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace mesh{

        //----------------------------------------------------------------------
        //template<typename MESH>
        //class GenerateBoundaryTriangulationFromConnectivity
        //{
        //    // Idea:
        //    // *  First read in a file which identifies a mesh boundary.
        //    //    The file could contain the surface connectivity
        //    //    (e.g., smf) and use mesh::Unstructured for that purpose.
        //    // *  Next, this object extracts the volume mesh's (MESH)
        //    //    element surfaces which match the surface mesh
        //    STATIC_ASSERT_MSG( sizeof( MESH ) == 0,
        //                       "Not implemented! "  );
        //};

        //----------------------------------------------------------------------
        template<typename GRID>
        class GenerateBoundaryTriangulationFromBox;

        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            //! Number of simplex elements on a face of the volume element
            template<base::Shape SHAPE>
            struct NumSimplicesPerElementSurface
            {
                static const unsigned value = 1;
            };

            template<>
            struct NumSimplicesPerElementSurface<base::HEX>
            {
                static const unsigned value = 2;
            };

            //------------------------------------------------------------------
            //! Convert direction and (begin/end) to a surface ID
            template<base::Shape SHAPE>
            unsigned convertToSurfaceID( const unsigned dir,
                                         const bool begin )
            {
                const unsigned dim = base::ShapeDim<SHAPE>::value;

                return ( static_cast<unsigned>( begin ) * dim +
                         dir );
            }

            //------------------------------------------------------------------
            template<base::Shape SHAPE>
            struct ExtractSurfaceSimplex
            {
                static const base::Shape surfShape = base::FaceShape<SHAPE>::value;
                
                static const unsigned numVertices =
                    base::NumNFaces<surfShape,base::VERTEX>::value;
                
                void operator()( const boost::array<unsigned,numVertices> & in,
                                 const unsigned numSimplex, 
                                 boost::array<unsigned,numVertices> & out ) const
                {
                    out = in;
                }
            };

            template<>
            struct ExtractSurfaceSimplex<base::HEX>
            {
                void operator()( const boost::array<unsigned,4> & in,
                                 const unsigned numSimplex,
                                 boost::array<unsigned,3> & out ) const
                {
                    if ( numSimplex == 0 )
                        out = {{ in[0], in[1], in[2] }};
                    else
                        out = {{ in[3], in[2], in[1] }};
                }
            };

        }
        
    }
}

//------------------------------------------------------------------------------
template<typename GRID>
class base::mesh::GenerateBoundaryTriangulationFromBox
{
public:
    //! template parameter: type of grid
    typedef GRID Grid;

    //! @name type of surface mesh deduced from the grid
    //@{
    typedef typename Grid::Element                    VolumeElement;
    typedef base::mesh::SurfaceElement<VolumeElement> SurfaceElement;
    typedef base::mesh::Unstructured<SurfaceElement>  SurfaceMesh;
    typedef typename SurfaceElement::Node             SurfaceNode;
    //@}

    //! @name Multi-Index stuff
    //@{
    typedef typename Grid::MultiIndex     MultiIndex;
    typedef typename Grid::MultiIndexType MultiIndexType;
    //@}

    //! Lookup of parameter space surface
    typedef base::mesh::ParameterSurface<VolumeElement::shape> ParameterSurface;

    //! @name Involved dimensions
    //@{
    static const unsigned volumeDim  = base::ShapeDim<VolumeElement::shape >::value;
    static const unsigned surfaceDim = base::ShapeDim<SurfaceElement::shape>::value;
    //@}

    //! Type of linear surface simplex shape function
    static const base::Shape surfaceSimplex = base::SimplexShape<surfaceDim>::value;
    typedef base::LagrangeShapeFun<1,surfaceSimplex>  LinearSimplexFun;

    static const unsigned numSurfSimplexVertices =
        base::NumNFaces<surfaceSimplex,base::VERTEX>::value;

    //! Type of surface shape function
    typedef typename SurfaceElement::GeomFun          SurfaceShapeFun;

    //! Type of linear volume shape function
    typedef base::LagrangeShapeFun<1,VolumeElement::shape> LinearVolumeFun;


    //--------------------------------------------------------------------------
    /**
     */
    void operator()( const Grid & grid,
                     const unsigned dir,
                     const bool begin,
                     SurfaceMesh & surfaceMesh ) const
    {
        //! Sanity check of direction
        VERIFY_MSG( (dir < Grid::dim), "Input argument for direction wrong" );

        //! Dimensions of the grid
        const MultiIndexType gridSizes = grid.gridSizes();

        //! Multi-index component to compare with
        const int mIndexComp = ( begin == false ? 0 : gridSizes[dir]-1 );

        //! Number of elements in the structured grid
        const std::size_t numElements = MultiIndex::length( gridSizes );

        //! Number of elements on the requested surface
        const std::size_t numElementsOnSurface = numElements / gridSizes[ dir ];

        //! Construct parametric geometry of the element
        const unsigned surfaceID =
            detail_::convertToSurfaceID<VolumeElement::shape>( dir, begin );

        //! List of indices of the parameter space vertices which form the surface
        typename ParameterSurface::Surface
            parameterSurfaceIndices = ParameterSurface::surfaceTable[ surfaceID ];

        //! Interpolation points of the surface shape function
        boost::array< typename SurfaceShapeFun::VecDim,
                      SurfaceShapeFun::numFun> surfaceSupportPoints;
        SurfaceShapeFun::supportPoints( surfaceSupportPoints );

        
        //! Vertices of the parameter volume
        boost::array< typename LinearVolumeFun::VecDim,
                      LinearVolumeFun::numFun> volumeVertices;
        LinearVolumeFun::supportPoints( volumeVertices );

        // --> Reorder for hiearchical ordering
        //     Then, the object ParameterFaces< shape, dim-1> can be used and
        //     ParamaterSurface discarded.

        //! Iterator to surface mesh elements
        typename SurfaceMesh::ElementPtrIter surfElemIter
            = surfaceMesh.elementsBegin();

        //! Linear shape function on the surface simplex
        LinearSimplexFun linearSimplexFun;

        //! Hexahedra will have two triangles per face
        const unsigned numSurfElementsPerElement = 
            detail_::NumSimplicesPerElementSurface<VolumeElement::shape>::value;

        //! Temporary storage of surface elements and nodes
        std::vector<SurfaceElement*> surfaceElements;
        surfaceElements.reserve( numElementsOnSurface * numSurfElementsPerElement );
        std::vector<SurfaceNode*>    surfaceNodes;

        //! node counter
        std::size_t nodeCtr = 0;

        //! Go through all elements
        for ( std::size_t e = 0; e < numElements; e ++ ) {

            //! Construct multi-index from linear counter
            const MultiIndexType eM = MultiIndex::wrap( e, gridSizes );

            //! Check if on requested boundary surface
            const bool onBoundary = ( eM[ dir ] == mIndexComp );

            //! If so, construct the surface element(s)
            if ( onBoundary ) {

                //! Get pointer to volume element
                VolumeElement * vep = grid.elementPtr( eM );

                //! Go through elements on the surface of the volume element
                for ( unsigned se = 0; se < numSurfElementsPerElement; se ++ ) {

                    //! Create new surface element
                    SurfaceElement * surfElem = new SurfaceElement;

                    //! Extract the surface simplex's vertices
                    boost::array<unsigned, numSurfSimplexVertices> surfSimplex;
                    detail_::ExtractSurfaceSimplex<VolumeElement::shape>()( parameterSurfaceIndices,
                                                                            se, surfSimplex );

                    //! Set volume element pointer
                    surfElem -> setVolumeElementPointer( vep );

                    //! Access to surface elements geometry nodes
                    typename SurfaceElement::NodePtrIter nodePtrIter =
                        surfElem -> nodesBegin();

                    //! Go through the parameter points of the surface element
                    typename SurfaceElement::ParamIter paramIter =
                        surfElem -> parametricBegin();
                    typename SurfaceElement::ParamIter paramEnd  =
                        surfElem -> parametricEnd();

                    for ( unsigned p = 0; paramIter != paramEnd;
                          ++paramIter, ++nodePtrIter, p++ ) {

                        //! ..
                        typename base::Vector<surfaceDim>::Type eta
                            = surfaceSupportPoints[ p ];

                        typename LinearSimplexFun::FunArray phi;
                        linearSimplexFun.evaluate( eta, phi );

                        typename base::Vector<volumeDim>::Type xi
                            = base::constantVector<volumeDim>( 0. );

                        for ( unsigned s = 0; s < phi.size(); s ++ ) {
                            xi += phi[s] * volumeVertices[ surfSimplex[s] ];

                        }

                        *paramIter = xi;

                        //! evaluate geometry at the parameter point
                        const typename VolumeElement::Node::VecDim x =
                            base::Geometry<VolumeElement>()( vep, xi );

                        //! convert to vector and pass to node
                        std::vector<double> xV( x.size() );
                        for ( int d = 0; d < x.size(); d ++ )
                            xV[d] = x[d];

                        SurfaceNode * surfNode = new SurfaceNode;
                        surfNode -> setX( xV.begin() );
                        surfNode -> setID( nodeCtr++ );
                        
                        *nodePtrIter = surfNode;

                        surfaceNodes.push_back( surfNode );
                                                
                    } // end loop over surface element's nodes

                    //! Store element pointer
                    surfaceElements.push_back( surfElem );
                    
                } // end loop over simplices per volume element
                
            } // end condition if volume element lies on requested bdry

        }// end loop over all volume elements

        surfaceMesh.allocate( surfaceNodes.size(), surfaceElements.size() );

        std::copy( surfaceNodes.begin(), surfaceNodes.end(),
                   surfaceMesh.nodesBegin() );
        
        std::copy( surfaceElements.begin(), surfaceElements.end(),
                   surfaceMesh.elementsBegin() );
        
        return;
    }
};

#endif
