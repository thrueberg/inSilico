//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Marching.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_marching_hpp
#define base_cut_marching_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <utility>
#include <map>
#include <bitset>
// boost includes
#include <boost/array.hpp>
// base/cut includes
#include <base/cut/DecomposeHyperCube.hpp>
#include <base/cut/MarchingUtils.hpp>
#include <base/cut/NonLinearCutCellUtils.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        //----------------------------------------------------------------------
        namespace detail_{
            template<unsigned DIM, bool ISSIMPLEX>
            struct Decompose;
            
            template<unsigned DIM>
            struct Decompose<DIM,true>
            {
                static const unsigned numSimplices = 1;
                static unsigned apply( const unsigned s, const unsigned v )
                {
                    return v;
                }
            };

            template<unsigned DIM>
            struct Decompose<DIM,false>
            {
                typedef base::cut::DecomposeHyperCube<DIM>   DHC;
                static const unsigned numSimplices = DHC::numSimplices;
                static unsigned apply( const unsigned s, const unsigned v )
                {
                    return  DHC::apply( s, v );
                }
            };
            
        }

        //----------------------------------------------------------------------
        template<base::Shape SHAPE> struct Marching;

        //----------------------------------------------------------------------
        template<base::Shape SHAPE, unsigned DEGREE>
        struct MarchingProxy;

        // Specialisation for linear elements
        template<base::Shape SHAPE>
        struct MarchingProxy<SHAPE,1>
            : base::cut::Marching<SHAPE>
        { };

        template<base::Shape SHAPE, unsigned DEGREE>
        struct GetInteriorSurface;


    }
}

//------------------------------------------------------------------------------
template<base::Shape SHAPE>
struct base::cut::Marching
{
    typedef base::LagrangeShapeFun<1,SHAPE> LinearLagrange;
    static const unsigned dim         = base::ShapeDim<SHAPE>::value;
    static const unsigned numVertices = LinearLagrange::numFun;

    typedef typename base::Vector<dim,double>::Type     VecDim;
    typedef typename base::cut::USimplex<dim  >::Type   VolSimplex;
    typedef typename base::cut::USimplex<dim-1>::Type   SurfSimplex;


    static void apply( const boost::array<double,numVertices>&   signedDistances,
                       const boost::array<unsigned,numVertices>& vertexIndices,
                       std::vector<VecDim>& nodes,
                       std::map<base::cut::Edge,unsigned>&  uniqueNodes,
                       std::vector<SurfSimplex>& surface,
                       std::vector<VolSimplex>&  volumeIn,
                       std::vector<VolSimplex>&  volumeOut,
                       const bool update = false )
    {
        typedef detail_::Decompose<dim,
                                   SHAPE==base::SimplexShape<dim>::value
                                   > Decompose;
            
        // go through all sub-simplices
        for ( unsigned s = 0; s < Decompose::numSimplices; s++ ) {

            // generate index set of volume simplex and the distances
            VolSimplex                        indexSimplex;
            typename DSimplex<dim>::Type distances;
            
            for ( unsigned v = 0; v < dim+1; v++ ) {
                const unsigned index = Decompose::apply( s, v );
                indexSimplex[v] = vertexIndices[ index ];
                if ( update )
                    distances[v] = signedDistances[v];
                else 
                    distances[v] = signedDistances[ indexSimplex[v] ];
            }

            // apply marching to volume simplex
            detail_::ApplyMarchingSimplex<dim,dim>::apply(
                indexSimplex, distances, nodes,
                uniqueNodes, surface,
                volumeIn, volumeOut );
        }

        return;
    }

};

//------------------------------------------------------------------------------
template<base::Shape SHAPE, unsigned DEGREE>
struct base::cut::MarchingProxy
{
    static const base::Shape shape  = SHAPE;
    static const unsigned    degree = DEGREE;
    static const unsigned       dim = base::ShapeDim<shape>::value;

    //! @name Simplex shapes
    //@{
    static const base::Shape  volSimplex = base::SimplexShape<dim>::value;
    static const base::Shape surfSimplex = base::SimplexShape<dim-1>::value;
    //@}

    //! Coordinate type
    typedef typename base::Vector<dim,double>::Type VecDim;

    //! @name Given element type of the cell
    //@{
    typedef base::fe::LagrangeElement<volSimplex,degree>  VolElement;
    static const unsigned numNodes =
        base::fe::LagrangeElement<SHAPE,degree>::numTotalDoFs;
    static const unsigned numVertices =
        base::NumNFaces<shape,base::VERTEX>::value;
    //@}


    //! @name Higher-order simplex elements
    //@{
    typedef base::LagrangeShapeFun<degree,volSimplex>  LagrangeVolSimplexFun;
    typedef base::LagrangeShapeFun<degree,surfSimplex> LagrangeSurfSimplexFun;
    static const unsigned nVolSimplexNodes  = LagrangeVolSimplexFun::numFun;
    static const unsigned nSurfSimplexNodes = LagrangeSurfSimplexFun::numFun;
    typedef boost::array<unsigned,nVolSimplexNodes>  VolumeSimplex;
    typedef boost::array<unsigned,nSurfSimplexNodes> SurfaceSimplex;
    //@}

    //! @name Linear simplex elements
    //@{
    static const unsigned nVolSimplexVerts  =
        base::NumNFaces<volSimplex,base::VERTEX>::value;
    static const unsigned nSurfSimplexVerts =
        base::NumNFaces<surfSimplex,base::VERTEX>::value;
    typedef boost::array<unsigned,nVolSimplexVerts>  LinearVolumeSimplex;
    typedef boost::array<unsigned,nSurfSimplexVerts> LinearSurfaceSimplex;
    //@}


    static void apply( const boost::array<double,numNodes>& signedDistances,
                       const boost::array<unsigned,numVertices>& vertexIndices,
                       std::vector<VecDim>& nodes,
                       std::map<base::cut::Edge,unsigned>& uniqueNodes,
                       std::vector<SurfaceSimplex>& surface,
                       std::vector<VolumeSimplex>&  volumeIn,
                       std::vector<VolumeSimplex>&  volumeOut,
                       const bool update = false )
    {
        //
        const unsigned maxIter = 10;
        const double   tolerance = 1.e-8;

        //
        const std::size_t numOldNodes = nodes.size() - numVertices;
        

        //----------------------------------------------------------------------
        // do first a linear marching operation
        std::vector<LinearSurfaceSimplex> surfaceLinear;
        std::vector<LinearVolumeSimplex>  volumeInLinear, volumeOutLinear;

        boost::array<double,numVertices> signedDistancesLinear;
        for ( unsigned v = 0; v < numVertices; v++ )
            signedDistancesLinear[v] = signedDistances[v];

        Marching<shape>::apply( signedDistancesLinear,
                                vertexIndices,
                                nodes,
                                uniqueNodes, surfaceLinear,
                                volumeInLinear, volumeOutLinear, update );

        //----------------------------------------------------------------------
        // non-linear update of the intersection points
        typename std::map<base::cut::Edge,unsigned>::iterator iIter = uniqueNodes.begin();
        typename std::map<base::cut::Edge,unsigned>::iterator iEnd  = uniqueNodes.end();
        for ( ; iIter != iEnd; ++iIter ) {
            // Index of the intersecting node
            const unsigned nodeNum = iIter -> second;
            // act only, if this node is new
            if ( nodeNum >= numOldNodes ) {
                // Indices of edge
                const unsigned v1 = (iIter -> first).first();
                const unsigned v2 = (iIter -> first).second();
                // Edge initial point and tangent, current point
                const VecDim xi1 = nodes[v1];
                const VecDim tan = nodes[v2] - xi1;
                VecDim        xi = nodes[nodeNum];
                // perform non-linear intersection
                nodes[nodeNum] =
                    NonLinearIntersection<shape,degree>::apply( xi1, tan, xi, 
                                                                signedDistances,
                                                                tolerance, maxIter );
            }
        }

        //----------------------------------------------------------------------
        // generate higher-order simplices
        
        // make a mesh object from the entire cell-mesh
        std::vector<LinearVolumeSimplex> volume;
        volume.insert( volume.end(), volumeInLinear.begin(),  volumeInLinear.end() );
        volume.insert( volume.end(), volumeOutLinear.begin(), volumeOutLinear.end() );

        typedef base::cut::SimplexMesh<dim,dim> LinearSimplexMesh;
        LinearSimplexMesh dummyMesh( nodes, volume );

        // find the surface of the inside mesh
        std::vector<std::pair<std::size_t,unsigned> > surfacePart, boundaryPart;
        base::cut::FindSurface<dim>::apply( nodes, volumeInLinear,
                                            surfaceLinear,
                                            surfacePart, boundaryPart );
        VERIFY_MSG( surfacePart.size() == surfaceLinear.size(),
                    "Could not find the surface elements" );

        // create a higher-order simplex mesh
        nodes.resize( numOldNodes );
        //nodes.clear();
        base::cut::MakeHigherOrderSimplices<dim,degree>::apply(
            dummyMesh, volumeInLinear.size(),
            nodes, volumeIn, volumeOut );

        // extract the surface
        base::cut::ExtractSurface<dim,degree>::apply(
            volumeIn, surfacePart, surface );

        //----------------------------------------------------------------------
        // update 'secondary' nodes of the surface

        // extract the boundary
        std::vector<SurfaceSimplex> boundaryIn;
        base::cut::ExtractSurface<dim,degree>::apply(
            volumeIn, boundaryPart, boundaryIn );


        // collect normal vectors for all surface nodes
        std::vector<VecDim> normals;
        base::cut::CreateNodeNormals<degree,dim>::apply( nodes,
                                                         surface,
                                                         normals );

        // modify the normals of the nodes along the surface boundary
        base::cut::ModifyNormalsAlongBoundary<degree,dim>::apply( nodes,
                                                                  boundaryIn,
                                                                  normals );

        // go through nodes of the surface
        for ( std::size_t s = 0; s < surface.size(); s++ ) {

            const SurfaceSimplex ss = surface[s];
            
            // go through secondary nodes
            for ( unsigned n = dim; n < ss.size(); n++ ) {
                // ID of node to update
                const unsigned nodeNum = ss[ n ];

                // node location and direction
                const VecDim node = nodes[ nodeNum ];
                const VecDim dir  = normals[ nodeNum ];

                // perform non-linear intersection
                nodes[nodeNum] =
                    NonLinearIntersection<shape,degree>::apply( node, dir, node, 
                                                                signedDistances,
                                                                tolerance, maxIter );
            }
        }
        
        return;
      
    }
};

//------------------------------------------------------------------------------
template<base::Shape SHAPE, unsigned DEGREE>
struct base::cut::GetInteriorSurface
{
    static const unsigned dim = base::ShapeDim<SHAPE>::value;

    typedef typename base::Vector<dim>::Type VecDim;
    
    //! @name Simplex shapes
    //@{
    static const base::Shape  volSimplex = base::SimplexShape<dim>::value;
    static const base::Shape surfSimplex = base::SimplexShape<dim-1>::value;
    //@}

    //! @name Higher-order simplex elements
    //@{
    typedef base::LagrangeShapeFun<DEGREE,volSimplex>  LagrangeVolSimplexFun;
    typedef base::LagrangeShapeFun<DEGREE,surfSimplex> LagrangeSurfSimplexFun;
    static const unsigned nVolSimplexNodes  = LagrangeVolSimplexFun::numFun;
    static const unsigned nSurfSimplexNodes = LagrangeSurfSimplexFun::numFun;
    typedef boost::array<unsigned,nVolSimplexNodes>  VolumeSimplex;
    typedef boost::array<unsigned,nSurfSimplexNodes> SurfaceSimplex;
    //@}

    //! @name Boundary of the mesh
    //@{
    typedef base::mesh::MeshBoundary MeshBoundary;
    typedef typename MeshBoundary::BoundaryElementContainer BEC;
    typedef typename MeshBoundary::BoundConstIter           BECIter;
    //@}

    //! @name Face extraction
    //@{
    typedef base::fe::LagrangeElement<volSimplex,DEGREE> GeomElement;
    static const base::NFace surface =
        base::ShapeSurface<GeomElement::shape>::value;
    typedef base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;
    //@}
    
    static void apply( const std::vector<VecDim>& nodes,
                       std::vector<SurfaceSimplex>& surface,
                       const std::vector<VolumeSimplex>& volumeIn,
                       const std::vector<VolumeSimplex>& volumeOut )
    {
        // make a mesh objects from the entire cell-mesh and the inside part
        std::vector<VolumeSimplex> volume;
        volume.insert( volume.end(),  volumeIn.begin(),  volumeIn.end() );
        volume.insert( volume.end(), volumeOut.begin(), volumeOut.end() );
        typedef base::cut::SimplexMesh<dim,dim,DEGREE> SimplexMesh;
        SimplexMesh meshAll( nodes, volume );
        SimplexMesh meshIn(  nodes, volumeIn );

        // extract boundary faces from volume mesh
        MeshBoundary boundaryAll, boundaryIn;
        boundaryAll.create( meshAll.elementsBegin(), meshAll.elementsEnd() );
        boundaryIn.create(  meshIn.elementsBegin(),  meshIn.elementsEnd() );

        // sort does not help for finding?
        //std::sort( boundaryAll.begin(), boundaryAll.end() );
        //std::sort( boundaryIn.begin(),  boundaryIn.end()  );
        // find the faces in boundaryIn, which are not in boundaryAll
        BEC difference;
        BECIter biIter = boundaryIn.begin();
        BECIter biEnd  = boundaryIn.end();
        for ( ; biIter != biEnd; ++biIter ) {
            BECIter find = std::find( boundaryAll.begin(), boundaryAll.end(), *biIter );
            if ( find == boundaryAll.end() ) difference.push_back( *biIter );
        }

        {
            std::cout << "BoundaryIn:  ";
            BECIter biIter = boundaryIn.begin();
            BECIter biEnd  = boundaryIn.end();
            for ( ; biIter != biEnd; ++biIter )
                std::cout << "(" << biIter->first << "," << biIter->second << "), ";
            std::cout << std::endl;
        }
        {
            std::cout << "BoundaryAll:  ";
            BECIter biIter = boundaryAll.begin();
            BECIter biEnd  = boundaryAll.end();
            for ( ; biIter != biEnd; ++biIter )
                std::cout << "(" << biIter->first << "," << biIter->second << "), ";
            std::cout << std::endl;
        }
        {
            std::cout << "Difference:  ";
            BECIter biIter = difference.begin();
            BECIter biEnd  = difference.end();
            for ( ; biIter != biEnd; ++biIter )
                std::cout << "(" << biIter->first << "," << biIter->second << "), ";
            std::cout << std::endl;
        }
        
        // extract the surface
        surface.clear();
        for ( std::size_t b = 0; b < difference.size(); b++ ) {

            const std::size_t elemNum = difference[b].first;
            const unsigned    faceNum = difference[b].second;

            std::vector<unsigned> faceIndices;
            FaceExtraction::apply( faceNum, faceIndices );


            typename SimplexMesh::Element* elemPtr =
                meshIn.elementPtr( elemNum );

            SurfaceSimplex surfaceSimplex;
            for ( std::size_t n = 0; n < faceIndices.size(); n++ ) {
                surfaceSimplex[n] =
                    static_cast<unsigned>(
                        elemPtr -> nodePtr( faceIndices[n] ) -> getID() );
            }

            surface.push_back( surfaceSimplex );
        }
        
        return;
    }
    
};


#endif

