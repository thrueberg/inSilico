//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   NonLinearCutCellUtils.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_cut_nonlinearcutcells_hpp
#define base_cut_nonlinearcutcells_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
#include <utility>
// boost includes
#include <boost/array.hpp>
// base  includes
#include <base/LagrangeShapeFun.hpp>
#include <base/fe/LagrangeElement.hpp>
#include <base/cut/MarchingUtils.hpp>
#include <base/cut/SimplexMesh.hpp>
#include <base/Unstructured.hpp>
#include <base/Field.hpp>
#include <base/fe/Basis.hpp>
#include <base/dof/generate.hpp>
#include <base/dof/location.hpp>
#include <base/mesh/MeshBoundary.hpp>


//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        //----------------------------------------------------------------------
        /** Compute the intersection point based on higher-order distances.
         *  Given the signed distances by means of
         *  \f[
         *       d(\xi) = \sum_K \phi_K(\xi) d_K
         *  \f]
         *  and a line
         *  \f[
         *        \xi(\eta) = \xi_0 + \eta  t
         *  \f]
         *  this function solves the nonlinear equation
         *  \f[
         *        d( \xi( \eta ) ) = 0
         *  \f]
         *  with the Newton method
         *  \f[
         *    \frac{\partial d}{\partial \xi}
         *    \frac{\partial \xi}{\partial \eta}|_{\eta=\eta^k}
         *    (\eta^{k+1} -\eta^k)=  -d( \xi( \eta^k) )
         *  \f]
         *
         *  \tparam SHAPE  Type of shape in which the solution is sought
         *  \tparam DEGREE Polynomial degree of approximation
         */
        template<base::Shape SHAPE, unsigned DEGREE>
        struct NonLinearIntersection
        {
            //! Coordiante type
            typedef typename base::Vector<base::ShapeDim<SHAPE>::value>::Type VecDim;

            /** Perform the Newton method
             *  \param[in] xi1       Starting point \f$ \xi_0 \f$ of line
             *  \param[in] tan       Tangent \f$ t \f$ of the line
             *  \param[in] xiInitial Initial value of the coordinate \f$ \xi \f$
             *  \param[in] signedDistances Array of nodal values \f$ d_K \f$
             *  \param[in] tolerance Stop tolerance of Newton method
             *  \param[in] maxIter   Maximal number of Newton iterations
             *  \return              Latest iterate of \f$ \xi \f$
             */
            template<typename DISTARRAY>
            static VecDim apply( const VecDim xi1,
                                 const VecDim tan,
                                 const VecDim xiInitial,
                                 const DISTARRAY& signedDistances,
                                 const double tolerance,
                                 const unsigned maxIter )
            {
                // initialise the solution iterate
                VecDim xi = xiInitial;
                    
                // Newton method to update the intersection point
                unsigned iter = 0;
                for ( ; iter <= maxIter; iter ++ ) {

                    // residual of the problem
                    const double residual =
                        detail_::InterpolatedDistance<SHAPE,DEGREE>::evaluate(
                            xi, signedDistances );

                    // first convergence test
                    if ( std::abs( residual ) < tolerance ) break;

                    // Consistent tangent of the Newton method
                    const double lhs =
                        base::dotProduct(
                            detail_::InterpolatedDistance<SHAPE,DEGREE>::gradient(
                                xi, signedDistances ),
                            tan );

                    // Solution increment
                    const double dEta = -residual / lhs;

                    // Update
                    xi += dEta * tan;

                    // Second convergence criterion
                    if ( std::abs( dEta ) < tolerance ) break;
                }

                if ( iter == maxIter ) std::cout <<  "Arsch" << std::endl;

                return xi;
            }

        };


        //----------------------------------------------------------------------
        /** Construct a higher-order simplex structure.
         *  Given a linear mesh of simplex elements, construct a mesh of higher
         *  order simplex elements. This is achieved by a temporary base::Field
         *  whose degrees of freedom represent the nodes of the new mesh.
         *  \tparam DIM    Spatial dimension of the problem
         *  \tparam DEGREE Polynomial degree of the higher-order elements
         */
        template<unsigned DIM, unsigned DEGREE>
        struct MakeHigherOrderSimplices
        {
            //! Shape of a Dim-simplex
            static const base::Shape simplex = base::SimplexShape<DIM>::value;

            //! Type of linear simplex mesh
            typedef          base::cut::SimplexMesh<DIM>       LinearSimplexMesh;
            //! Coordinate type
            typedef typename LinearSimplexMesh::VecDim           VecDim;

            //! Basis and definition of dummy FE Field
            typedef base::fe::Basis<simplex,DEGREE>              FEBasis;
            typedef base::Field<FEBasis,1>                       DummyField;

            //! Shape function gives rise to the size of a higher-order simplex
            typedef base::LagrangeShapeFun<DEGREE,simplex>       SFun;
            typedef boost::array<unsigned,SFun::numFun>          VolumeSimplex;

            static void apply( const LinearSimplexMesh& linearSimplexMesh,
                               const std::size_t inOut,
                               std::vector<VecDim>& nodes,
                               std::vector<VolumeSimplex>& volumeIn,
                               std::vector<VolumeSimplex>& volumeOut )
            {
                // Dummy Field with the desired order
                DummyField dummyField;
                base::dof::generate<FEBasis>( linearSimplexMesh, dummyField );
                std::vector<std::pair<std::size_t,VecDim> > location;
                base::dof::associateLocation( dummyField, location );

                std::vector<VecDim> tmpNodes;

                // Location of the DoFs gives the nodes
                for ( std::size_t n = 0; n < location.size(); n++ ) {
                    const VecDim xi =
                        base::Geometry<typename LinearSimplexMesh::Element>()(
                            linearSimplexMesh.elementPtr( location[n].first ),
                            location[n].second );
                    nodes.push_back( xi );
                }

                // The DoF connectivity gives the higher-order mesh connectivity
                typename DummyField::ElementPtrIter fIter = dummyField.elementsBegin();
                typename DummyField::ElementPtrIter fEnd  = dummyField.elementsEnd();
                for ( std::size_t e = 0;  fIter != fEnd; e++, ++fIter ) {
                    VolumeSimplex vs;
                    typename DummyField::Element::DoFPtrConstIter dIter = (*fIter) -> doFsBegin();
                    for ( std::size_t v = 0; v < vs.size(); v++, ++dIter ) {
                        const unsigned vIndex = static_cast<unsigned>( (*dIter) -> getID() );
                        vs[v] = vIndex;
                    }

                    // choose between inside and outside volume
                    if ( e < inOut )
                        volumeIn.push_back( vs );
                    else
                        volumeOut.push_back( vs );
                }

                return;
            }

        };

        //----------------------------------------------------------------------
        /** Find the surface portion of the boundary of the mesh.
         *  Generate a mesh from the inside volume portion and create its
         *  boundary. The faces of that boundary which conincide with the
         *  surface and the ones which do not are stored in seperated
         *  containers.
         *  \tparam DIM  Spatial dimension of the problem
         */
        template<unsigned DIM>
        struct FindSurface
        {
            typedef typename base::Vector<DIM>::Type VecDim;
            typedef boost::array<unsigned,DIM+1>     LinearVolumeSimplex;
            typedef boost::array<unsigned,DIM>       LinearSurfaceSimplex;

            typedef base::cut::SimplexMesh<DIM>      LinearSimplexMesh;

            typedef base::fe::LagrangeElement<LinearSimplexMesh::Element::shape,
                                              1> GeomElement;
            static const base::NFace surface =
                base::ShapeSurface<GeomElement::shape>::value;

            typedef
            base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;

            static void apply( const std::vector<VecDim>& nodes,
                               const std::vector<LinearVolumeSimplex>& elementsIn,
                               const std::vector<LinearSurfaceSimplex>& surface,
                               std::vector<std::pair<std::size_t,unsigned> >&  surfacePart,
                               std::vector<std::pair<std::size_t,unsigned> >& boundaryPart )
            {
                // boundary of the mesh inside the domain
                LinearSimplexMesh dummyMeshIn( nodes, elementsIn );
                typedef base::mesh::MeshBoundary MeshBoundary;
                MeshBoundary meshBoundary;
                meshBoundary.create( dummyMeshIn.elementsBegin(),
                                     dummyMeshIn.elementsEnd() );

                typename MeshBoundary::BoundConstIter bIter = meshBoundary.begin();
                typename MeshBoundary::BoundConstIter bEnd  = meshBoundary.end();
                for ( ; bIter != bEnd; ++bIter ) {
                    const std::size_t elemNum = bIter -> first;
                    const unsigned    faceNum = bIter -> second;

                    std::vector<unsigned> faceIndices;
                    FaceExtraction::apply( faceNum, faceIndices );

                    typename LinearSimplexMesh::Element* volElemPtr =
                        dummyMeshIn.elementPtr( elemNum );

                    LinearSurfaceSimplex boundarySimplex;
                    for ( std::size_t n = 0; n < faceIndices.size(); n++ )
                        boundarySimplex[n] =
                            static_cast<unsigned>(
                                volElemPtr -> nodePtr( faceIndices[n] ) -> getID() );

                    // find in given surface
                    for ( std::size_t s = 0; s < surface.size(); s++ ) {

                        if ( surface[s] == boundarySimplex ) 
                            surfacePart.push_back( std::make_pair( elemNum, faceNum ) );
                        else
                            boundaryPart.push_back( std::make_pair( elemNum, faceNum ) );
                    }
                }

                return;
            }
            

        };

        //----------------------------------------------------------------------
        /** Extract the higher-order simplices that form the surface.
         */
        template<unsigned DIM, unsigned DEGREE>
        struct ExtractSurface
        {
            static const base::Shape  volSimplex = base::SimplexShape<DIM>::value;
            static const base::Shape surfSimplex = base::SimplexShape<DIM-1>::value;
            
            typedef base::LagrangeShapeFun<DEGREE,volSimplex>       VolSFun;
            typedef boost::array<unsigned,VolSFun::numFun>          VolumeSimplex;

            typedef base::LagrangeShapeFun<DEGREE,surfSimplex>      SurfSFun;
            typedef boost::array<unsigned,SurfSFun::numFun>         SurfaceSimplex;
            
            typedef base::fe::LagrangeElement<volSimplex,DEGREE> GeomElement;
            static const base::NFace surface = base::ShapeSurface<volSimplex>::value;

            typedef
            base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;

            static void apply( const std::vector<VolumeSimplex>& volumeIn,
                               const std::vector<std::pair<std::size_t,unsigned> >&
                               surfacePart,
                               std::vector<SurfaceSimplex>& surface )
            {
                typename std::vector<std::pair<std::size_t,unsigned> >::const_iterator
                    sIter = surfacePart.begin();
                typename std::vector<std::pair<std::size_t,unsigned> >::const_iterator
                    sEnd  = surfacePart.end();
                for ( ; sIter != sEnd; ++sIter ) {

                    const std::size_t elemNum = sIter -> first;
                    const unsigned    faceNum = sIter -> second;

                    std::vector<unsigned> faceIndices;
                    FaceExtraction::apply( faceNum, faceIndices );

                    const VolumeSimplex vs = volumeIn[ elemNum ];

                    SurfaceSimplex surfaceSimplex;
                    for ( std::size_t n = 0; n < faceIndices.size(); n++ )
                        surfaceSimplex[n] = vs[ faceIndices[n] ];
                    
                    surface.push_back( surfaceSimplex );
                }

                return;
            }
            
        };

        //----------------------------------------------------------------------
        //! Compute the (un-normalised) outward normal vector to a simplex
        template<unsigned DEGREE, unsigned DIM>
        struct SimplexNormal;

        //! \cond SKIPDOX
        template<unsigned DEGREE>
        struct SimplexNormal<DEGREE,1>
        {
            typedef base::Vector<1>::Type Vec1;
            typedef base::LagrangeShapeFun<DEGREE,base::POINT> SurfSFun;
            typedef boost::array<unsigned,SurfSFun::numFun>   SurfaceSimplex;

            static Vec1 apply( const std::vector<Vec1>& nodes,
                               const SurfaceSimplex& ss )
            {
                const Vec1 tan = nodes[ ss[0] ];
                return tan;
            }
        };

        template<unsigned DEGREE>
        struct SimplexNormal<DEGREE,2>
        {
            typedef base::Vector<2>::Type Vec2;
            typedef base::LagrangeShapeFun<DEGREE,base::LINE> SurfSFun;
            typedef boost::array<unsigned,SurfSFun::numFun>   SurfaceSimplex;

            static Vec2 apply( const std::vector<Vec2>& nodes,
                               const SurfaceSimplex& ss )
            {
                const Vec2 tan = nodes[ ss[1] ] - nodes[ ss[0] ];
                Vec2 normal;
                normal[0] = -tan[1];
                normal[1] =  tan[0];
                return normal;
            }
        };

        template<unsigned DEGREE>
        struct SimplexNormal<DEGREE,3>
        {
            typedef base::Vector<3>::Type Vec3;
            typedef base::LagrangeShapeFun<DEGREE,base::TRI>  SurfSFun;
            typedef boost::array<unsigned,SurfSFun::numFun>   SurfaceSimplex;

            static Vec3 apply( const std::vector<Vec3>& nodes,
                               const SurfaceSimplex& ss )
            {
                const Vec3 tan1 = nodes[ ss[1] ] - nodes[ ss[0] ];
                const Vec3 tan2 = nodes[ ss[2] ] - nodes[ ss[0] ];
                return base::crossProduct( tan1, tan2 );
            }
        };
        //! \endcond
        
        //----------------------------------------------------------------------
        /** For all nodes on the immersed surface, compute a normal vector.
         *  The non-linear update of the 'secondary' surface nodes shall take
         *  place in the direction normal to the surface. Therefore, the
         *  piece-wise constant simplex normals are registered at the nodes
         *  and divided by the sum of the area. This gives area-averaged vertex
         *  normals.
         *  \tparam DEGREE  Polynomial degree of the surface simplices
         *  \tparam DIM     Spatial dimension of the problem
         */
        template<unsigned DEGREE, unsigned DIM>
        struct CreateNodeNormals
        {

            typedef typename base::Vector<DIM>::Type VecDim;
            static const base::Shape surfSimplex = base::SimplexShape<DIM-1>::value;
            
            typedef base::LagrangeShapeFun<DEGREE,surfSimplex>      SurfSFun;
            typedef boost::array<unsigned,SurfSFun::numFun>         SurfaceSimplex;

            static void apply( const std::vector<VecDim>& nodes,
                               const std::vector<SurfaceSimplex>& surface,
                               std::vector<VecDim>& normals )
            {
                normals.resize( nodes.size(), base::constantVector<DIM>( 0. ) );

                std::vector<double> areas( nodes.size(), 0. );

                for ( std::size_t s = 0; s < surface.size(); s++ ) {

                    const VecDim n = 
                        base::cut::SimplexNormal<DEGREE,DIM>::apply( nodes,
                                                                     surface[s] );
                    const double area = n.norm();
                    
                    for ( unsigned v = 0; v < SurfSFun::numFun; v++ ) {
                        normals[ surface[s][v] ] += n;
                        areas[   surface[s][v] ] += area;
                    }
                }

                for ( std::size_t n = 0; n < normals.size(); n++ ) {
                    const double area = areas[n];
                    if ( area > 0. ) normals[n] /= area;
                }
                
            }

        };

        //----------------------------------------------------------------------
        //! In 1 and 2D, there are no secondary nodes on the cell boundary
        template<unsigned DEGREE,unsigned DIM>
        struct ModifyNormalsAlongBoundary
        {
            template<typename DUMMY1, typename DUMMY2>
            static void apply( const DUMMY1& dummy1,
                               const DUMMY2& dummy2,
                               DUMMY1& dumm3 )
            {
                return; 
            }
        };
        
        /** Remove outward cell-normal component for boundary-of-surface nodes.
         *  In 3D, the vertex normal of secondary node on the boundary of the
         *  cell might point out of the face on which the load nodes. This could
         *  lead to the undesirable effect that the node lies outside of the
         *  cell. Therefore, the normal vectors of the secondary nodes along the
         *  line which conincides with the cell boundary are modified, such that
         *  they lie in the plane of the boundary of the cell.
         *  The cell boundary surface is created outside of this object and
         *  yields the tangent plane for this modification.
         *  \tparam DEGREE  Polynomial degree of the geometry approximation
         *  \tparam DIM     Dimension of the problem
         */
        template<unsigned DEGREE>
        struct ModifyNormalsAlongBoundary<DEGREE,3>
        {
            typedef base::Vector<3>::Type Vec3;
            static const base::Shape surfSimplex = base::TRI;
            
            typedef base::LagrangeShapeFun<DEGREE,surfSimplex>      SurfSFun;
            typedef boost::array<unsigned,SurfSFun::numFun>         SurfaceSimplex;

            typedef base::fe::LagrangeElement<surfSimplex,DEGREE> GeomElement;
            static const base::NFace surface = base::EDGE;

            typedef
            base::fe::FaceExtraction<GeomElement,surface> FaceExtraction;
            
            static void apply( const std::vector<Vec3>& nodes, 
                               const std::vector<SurfaceSimplex>& boundaryIn,
                               std::vector<Vec3>& normals )
            {
                // create mesh from surface elements
                typedef base::cut::SimplexMesh<3,2> LinearSimplexMesh;
                LinearSimplexMesh boundaryMesh( nodes, boundaryIn );

                // create the boundary of that surface mesh
                typedef base::mesh::MeshBoundary MeshBoundary;
                MeshBoundary meshBoundary;
                meshBoundary.create( boundaryMesh.elementsBegin(),
                                     boundaryMesh.elementsEnd() );

                // extract the faces
                typename MeshBoundary::BoundConstIter bIter = meshBoundary.begin();
                typename MeshBoundary::BoundConstIter bEnd  = meshBoundary.end();
                for ( ; bIter != bEnd; ++bIter ) {

                    const std::size_t elemNum = bIter -> first;
                    const unsigned    faceNum = bIter -> second;

                    // normal to the boundary element
                    const Vec3 normal2 = 
                        base::cut::SimplexNormal<DEGREE,3>::apply( nodes,
                                                                   boundaryIn[elemNum] );

                    std::vector<unsigned> faceIndices;
                    FaceExtraction::apply( faceNum, faceIndices );

                    typename LinearSimplexMesh::Element* surfElemPtr =
                        boundaryMesh.elementPtr( elemNum );

                    for ( unsigned d = 0; d < DEGREE-1; d++ ) {
                        // index of the secondary node on boundary line
                        const std::size_t nodeID =
                            surfElemPtr -> nodePtr( faceIndices[d+1] ) -> getID();
                        
                        // its normal vector
                        const Vec3 normal1 = normals[ nodeID ];

                        // outward pointing part
                        const Vec3 outward =
                            base::dotProduct( normal1, normal2 ) * normal2;

                        // remove this part
                        normals[ nodeID ] -= outward;
                    }

                }

                return;
            }

        };

    }
}

#endif
