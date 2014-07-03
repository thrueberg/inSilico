#ifndef typesandttributes
#define typesandttributes

#include <base/shape.hpp>
#include <base/Unstructured.hpp>

#include <base/mesh/MeshBoundary.hpp>
#include <base/mesh/generateBoundaryMesh.hpp>
#include <base/cut/generateSurfaceMesh.hpp>

//------------------------------------------------------------------------------
template<unsigned DIM, bool SIMPLEX>
struct ElemShape
{
    static const base::Shape value = base::HyperCubeShape<DIM>::value;
};

template<unsigned DIM>
struct ElemShape<DIM,true>
{
    static const base::Shape value = base::SimplexShape<DIM>::value;
};

//------------------------------------------------------------------------------
template<unsigned DIM, bool SIMPLEX>
struct TypesAndAttributes
{
    // template parameter
    static const unsigned                  dim = DIM;
    static const bool                  simplex = SIMPLEX;
    
    // basic attributes of the computation
    static const unsigned             geomDeg  = 1;
    static const base::Shape             shape = ElemShape<dim,simplex>::value;
    static const base::Shape         surfShape = base::SimplexShape<dim-1>::value;
    static const unsigned    kernelDegEstimate = 5;
    static const unsigned              doFSize = dim;
    static const unsigned                nHist = 1;

    //--------------------------------------------------------------------------
    // Mesh types
    typedef          base::Unstructured<shape,geomDeg>  Mesh;
    typedef typename Mesh::Node::VecDim                 VecDim;
    typedef typename base::mesh::BoundaryMeshBinder<Mesh,true>::Type BoundaryMesh;
    typedef typename base::cut::SurfaceMeshBinder<Mesh>::SurfaceMesh SurfaceMesh;

    //--------------------------------------------------------------------------
    // Cell types
    typedef base::cut::Cell<shape>     Cell;
    typedef base::cut::Cell<surfShape> SurfCell;

    //--------------------------------------------------------------------------
    // Quadrature types
    typedef base::cut::Quadrature<kernelDegEstimate,shape>     CutQuadrature;
    typedef base::SurfaceQuadrature<kernelDegEstimate,shape>   SurfaceQuadrature;
    typedef base::cut::Quadrature<kernelDegEstimate,surfShape> SurfCutQuadrature;
};

#endif
