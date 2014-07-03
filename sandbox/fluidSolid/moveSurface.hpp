#ifndef movesurface_h
#define movesurface_h

#include <vector>

#include <base/linearAlgebra.hpp>
#include <base/cut/generateSurfaceMesh.hpp>
#include <base/post/evaluateField.hpp>


//------------------------------------------------------------------------------
/** Move the immersed surface by evaluation of the velocity increment.
 *  The new surface coordinate is computed as
 *  \f[
 *      X_{n+1} = X_n + v_n \Delta t + (\alpha-1) (X_n - c_n) 
 *  \f]
 *  Here, \f$ v_n \Delta t \f$ is the current displacement increment (velocity
 *  times time step size) and there is a correction term which ensures that the
 *  transition between implicit and explicit surface representations prior to
 *  the call of this function is without change in confined volume.
 *  
 */
template<typename MESH, typename FIELD>
void moveSurface( const MESH& mesh,
                  const FIELD& velocity,
                  typename base::cut::SurfaceMeshBinder<MESH>::SurfaceMesh &surfaceMesh,
                  const double dt, 
                  const double factor,
                  const typename MESH::Node::VecDim &centroid )
{
    typedef typename base::cut::SurfaceMeshBinder<MESH>::SurfaceMesh SurfaceMesh;
    typedef typename base::Vector<MESH::Node::dim>::Type             VecDim;
    
    // go through all elements of the surface mesh
    typename SurfaceMesh::ElementPtrIter eIter = surfaceMesh.elementsBegin();
    typename SurfaceMesh::ElementPtrIter eEnd  = surfaceMesh.elementsEnd();
    for ( ; eIter != eEnd; ++eIter )
    {
        const std::size_t elemID = (*eIter) -> getID();
            
        const typename MESH::Element*  geomElem  = mesh.elementPtr( elemID );
        const typename FIELD::Element* velocElem = velocity.elementPtr( elemID );

        typename SurfaceMesh::Element::ParamIter   pIter = (*eIter) -> parametricBegin();
        typename SurfaceMesh::Element::ParamIter   pEnd  = (*eIter) -> parametricEnd();
        typename SurfaceMesh::Element::NodePtrIter nIter = (*eIter) -> nodesBegin();
        for ( ; pIter != pEnd; ++pIter, ++nIter ) {

            // evaluate time step times velocity at parameter coordinate
            const VecDim du = dt * 
                base::post::evaluateFieldHistory<0>( geomElem, velocElem, *pIter );

            // get old location
            VecDim x;
            (*nIter) -> getX( &(x[0]) );
            
            // add value to nodal coordinate
            //const VecDim xNew = x + du;
            const VecDim xNew = centroid + factor * (x-centroid) + du;
            
            (*nIter) -> setX( &(xNew[0]) );
        }
    }

    return;
}

//------------------------------------------------------------------------------
// Rescale surface coordinates in order to maintain a given volume
template<typename SMESH>
void rescaleSurface( SMESH &surfaceMesh,
                     const double initialVolume,
                     const double enclosed, 
                     const typename SMESH::Node::VecDim &centroid )
{
    typedef typename base::Vector<SMESH::Node::dim>::Type             VecDim;

    const double factor =
        std::pow( initialVolume/enclosed,
                  1./static_cast<double>( SMESH::Node::dim ) );
    
    // go through all elements of the surface mesh
    typename SMESH::ElementPtrIter eIter = surfaceMesh.elementsBegin();
    typename SMESH::ElementPtrIter eEnd  = surfaceMesh.elementsEnd();
    for ( ; eIter != eEnd; ++eIter )
    {
        typename SMESH::Element::ParamIter   pIter = (*eIter) -> parametricBegin();
        typename SMESH::Element::ParamIter   pEnd  = (*eIter) -> parametricEnd();
        typename SMESH::Element::NodePtrIter nIter = (*eIter) -> nodesBegin();
        for ( ; pIter != pEnd; ++pIter, ++nIter ) {

            // get old location
            VecDim x;
            (*nIter) -> getX( &(x[0]) );
            
            // add value to nodal coordinate
            const VecDim xNew = centroid + factor * (x-centroid);
            
            (*nIter) -> setX( &(xNew[0]) );
        }
    }

    return;
}

#endif
