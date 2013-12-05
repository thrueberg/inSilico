template<typename MESH, typename FIELD>
void moveSurface( const MESH& mesh,
                  const FIELD& displacement,
                  typename base::cut::SurfaceMeshBinder<MESH>::SurfaceMesh &surfaceMesh )
{
    typedef typename base::cut::SurfaceMeshBinder<MESH>::SurfaceMesh SurfaceMesh;
    typedef typename base::Vector<MESH::Node::dim>::Type             VecDim;
    
    // go through all elements of the surface mesh
    typename SurfaceMesh::ElementPtrIter eIter = surfaceMesh.elementsBegin();
    typename SurfaceMesh::ElementPtrIter eEnd  = surfaceMesh.elementsEnd();
    for ( ; eIter != eEnd; ++eIter )
    {
        const std::size_t elemID = (*eIter) -> getID();
            
        const typename MESH::Element*  geomElem = mesh.elementPtr( elemID );
        const typename FIELD::Element* dispElem = displacement.elementPtr( elemID );

        typename SurfaceMesh::Element::ParamIter   pIter = (*eIter) -> parametricBegin();
        typename SurfaceMesh::Element::ParamIter   pEnd  = (*eIter) -> parametricEnd();
        typename SurfaceMesh::Element::NodePtrIter nIter = (*eIter) -> nodesBegin();
        for ( ; pIter != pEnd; ++pIter, ++nIter ) {

            // evaluate displacement at parameter coordinate
            const VecDim u =
                base::post::evaluateFieldHistory<0>( geomElem, dispElem, *pIter ) -
                base::post::evaluateFieldHistory<1>( geomElem, dispElem, *pIter );

            // add value to nodal coordinate
            VecDim X;
            (*nIter) -> getX( &(X[0]) );
            
            const VecDim x = X + u;
            (*nIter) -> setX( &(x[0]) );
        }
    }

    return;
}
