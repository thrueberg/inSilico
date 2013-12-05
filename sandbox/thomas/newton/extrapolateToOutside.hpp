

template<typename MESH, typename FIELD>
void extrapolateToOutside( const MESH& mesh,
                           FIELD& field,
                           const double length,
                           const std::vector< std::pair<std::size_t,
                           typename FIELD::Element::FEFun::VecDim> >&
                           doFLocation,
                           const std::vector<base::cut::LevelSet<MESH::Node::dim> >&
                           levelSet )
{
    typedef base::GeomTraits<typename MESH::Element> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;
    typedef typename base::Vector<FIELD::DegreeOfFreedom::size>::Type VecDoF;

    // go through field's dofs
    typename FIELD::DoFPtrIter fIter = field.doFsBegin();
    typename FIELD::DoFPtrIter fEnd  = field.doFsEnd();
    for ( ; fIter != fEnd; ++fIter ) {
        
        const std::size_t doFID  = (*fIter) -> getID();
        const std::size_t elemID = doFLocation[ doFID ].first;
        const LocalVecDim xi     = doFLocation[ doFID ].second;
        const typename MESH::Element* geomEp = mesh.elementPtr( elemID );
        
        const double signedDist =
            base::cut::signedDistance( geomEp, xi, levelSet );

        // Inactive DoFs with negative distance are given extrapolated values
        if ( not( (*fIter) -> isActive(0) ) and (signedDist < 0.) ){

            // location of the DoF
            const GlobalVecDim doFPoint =
                base::Geometry<typename MESH::Element>()( geomEp, xi );

            // location of closest point on surface
            const GlobalVecDim closestPoint = base::cut::closestPoint( geomEp, xi, levelSet );
                    
            // ghost point used for extrapolation
            const GlobalVecDim ghostPoint = closestPoint +
                (length / signedDist) * (doFPoint - closestPoint);

            // find closest and ghost points in the mesh
            std::pair<std::size_t,LocalVecDim> aux1, aux2;
            
            const bool found1 = 
                base::post::findLocationInMesh( mesh, closestPoint, 1.e-6, 10, aux1 );
            VERIFY( found1 );
                    
            const bool found2 = 
                base::post::findLocationInMesh( mesh, ghostPoint,   1.e-6, 10, aux2 );
            VERIFY( found2 );

            // evaluate data at closest point
            const typename MESH::Element*  gElem1 = mesh.elementPtr(  aux1.first );
            const typename FIELD::Element* fElem1 = field.elementPtr( aux1.first );
                    
            const VecDoF uClosest = base::post::evaluateField( gElem1, fElem1,
                                                               aux1.second );
            // evaluate data at ghost point
            const typename MESH::Element*  gElem2 = mesh.elementPtr(  aux2.first );
            const typename FIELD::Element* fElem2 = field.elementPtr( aux2.first );

            const VecDoF uGhost   = base::post::evaluateField( gElem2, fElem2,
                                                               aux2.second );
            
            // compute extrapolation
            const VecDoF uExtra =
                uClosest + (signedDist/length) * (uGhost - uClosest);

            // set values 
            for ( unsigned d = 0; d < FIELD::DegreeOfFreedom::size; d++ ) {
                
                (*fIter) ->template setHistoryValue<0>( d, uExtra[d] );
                        
            } // all components
            
        } // check if inactive and outside
        
    } // loop over field dofs

    return;
}
