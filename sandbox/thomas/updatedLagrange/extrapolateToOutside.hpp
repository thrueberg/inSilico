#ifndef extrapolatetofictitious_h
#define extrapolatetofictitious_h

#include <vector>
#include <utility>

#include <base/linearAlgebra.hpp>
#include <base/verify.hpp>
#include <base/geometry.hpp>

#include <base/cut/LevelSet.hpp>
#include <base/post/findLocation.hpp>


namespace detail_{

    //--------------------------------------------------------------------------
    // Set history value of DoF to extrapolated
    template<int H>
    struct SetExtrapolatedHistoryValue
    {
        template<typename GEOMELEM, typename FIELDELEM>
        static void apply( const GEOMELEM*  geomEPtr1,  const GEOMELEM*  geomEPtr2,
                           const FIELDELEM* fieldEPtr1, const FIELDELEM* fieldEPtr2,
                           const typename FIELDELEM::FEFun::VecDim& xi1,
                           const typename FIELDELEM::FEFun::VecDim& xi2,
                           const double factor, 
                           typename FIELDELEM::DegreeOfFreedom* doFPtr )
        {
            typedef typename
                base::Vector<FIELDELEM::DegreeOfFreedom::size,
                             base::number>::Type VecDoF;

            // Value at closest point
            const VecDoF uClosest =
                base::post::evaluateFieldHistory<H>( geomEPtr1, fieldEPtr1, xi1 );
            // Value at ghost point
            const VecDoF uGhost   =
                base::post::evaluateFieldHistory<H>( geomEPtr2, fieldEPtr2, xi2 );

            // Extrapolated value
            const VecDoF uExtra =
                uClosest + factor * (uGhost - uClosest);

            // Set dof to value
            for ( unsigned d = 0; d < FIELDELEM::DegreeOfFreedom::size; d++ )
                doFPtr ->template setHistoryValue<H>( d, uExtra[d] );

            // call recursively in order to extrapolate the entire history
            SetExtrapolatedHistoryValue<H-1>::apply( geomEPtr1, geomEPtr2,
                                                     fieldEPtr1, fieldEPtr2,
                                                     xi1, xi2, factor,
                                                     doFPtr );
        }
    };

    //--------------------------------------------------------------------------
    // Truncated template recursion
    template<>
    struct SetExtrapolatedHistoryValue<-1>
    {
        template<typename GEOMELEM, typename FIELDELEM>
        static void apply( const GEOMELEM*  geomEPtr1,  const GEOMELEM*  geomEPtr2,
                           const FIELDELEM* fieldEPtr1, const FIELDELEM* fieldEPtr2,
                           const typename FIELDELEM::FEFun::VecDim& xi1,
                           const typename FIELDELEM::FEFun::VecDim& xi2,
                           const double factor, 
                           typename FIELDELEM::DegreeOfFreedom* doFPtr )
        {
            return;
        }        
    };
}

//------------------------------------------------------------------------------
/** The degrees of freedom on the ficitious side of the domain are given
 *  values that are extrapolated from the physical side of the domain.
 *  Here, a linear extrapolation is carried out by using the solution at the
 *  closest point on the domain boundary and the solution at an additional
 *  'ghost' point in the domain interior. Extrapolation of these two values
 *  gives the value to which the inactive DoF of the ficitious domain is set.
 *  \tparam MESH  Type of geometry representation
 *  \tparam FIELD Type of field to extrapolate
 */
template<typename MESH, typename FIELD>
void extrapolateToFictitious(
    const MESH& mesh, FIELD& field,
    const double length,
    const std::vector< std::pair<std::size_t,
    typename FIELD::Element::FEFun::VecDim> >&
    doFLocation,
    const std::vector<base::cut::LevelSet<MESH::Node::dim> >& levelSet,
    const BoundingBox<MESH::Node::dim>& bbox,
    const bool outside = true )
{
    // Convenience typedefs
    typedef base::GeomTraits<typename MESH::Element> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;

    // go through field's dofs
    typename FIELD::DoFPtrIter fIter = field.doFsBegin();
    typename FIELD::DoFPtrIter fEnd  = field.doFsEnd();
    for ( ; fIter != fEnd; ++fIter ) {
        
        const std::size_t doFID  = (*fIter) -> getID();
        const std::size_t elemID = doFLocation[ doFID ].first;
        const LocalVecDim xi     = doFLocation[ doFID ].second;
        const typename MESH::Element* geomEp = mesh.elementPtr( elemID );

        // compute signed distance for specific DoF
        const double signedDist =
            base::cut::signedDistance( geomEp, xi, levelSet );

        // decide if that DoF lies in the fictitious domain
        const bool isFictitious =
            ( outside ? (signedDist < 0.) : (signedDist > 0. ) );

        // Inactive DoFs with are given extrapolated values
        if ( not( (*fIter) -> isActive(0) ) and isFictitious ){
        
            // location of the DoF
            const GlobalVecDim doFPoint =
                base::Geometry<typename MESH::Element>()( geomEp, xi );

            // location of closest point on surface
            GlobalVecDim closestPoint =
                base::cut::closestPoint( geomEp, xi, levelSet );
            closestPoint = bbox.snap( closestPoint );
                    
            // ghost point used for extrapolation
            GlobalVecDim ghostPoint = closestPoint +
                (length / signedDist) * (doFPoint - closestPoint);
            ghostPoint = bbox.snap( ghostPoint );

            const double newLength = (ghostPoint - closestPoint).norm();

            // find closest and ghost points in the mesh
            std::pair<std::size_t,LocalVecDim> aux1, aux2;
            
            const bool found1 = 
                base::post::findLocationInMesh( mesh, closestPoint, 1.e-6, 10, aux1 );
            const bool found2 = 
                base::post::findLocationInMesh( mesh, ghostPoint,   1.e-6, 10, aux2 );
            
            VERIFY_MSG( found1,
                        x2s("Cannot find closest point ") +
                        x2s( closestPoint.transpose()) +
                        x2s(" in mesh" ) );
                    
            VERIFY_MSG( found2,
                        x2s("Cannot find ghost point ") +
                        x2s( ghostPoint.transpose()) +
                        x2s(" in mesh" ) );


            // evaluate data at closest point
            const typename MESH::Element*  gElem1 = mesh.elementPtr(  aux1.first );
            const typename FIELD::Element* fElem1 = field.elementPtr( aux1.first );
                    
            // evaluate data at ghost point
            const typename MESH::Element*  gElem2 = mesh.elementPtr(  aux2.first );
            const typename FIELD::Element* fElem2 = field.elementPtr( aux2.first );

            // check location of ghost point
            const double distGhost =
                base::cut::signedDistance( gElem2, aux2.second, levelSet );
            const bool ghostIsFictitious =
                ( outside ? (distGhost < 0.) : (distGhost > 0.) );

            // set multiplier to zero for ghost point in ficitious domain
            const double factor = ( ghostIsFictitious ? 0. : (signedDist/newLength) );

            // Call Helper object which sets the entire DoF history
            detail_::SetExtrapolatedHistoryValue<FIELD::DegreeOfFreedom::nHist>::apply(
                gElem1, gElem2, fElem1, fElem2, aux1.second, aux2.second,
                factor, *fIter );

            
        } // check if inactive and outside
        
    } // loop over field dofs

    return;
}


#endif
