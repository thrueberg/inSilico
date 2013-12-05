#ifndef findlocation_hpp
#define findlocation_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <algorithm>
#include <utility>
// boost includes
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
// base includes
#include <base/linearAlgebra.hpp>
#include <base/geometry.hpp>
#include <base/post/evaluateField.hpp>
#include <solid/Deformation.hpp>

//------------------------------------------------------------------------------
template<typename MESH,typename FIELD>
bool findPreviousLocationInMesh( const MESH&  mesh,
                                 const FIELD& displacement,
                                 const typename MESH::Node::VecDim& x,
                                 const double tolerance, const unsigned maxIter,
                                 std::pair<std::size_t,
                                 typename MESH::Element::GeomFun::VecDim>& result );

//------------------------------------------------------------------------------
template<typename GEOMELEM, typename FIELDELEM>
std::pair<typename GEOMELEM::GeomFun::VecDim,bool>
findPreviousLocationInElement( const GEOMELEM* geomEp, const FIELDELEM* fieldEp,
                               const typename GEOMELEM::Node::VecDim& x,
                               const double tolerance, const unsigned maxIter );

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
void findPreviousDoFLocations( const MESH& mesh, const FIELD& displacement,
                               const std::vector<double>& supports,
                               const std::vector<std::pair<std::size_t,
                               typename MESH::Element::GeomFun::VecDim> >&
                               doFLocation,
                               std::vector<std::pair<std::size_t,
                               typename MESH::Element::GeomFun::VecDim> >&
                               previousDoFLocation,
                               const double supportThreshold,
                               const double tolerance,
                               const unsigned maxIter )
{
    typedef base::GeomTraits<typename MESH::Element> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;
    
    std::vector<bool> marker( doFLocation.size(), false );
    previousDoFLocation.resize( doFLocation.size() );

    typename FIELD::DoFPtrConstIter dIter = displacement.doFsBegin();
    typename FIELD::DoFPtrConstIter dEnd  = displacement.doFsEnd();
    for ( ; dIter != dEnd; ++dIter ) {

        const std::size_t doFID = (*dIter) -> getID();

        // nothing to be done for inactive DoF
        if (  supports[doFID] < supportThreshold ) marker[doFID] = true;

        // only if DoF has not yet been considered
        if ( not marker[ doFID ] ) {

            const std::size_t elemID = doFLocation[ doFID ].first;
            const LocalVecDim currXi = doFLocation[ doFID ].second;
            const GlobalVecDim currX =
                base::Geometry<typename MESH::Element>()( mesh.elementPtr( elemID ),
                                                          currXi );

            std::pair<std::size_t,LocalVecDim> prev;
            const bool found =
                findPreviousLocationInMesh( mesh, displacement,
                                            currX, tolerance, maxIter, prev );

            VERIFY_MSG( found, "Could not find location" );

            previousDoFLocation[ doFID ] = prev;
        } // check if already considered
    } // loop over DoFs

    return;
}

//------------------------------------------------------------------------------
namespace detail_{
    
    //! Helper to sort the vector of pairs of distance and element pointer
    struct CompareDistancePairs
        : public boost::function<bool (const std::pair<double,std::size_t>&,
                                       const std::pair<double,std::size_t>& )>
    {
        bool operator()( const std::pair<double,std::size_t>& a,
                         const std::pair<double,std::size_t>& b ) const
        {
            return (a.first < b.first);
        }
    };
    
} // end detail_

//------------------------------------------------------------------------------
template<typename MESH, typename FIELD>
bool findPreviousLocationInMesh( const MESH&  mesh,
                                 const FIELD& displacement,
                                 const typename MESH::Node::VecDim& x,
                                 const double tolerance, const unsigned maxIter,
                                 std::pair<std::size_t, 
                                 typename MESH::Element::GeomFun::VecDim> & result )
{
    typedef typename MESH::Element GeomElement;
    typedef typename FIELD::Element FieldElement;
    
    typedef typename base::GeomTraits<GeomElement> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;

    // collect elements with the distance to centroid
    typedef std::vector< std::pair<double,std::size_t> > DistanceElementPairVector;
    DistanceElementPairVector distanceElementPairs;
    typename MESH::ElementPtrConstIter geIter = mesh.elementsBegin();
    typename MESH::ElementPtrConstIter geLast = mesh.elementsEnd();
    typename FIELD::ElementPtrConstIter feIter = displacement.elementsBegin();
    for ( ; geIter != geLast; ++geIter, ++feIter ) {

        const std::size_t elemID = (*geIter) -> getID();

        const GlobalVecDim centroid =
            base::Geometry<GeomElement>()( *geIter,
                                           base::ShapeCentroid<GeomElement::shape>::apply() );

        const GlobalVecDim u =
            base::post::evaluateField( *geIter, *feIter,
                                       base::ShapeCentroid<GeomElement::shape>::apply() );

        const double distance = ( (centroid+u) - x).norm();
            
        distanceElementPairs.push_back( std::make_pair( distance, elemID ) );
    }

    // sort this vector
    ::detail_::CompareDistancePairs cdp;
    std::sort( distanceElementPairs.begin(), distanceElementPairs.end(), cdp );

    // go from beginning to end through this vector
    typename DistanceElementPairVector::iterator fIter = distanceElementPairs.begin();
    typename DistanceElementPairVector::iterator fLast = distanceElementPairs.end();
    for ( ; fIter != fLast; ++fIter ) {

        const std::size_t elemID = fIter -> second;
        const GeomElement*  geomEp  = mesh.elementPtr( elemID );
        const FieldElement* fieldEp = displacement.elementPtr( elemID );

        // check this element
        const std::pair<LocalVecDim,bool> trial =
            findPreviousLocationInElement( geomEp, fieldEp, x,
                                           tolerance, maxIter );

        // in case of success return element pointer and local coordinate
        if ( trial.second ) {
            result = std::make_pair( elemID, trial.first );
            return true;
        }
    }

    // if this point is reached, the point has not been found in the entire mesh
    result = std::make_pair( base::invalidInt,
                             base::invalidVector<GT::localDim>() );
    return false;
}

//------------------------------------------------------------------------------
template<typename GEOMELEM, typename FIELDELEM>
std::pair<typename GEOMELEM::GeomFun::VecDim,bool>
findPreviousLocationInElement( const GEOMELEM* geomEp,
                               const FIELDELEM* fieldEp,
                               const typename GEOMELEM::Node::VecDim& x,
                               const double tolerance,
                               const unsigned maxIter )
{
    typedef typename base::GeomTraits<GEOMELEM> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;

    // initial guess
    LocalVecDim xi = base::ShapeCentroid<GEOMELEM::shape>::apply();

    // flag if point has been found
    bool success = false;

    // store previous residual to check if it increases
    double prevResidual = std::numeric_limits<double>::max();

    // Newton iteration
    for ( unsigned iter = 0; iter < maxIter; iter++ ) {

        // Right hand side
        const GlobalVecDim xOfXi = base::Geometry<GEOMELEM>()( geomEp, xi );
        const GlobalVecDim du     =
            base::post::evaluateFieldHistory<0>( geomEp, fieldEp, xi ) -
            base::post::evaluateFieldHistory<1>( geomEp, fieldEp, xi );
        const GlobalVecDim rhs = x - xOfXi - du;

        // already close enough, quit
        const double residual = rhs.norm();
        if ( residual < tolerance )
        {
            success = true;
            break;
        }

        // after two iterations, the residual shall not grow, otherwise quit
        if ( ( iter > 1 ) and ( residual > prevResidual ) )
        {
            success = false;
            break;
        }

        // store new residual
        prevResidual = residual;

        // Ideally, get here the contra-variant basis of the current configuration

        // Jacobi matrix
        typedef typename base::JacobiMatrix<GEOMELEM>::result_type MatDxD;
        const MatDxD J = base::JacobiMatrix<GEOMELEM>()( geomEp, xi );
        // Displacement gradients
        const MatDxD dU1dX =
            base::post::evaluateFieldGradientHistory<0>( geomEp, fieldEp, xi );
        const MatDxD dU2dX =
            base::post::evaluateFieldGradientHistory<1>( geomEp, fieldEp, xi );
        
        // combine and invert
        MatDxD aux;
        for ( unsigned d1 = 0; d1 < GT::globalDim; d1++ )
            for ( unsigned d2 = 0; d2 < GT::localDim; d2++ )
                aux(d1,d2) = (d1==d2 ? 1. : 0.) + dU1dX(d1,d2) - dU2dX(d1,d2);
        aux *= J;

        // solve system
        LocalVecDim dXi;
        dXi.noalias() = aux.inverse() * rhs;

        // check if inside the reference element, if not quit
        if ( not base::InsideShape<GEOMELEM::shape>::apply( xi + dXi, tolerance ) )
        {
            success = false;
            break;
        }

        // update the coordinate
        xi += dXi;

    }

    // return the latest local coordinate and a success flag
    return std::make_pair( xi, success );
}

#endif
