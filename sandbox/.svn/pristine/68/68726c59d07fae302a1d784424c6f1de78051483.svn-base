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
#include <base/auxi/BoundingBox.hpp>

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
                               const double tolerance,
                               const unsigned maxIter,
                               const bool verbose = false );

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
                               const base::auxi::BoundingBox<MESH::Node::dim>& bbox, 
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
        //if ( (*dIter) -> isConstrained(0) )        marker[doFID] = true;

        // only if DoF has not yet been considered
        if ( not marker[ doFID ] ) {

            const std::size_t elemID = doFLocation[ doFID ].first;
            const LocalVecDim currXi = doFLocation[ doFID ].second;
            const GlobalVecDim currX =
                base::Geometry<typename MESH::Element>()( mesh.elementPtr( elemID ),
                                                          currXi );

            // do not search for nodes which are on the box boundary
            if ( bbox.isOnAnyBoundary( currX, 1.e-5 ) ) {
                previousDoFLocation[ doFID ] = doFLocation[ doFID ];
            }
            else {

                std::pair<std::size_t,LocalVecDim> prev;
                const bool found =
                    findPreviousLocationInMesh( mesh, displacement,
                                                currX, tolerance, maxIter, prev );

                VERIFY_MSG( found,
                            x2s( "Could not find location of DoF ") +
                            x2s( doFID ) + x2s( " at x=" ) +
                            x2s( currX.transpose() ) + ", " +
                            x2s( supports[doFID] ) );

                previousDoFLocation[ doFID ] = prev;
            }
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
            base::Geometry<GeomElement>()(
                *geIter,
                base::ShapeCentroid<GeomElement::shape>::apply() );

        const GlobalVecDim du =
            base::post::evaluateFieldHistory<0>(
                *geIter, *feIter,
                base::ShapeCentroid<GeomElement::shape>::apply() );
        
        const double distance = ( (centroid+du) - x).norm();
        
        distanceElementPairs.push_back( std::make_pair( distance, elemID ) );
    }
    
    // sort this vector
    ::detail_::CompareDistancePairs cdp;
    std::sort( distanceElementPairs.begin(), distanceElementPairs.end(), cdp );

    // go from beginning to end through this vector
    const std::size_t numCheck = 10; //( distanceElementPairs.size() );
    for ( std::size_t i = 0; i < numCheck; i++ ) {

        const std::size_t elemID    = distanceElementPairs[i].second;
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

    // Catch
    // { 
    //     result = std::make_pair( distanceElementPairs[0].second,
    //                              base::ShapeCentroid<GeomElement::shape>::apply() );
    //     return true;
    // }
    for ( std::size_t i = 0; i < numCheck; i++ ) {

        const std::size_t elemID    = distanceElementPairs[i].second;
        const GeomElement*  geomEp  = mesh.elementPtr( elemID );
        const FieldElement* fieldEp = displacement.elementPtr( elemID );

        std::cout << "ID: " << distanceElementPairs[i].second
                  << " Dist: " << distanceElementPairs[i].first << "\n";

        // std::cout << " Trying to find (" << x.transpose()
        //           << ") in \n";
        // for ( unsigned i = 0; i < 4; i++ ) {
        //     typename MESH::Node::VecDim x, u;
        //     (*geomEp) -> nodePtr(i) -> getX( &(x[0]) );
        // 
        //     for ( unsigned d = 0; d < 2; d++ )
        //         (*fieldEp) -> doFPtr(i) -> getValue(
        // }
        
        // check this element
        findPreviousLocationInElement( geomEp,
                                       fieldEp,
                                       x,
                                       tolerance, maxIter, true );

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
                               const unsigned maxIter,
                               const bool verbose )
{
    if ( verbose) std::cout << std::endl;
    
    typedef typename base::GeomTraits<GEOMELEM> GT;
    typedef typename GT::GlobalVecDim GlobalVecDim;
    typedef typename GT::LocalVecDim  LocalVecDim;

    // initial guess
    LocalVecDim xi = base::ShapeCentroid<GEOMELEM::shape>::apply();

    // store previous residual to check if it increases
    double prevResidual = std::numeric_limits<double>::max();

    // Newton iteration
    for ( unsigned iter = 0; iter < maxIter; iter++ ) {

        // Right hand side
        const GlobalVecDim xOfXi = base::Geometry<GEOMELEM>()( geomEp, xi );
        const GlobalVecDim du     =
            base::post::evaluateFieldHistory<0>( geomEp, fieldEp, xi );

        const GlobalVecDim rhs = x - xOfXi - du;

        // residual norm
        const double residual = rhs.norm();

        if ( verbose ) {
            std::cout << "Iter-" << iter << ": |(" << rhs.transpose()
                      << ")| = " << residual << std::endl;
        }

        // successful exit
        if ( residual < tolerance )
        {
            return std::make_pair( xi, true );
        }

        // after two iterations, the residual has to decrease
        if ( ( iter > 1 ) and
             ( residual / prevResidual >= 1. - tolerance ) ) {
            return std::make_pair( xi, false );
        }

        // store new residual
        prevResidual = residual;

        // Ideally, get here the contra-variant basis of the current configuration

        // Jacobi matrix
        typedef typename base::JacobiMatrix<GEOMELEM>::result_type MatDxD;
        const MatDxD J = base::JacobiMatrix<GEOMELEM>()( geomEp, xi );
        // Displacement gradient
        const MatDxD dUdX =
            base::post::evaluateFieldGradientHistory<0>( geomEp, fieldEp, xi );
        
        // combine and invert (note the transposed gradients!)
        MatDxD aux;
        for ( unsigned d1 = 0; d1 < GT::globalDim; d1++ )
            for ( unsigned d2 = 0; d2 < GT::localDim; d2++ )
                aux(d1,d2) = (d1==d2 ? 1. : 0.) + dUdX(d2,d1); 
        aux *= J;

        // solve system
        LocalVecDim dXi;
        dXi.noalias() = aux.inverse() * rhs;

        // update the coordinate
        xi += dXi;

        // range check: if new coordinate is 'miles away', quit
        if ( not base::InsideShape<GEOMELEM::shape>::apply( xi, 1.0 ) )
        {
            return std::make_pair( xi, false );
        }
        
        // range chack: if new coordinate is outside, snap to the ref shape
        if ( not base::InsideShape<GEOMELEM::shape>::apply( xi, tolerance ) )
        {
            xi = base::SnapToShape<GEOMELEM::shape>::apply( xi );
        }

    }

    // return the latest local coordinate and a success flag
    return std::make_pair( xi, false );
}

#endif
