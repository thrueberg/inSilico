//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   TransferSurfaceDatum.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_transfersurfacedatum_hpp
#define base_cut_transfersurfacedatum_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
// boost includes
#include <boost/function.hpp>
// base  includes
#include <base/linearAlgebra.hpp>
// base/cut includes
#include <base/cut/LevelSet.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{


        template<typename SURFACEMESH, typename SURFACEFIELD,
                 typename DOMAINELEMENT>
        class TransferSurfaceDatum;

    }
}

//------------------------------------------------------------------------------
/** Transfer a datum from an immersed surface to the domain.
 *  In Nitsche's method, or a pure penalty method, the Dirichlet datum of an
 *  immersed surface is applied in a weak form, leading to the right-hand-side
 *  terms
 *  \f[
 *       -\int_\Gamma B(v) g d s + \gamma \int_\Gamma v g d s
 *  \f]
 *  the former being the variationally consistent part of Nitsche's  method and
 *  the latter the penalty or stabilisation term. Moreover, \f$ v \f$ is the
 *  test field, \f$ B \f$ a model-dependent boundary operator, and \f$ g \f$ is
 *  the Dirichlet datum of the surface.
 *
 *  In order to decouple the domain and surface solvers, the given datum is
 *  approximated by
 *  \f[
 *      \tilde{g}(x) = \sum_K \phi_K(x) g(x^*_K)
 *  \f]
 *  where \f$ \phi^K \f$ are the shape functions of the geometry approximation
 *  and \f$ x^*_K \f$ is the surface point closest to the node \f$ x_K \f$.
 *
 *  This object is given the surface mesh and the field (living on the surface)
 *  which is considered the Dirichlet datum for the domain. Using the level
 *  set data for every domain node, the above interpolation can be realised.
 *
 *  \tparam SURFACEMESH   Type of surface mesh
 *  \tparam SURFACEFIELD  Type of field on surface (Dirichlet datum)
 *  \tparam DOMAINELEMENT Type of element in the domain
 */
template<typename SURFACEMESH, typename SURFACEFIELD, typename DOMAINELEMENT>
class base::cut::TransferSurfaceDatum
    : public boost::function<
    typename base::Vector<SURFACEFIELD::DegreeOfFreedom::
                          size>::Type( const DOMAINELEMENT*,
                                       const typename DOMAINELEMENT::GeomFun::VecDim& )>
        
{
public:
    //! @name Template parameter
    //@{
    typedef SURFACEMESH   SurfaceMesh;
    typedef SURFACEFIELD  SurfaceField;
    typedef DOMAINELEMENT DomainElement;
    //@}

    //! @name Used attributes
    //@{
    static const unsigned dim     = DOMAINELEMENT::dim;
    static const unsigned doFSize = SurfaceField::DegreeOfFreedom::size;
    //@}

    //! @name Convenience typedefs
    //@{
    typedef base::cut::LevelSet<dim>                      LevelSet;
    typedef typename base::Vector<doFSize>::Type          VecDof;
    typedef typename DomainElement::GeomFun::VecDim       DomainVecDim;
    typedef typename SurfaceField::Element::FEFun::VecDim SurfaceVecDim;
    //@}

    //! Constructor to set the object's references
    TransferSurfaceDatum( const SurfaceMesh&  surfaceMesh,
                          const SurfaceField& surfaceField,
                          const std::vector<LevelSet>& levelSet )
        : surfaceMesh_(  surfaceMesh ),
          surfaceField_( surfaceField ),
          levelSet_(     levelSet ) { }


    //--------------------------------------------------------------------------
    /** Main function call: interpolate surface datum inside given element.
     *
     *  \param[in] domainEp Pointer to domain element
     *  \param[in] xi       Local interpolation coordinate in element
     *  \return             Interpolated surface datum
     */
    VecDof operator()( const DomainElement* domainEp,
                       const DomainVecDim&  xi ) const
    {
        // result container
        VecDof result = base::constantVector<doFSize>( 0. );

        // evaluate geometry function
        typename DomainElement::GeomFun::FunArray geomFun;
        (domainEp ->  geomFun()).evaluate( domainEp, xi, geomFun );


        // iterate over domain element's nodes
        typename DomainElement::NodePtrConstIter nodeIter =
            domainEp -> nodesBegin();
        typename DomainElement::NodePtrConstIter nodeEnd =
            domainEp -> nodesEnd();
        for ( unsigned n = 0; nodeIter != nodeEnd; ++nodeIter, n++ ) {

            // get node ID
            const std::size_t nodeID = (*nodeIter) -> getID();

            // get closest element ID
            const std::size_t elemID = levelSet_[nodeID].getClosestElement();

            // get closest element's local coordinate
            const SurfaceVecDim eta  = levelSet_[nodeID].getClosestLocalCoordinate();

            // evaluate surface datum
            const VecDof aux =
                base::post::evaluateField( surfaceMesh_.elementPtr(  elemID ),
                                           surfaceField_.elementPtr( elemID ),
                                           eta );

            // linear combination
            result += aux * geomFun[n];

        }
        return result;
    }


private:
    const SurfaceMesh&           surfaceMesh_;  //!< Surface mesh
    const SurfaceField&          surfaceField_; //!< Surface field 
    const std::vector<LevelSet>& levelSet_;     //!< Level set data

};


#endif
