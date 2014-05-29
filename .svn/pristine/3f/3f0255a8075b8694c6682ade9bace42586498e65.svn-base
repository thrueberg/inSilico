//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ComputeSurfaceForces.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_computesurfaceforces_hpp
#define base_cut_computesurfaceforces_hpp

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


        template<typename SURFACEMESH, typename FORCEFIELD,
                 typename SURFACEQUADRATURE, typename SURFFIELDTUPLE,
                 typename TRACTIONFUN>
        class ComputeSurfaceForces;

        namespace detail_{
            
            //------------------------------------------------------------------
            //! General case: multiplication of a tensor and the normal vector
            template<unsigned DIM, unsigned DOFSIZE>
            struct TractionTraits
            {
                typedef typename base::Vector<DIM>::Type     VecDim;
                typedef typename base::Vector<DOFSIZE>::Type VecDoF;

                template<typename TENSOR>
                static VecDoF apply( const TENSOR& tensor,
                                     const VecDim& normal )
                {
                    VecDoF traction = base::constantVector<DOFSIZE>( 0. );
                    for ( unsigned i = 0; i < DOFSIZE; i++ )
                        for ( unsigned j = 0; j < DIM; j++ )
                            traction[i] += tensor( i, j ) * normal[j];

                    return traction;
                }
            };
            
            //------------------------------------------------------------------
            //! Scalar case: dot-product of gradient with normal vector
            template<unsigned DIM>
            struct TractionTraits<DIM,1>
            {
                typedef typename base::Vector<DIM>::Type  VecDim;
                typedef typename base::Vector<1>::Type    VecDoF;

                template<typename TENSOR>
                static VecDoF apply( const TENSOR& tensor,
                                     const VecDim& normal )
                {
                    VecDoF traction = base::constantVector<1>( 0. );
                    for ( unsigned j = 0; j < DIM; j++ )
                        traction[0] += tensor[j] * normal[j];
                    
                    return traction;
                }
            };

            //------------------------------------------------------------------
            //! Define the type of function, the TractionKernel constitutes
            template<typename SURFFIELDTUPLE, unsigned DOFSIZE>
            struct FunTraits
            {
                //! local evaluation coordinate
                typedef typename
                base::GeomTraits<typename SURFFIELDTUPLE::GeomElement>::LocalVecDim
                LocalVecDim;

                //! Result type
                typedef typename base::Vector<DOFSIZE>::Type VecDoF;

                //! Function type
                typedef
                boost::function<void( const SURFFIELDTUPLE&,
                                      const LocalVecDim&, const double,
                                      std::vector<VecDoF>& ) > Type;
                
            };

            //------------------------------------------------------------------
            template<unsigned IDIM, unsigned LDIM>
            struct SurfaceFieldPolicy
            {
                typedef typename base::Vector<IDIM>::Type VecIDim;
                typedef typename base::Vector<LDIM>::Type VecLDim;
                
                template<typename SELEMENT>
                static VecIDim mapCoordinate( const SELEMENT* sElemPtr, const VecLDim eta )
                {
                    return sElemPtr -> localDomainCoordinate( eta );
                }

                template<typename SELEMENT>
                static typename SELEMENT::DomainElement* geomElemPtr( SELEMENT* sElemPtr )
                {
                    return sElemPtr -> getDomainElementPointer();
                }
            };

            template<unsigned DIM>
            struct SurfaceFieldPolicy<DIM,DIM>
            {
                typedef typename base::Vector<DIM>::Type   VecDim;
                
                template<typename SELEMENT>
                static VecDim mapCoordinate( const SELEMENT* sElemPtr, const VecDim eta )
                {
                    return eta;
                }

                template<typename SELEMENT>
                static SELEMENT* geomElemPtr( SELEMENT* sElemPtr )
                {
                    return sElemPtr;
                }

            };

            //------------------------------------------------------------------
            //! Kernel to compute the nodal domain forces
            template<typename SURFFIELDTUPLE, typename TRACTIONFUN,
                     unsigned DOFSIZE>
            class TractionKernel
                : public FunTraits<SURFFIELDTUPLE,DOFSIZE>::Type
            {                                   
            public:
                //! @name Template parameter
                //@{
                typedef SURFFIELDTUPLE    SurfaceFieldTuple;
                typedef TRACTIONFUN       TractionFun;
                static const unsigned doFSize = DOFSIZE;
                //@}

                //! Type of surface element
                typedef typename SurfaceFieldTuple::GeomElement SurfaceElement;
                //! Type of domain element
                typedef typename SurfaceElement::DomainElement DomainElement;

                //! Local evaluation coordinate
                typedef typename base::GeomTraits<SurfaceElement>::LocalVecDim
                LocalVecDim;

                //! spatial dimension
                static const unsigned dim =
                    base::GeomTraits<SurfaceElement>::globalDim;
                
                //! Local evaluation coordinate
                typedef typename base::GeomTraits<SurfaceElement>::GlobalVecDim
                GlobalVecDim;

                //! Traction vector
                typedef typename base::Vector<doFSize>::Type VecDoF;

                //! The type of domain tuple used for the stress computation
                typedef typename
                base::asmb::DomainFieldElementPointerTuple<SurfaceFieldTuple>::Type
                DomainFieldTuple;

                //! Constructor with access to stress evaluation function
                TractionKernel( const TractionFun& tractionFun )
                    : tractionFun_( tractionFun ) { }
                
                //--------------------------------------------------------------
                //! Function call operator given to quadrature
                void operator()( const SurfaceFieldTuple&   sft,
                                 const LocalVecDim&         eta,
                                 const double               weight,
                                 std::vector<VecDoF>&       result ) const
                {
                    // extract surface geometry element of the tuple
                    const SurfaceElement* surfEp  = sft.geomElementPtr();

                    // Get surface metric and normal
                    GlobalVecDim normal;
                    const double detG =
                        base::SurfaceNormal<SurfaceElement>()( surfEp,
                                                               eta, normal );
                    
                    // Get local domain coordinate
                    typename DomainElement::GeomFun::VecDim xi =
                        surfEp -> localDomainCoordinate( eta );

                    // Pointer to domain element
                    const DomainElement* domainEp =
                        surfEp -> getDomainElementPointer();

                    // Evaluate the shape functions of the domain element
                    typename DomainElement::GeomFun::FunArray geomFunValues;
                    (domainEp -> geomFun()).evaluate( domainEp, xi, geomFunValues );
                    
                    // make tuple out of domain field
                    const DomainFieldTuple dft =
                        base::asmb::DomainFieldElementPointerTuple<SurfaceFieldTuple>::
                        convert( sft );

                    // compute traction by multiplication with normal vector
                    const VecDoF traction = tractionFun_( dft, xi, normal );
                    
                    // add to result container the weighted tractions
                    for ( unsigned s = 0; s < geomFunValues.size(); s++ )
                        result[s] += traction * weight * detG * geomFunValues[s];
                }
                
            private:
                const TractionFun& tractionFun_;
            };
            
        }// end namespace detail_
        
    }
}

//------------------------------------------------------------------------------
/**  Computation of consistent surface forces.
 *   In order to ensure stability of a fluid-structure coupling the transfer of
 *   forces from the domain to the surface (fluid to structure) has to have
 *   mapping properties which are ideally the transpose of the transfer of the
 *   surface datum to the domain (structure to fluid). Here, the approach as
 *   proposed by Farhat, Lesoinne and LeTalle, CMAME 1998, is implemented.
 *   This approache consists of two steps: (1) computation of local reaction
 *   forces in the fluid and (2) transfer of these forces to the surface.
 *
 *   For the first step, the following forces are computed
 *   \f[
 *       F_K = \int_{\Gamma^f} t(x) \phi^f_K (x) d s
 *   \f]
 *   with the traction function \f$ t(x) = \sigma \cdot n \f$ and the shape
 *   function of the domain element. The superscripts \f$ f \f$ indicate that
 *   the interface and the shape functions are considered on the fluid side
 *   (or in general the domain side) only. Now, every fluid 'node' is equipped
 *   with a nodal force \f$ F_K \f$.
 *
 *   In the second step, these nodal forces are transferred to the surface
 *   (structure \f$ s \f$) by using the surface shape function values at the
 *   surface points closest to the domain nodes. Therefore, the surface node
 *   \f$ j \f$ receives the force
 *   \f[
 *        F_j = \sum F_K \phi^s_j (x^*_K)
 *   \f]
 *   where \f$ x^*_K \f$ is the surface point closest to the K-th node at
 *   \f$ x_K \f$.
 *
 *   \tparam SURFACEMESH        Mesh of the surface of immersed domain
 *   \tparam FORCEFIELD         Field receiving the computed forces
 *   \tparam SURFACEQUADRATURE  Quadrature along the domain interface
 *   \tparam SURFACEFIELDTUPLE  Tuple of surface mesh and field
 *                              element pointers
 *   \tparam TRACTIONFUN        Evaluation function for the traction
 */
template<typename SURFACEMESH, typename FORCEFIELD,
         typename SURFACEQUADRATURE, typename SURFFIELDTUPLE,
         typename TRACTIONFUN>
class base::cut::ComputeSurfaceForces
    : boost::noncopyable
{
public:
    //! @name Template parameter
    //@{
    typedef SURFACEMESH       SurfaceMesh;
    typedef FORCEFIELD        ForceField;
    typedef SURFACEQUADRATURE SurfaceQuadrature;
    typedef SURFFIELDTUPLE    SurfaceFieldTuple;
    typedef TRACTIONFUN       TractionFun;
    //@}

    //! @name Involved element types
    //@{
    typedef typename SurfaceFieldTuple::GeomElement ImmersedElement;
    typedef typename ImmersedElement::DomainElement DomainElement;
    typedef typename SurfaceMesh::Element           SurfaceElement;
    typedef typename ForceField::Element            ForceElement;
    //@}

    //! @name Attributes
    //@{
    static const unsigned globalDim      = base::GeomTraits<DomainElement>::globalDim;
    static const unsigned localDim       = SurfaceElement::dim;
    static const unsigned immersedDim    = ForceElement::dim;
    static const unsigned numDomainNodes = DomainElement::numNodes;
    static const unsigned doFSize        = ForceElement::DegreeOfFreedom::size;
    //@}

    //!
    typedef typename detail_::SurfaceFieldPolicy<immersedDim,localDim> SFP;

    //! Force vector
    typedef typename base::Vector<doFSize>::Type VecDoF;

    //! Type of level set
    typedef typename base::cut::LevelSet<globalDim> LevelSet;
    
    //! Type of kernel for the traction computation
    typedef detail_::TractionKernel<SurfaceFieldTuple,TractionFun,doFSize>
    TractionKernel;

    typedef base::asmb::SurfaceFieldBinder<SurfaceMesh,ForceField> SFB;

    //! Constructor 
    ComputeSurfaceForces( SurfaceMesh&                 surfaceMesh,
                          ForceField&                  forceField,
                          const SurfaceQuadrature&     surfaceQuadrature,
                          const std::vector<LevelSet>& levelSet,
                          const TractionFun&           tractionFun,
                          const bool plusMinus = true )
        : surfaceFieldBinder_( surfaceMesh, forceField ),
          surfaceQuadrature_( surfaceQuadrature ),
          levelSet_(          levelSet ),
          tractionKernel_(    tractionFun ),
          sign_( plusMinus ? +1.0 : -1.0 )
    {
    }

    //------------------------------------------------------------------------------
    /** Main function call operator which computes the nodal forces.
     *
     *  \tparam surfaceFieldTuple  Tuple of field element pointers
     *  \return                    Sum of computed forces
     */
    VecDoF operator()( const SurfaceFieldTuple& surfaceFieldTuple )
    {
        // result container for users interested in the sum of forces
        VecDoF sumOfForces = base::constantVector<doFSize>( 0. );

        // get domain element in order to access the nodes
        const DomainElement* domainEp =
            (surfaceFieldTuple.geomElementPtr()) -> getDomainElementPointer();
        
        // iterators to nodes of the domain element
        typename DomainElement::NodePtrConstIter nIter = domainEp -> nodesBegin();
        typename DomainElement::NodePtrConstIter nLast = domainEp -> nodesEnd();
        
        // allocate space for nodal forces of domain 
        std::vector<VecDoF> domainNodalForces;
        for ( ; nIter != nLast; ++nIter )
            domainNodalForces.push_back( base::constantVector<doFSize>( 0. ) );

        // use force kernel to compute nodal forces
        surfaceQuadrature_.apply( tractionKernel_, surfaceFieldTuple, domainNodalForces );
        
        // go through the domain element's nodes
        nIter = domainEp -> nodesBegin();
        for ( unsigned n = 0; nIter != nLast; ++nIter, n++ ) {

            // nodal force of this node
            const VecDoF nodalForce = domainNodalForces[n];

            // node ID
            const std::size_t nodeID = (*nIter) -> getID();

            // get level set data: closest element ID and local coordinates
            const std::size_t closestElementID =
                levelSet_[nodeID].getClosestElement();
            const typename LevelSet::LocalVecDim eta =
                levelSet_[nodeID].getClosestLocalCoordinate();

            // evaluate the shape functions of the surface element
            const typename SFB::ElementPtrTuple tmp
                = surfaceFieldBinder_.elementPtr( closestElementID );
            
            const SurfaceElement* sElemPtr = tmp.geomElementPtr();
            const ForceElement*   fElemPtr = tmp.testElementPtr();

            // map coordinate (if necessary)
            const typename ForceElement::FEFun::VecDim xi
                = SFP::mapCoordinate( sElemPtr, eta );
            
            // evaluate the fe-functions of the force element
            typename ForceElement::FEFun::FunArray surfFun;
            (fElemPtr -> fEFun() ).evaluate( SFP::geomElemPtr( sElemPtr ),
                                             xi, surfFun );
                                               
            // go through degrees of freedom of surface element
            typename ForceElement::DoFPtrConstIter dIter = fElemPtr -> doFsBegin();
            typename ForceElement::DoFPtrConstIter dLast = fElemPtr -> doFsEnd();
            for ( unsigned s = 0; dIter != dLast; ++dIter, s++ ) {

                for ( unsigned d = 0; d < doFSize; d++ ) {
                    const number oldValue = (*dIter) -> getValue( d );
                    const number newValue = oldValue + sign_ * surfFun[s] * nodalForce[d];
                    (*dIter) -> setValue( d, newValue );
                }

                sumOfForces += sign_ * surfFun[s] * nodalForce;
            }
        }
        
        return sumOfForces;
    }
    
private:
    SFB                          surfaceFieldBinder_;
    const SurfaceQuadrature&     surfaceQuadrature_;//!< Quadrature at interface
    const std::vector<LevelSet>& levelSet_;         //!< Level set data of nodes
    const TractionKernel         tractionKernel_;   //!< Traction kernel
    const double                 sign_;             //!< Sign change if wanted
};

#endif
