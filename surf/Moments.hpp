//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Moments.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef surf_moments_hpp
#define surf_moments_hpp

//------------------------------------------------------------------------------
//
#include <boost/bind.hpp>
// base includes
#include <base/verify.hpp>
#include <base/asmb/SimpleIntegrator.hpp>
#include <base/asmb/SurfaceFieldBinder.hpp>

//------------------------------------------------------------------------------
namespace surf{

    template<typename FIELDTUPLE>
    class Moments;

    //--------------------------------------------------------------------------
    /** Compute the volume enclosed by surface mesh.
     *  \tparam SURFACEMESH Type of surface mesh
     *  \tparam QUADRATURE Surface quadrature
     */
    template<typename SURFACEMESH, typename QUADRATURE>
    double enclosedVolume( const SURFACEMESH& surfaceMesh, 
                           const QUADRATURE&  quadrature )
    {
        // sum up the volume
        double result = 0.;

        typedef base::asmb::SurfaceFieldBinder<const SURFACEMESH>  SurfaceBinder;
        typedef typename SurfaceBinder::template TupleBinder<>::Type STB;
        
        SurfaceBinder surfaceBinder( surfaceMesh );
        // perform integration over surface
        base::asmb::simplyIntegrate<STB>(
            quadrature, result, surfaceBinder,
            boost::bind( &surf::Moments<typename STB::Tuple>::enclosedVolume,
                         _1, _2, _3, _4 ) );
        
        return result;
    }

    //--------------------------------------------------------------------------
    /** Compute the volume moment as enclosed by a surface mesh.
     *  \tparam SURFACEMESH Type of surface mesh
     *  \tparam QUADRATURE Surface quadrature
     */
    template<typename SURFACEMESH, typename QUADRATURE>
    typename base::Vector<SURFACEMESH::Node::dim,double>::Type
    enclosedVolumeMoment( const SURFACEMESH& surfaceMesh, 
                          const QUADRATURE&  quadrature )
    {
        static const unsigned dim = SURFACEMESH::Node::dim;
        // 
        typename base::Vector<dim>::Type moment = base::constantVector<dim>(0.);
        
        typedef base::asmb::SurfaceFieldBinder<const SURFACEMESH>  SurfaceBinder;
        typedef typename SurfaceBinder::template TupleBinder<>::Type STB;
        SurfaceBinder surfaceBinder( surfaceMesh );
        
        base::asmb::simplyIntegrate<STB>(
            quadrature, moment, surfaceBinder,
            boost::bind( &surf::Moments<typename STB::Tuple>::volumeMoment,
                         _1, _2, _3, _4 ) );

        return moment;
    }

    //--------------------------------------------------------------------------
    /** Compute the volume moment as enclosed by a surface mesh.
     *  \tparam SURFACEMESH Type of surface mesh
     *  \tparam QUADRATURE Surface quadrature
     */
    template<typename SURFACEMESH, typename QUADRATURE>
    typename base::Matrix<SURFACEMESH::Node::dim,
                          SURFACEMESH::Node::dim, double>::Type
    enclosedVolumeSecondMoment( const SURFACEMESH& surfaceMesh, 
                                const QUADRATURE&  quadrature )
    {
        static const unsigned dim = SURFACEMESH::Node::dim;
        // 
        typename base::Matrix<dim,dim>::Type moment =
            base::constantMatrix<dim,dim>(0.);
        
        typedef base::asmb::SurfaceFieldBinder<const SURFACEMESH>  SurfaceBinder;
        typedef typename SurfaceBinder::template TupleBinder<>::Type STB;
        SurfaceBinder surfaceBinder( surfaceMesh );

        base::asmb::simplyIntegrate<STB>(
            quadrature, moment, surfaceBinder,
            boost::bind( &surf::Moments<typename STB::Tuple>::secondVolumeMoment,
                         _1, _2, _3, _4 ) );

        return moment;
    }

    //--------------------------------------------------------------------------
    /** Compute the volume and the centroid
     *  \tparam SURFACEMESH Type of surface mesh
     *  \tparam QUADRATURE Surface quadrature
     */
    template<typename SURFACEMESH, typename QUADRATURE>
    double volumeAndCentroid(
        const SURFACEMESH& surfaceMesh, 
        const QUADRATURE&  quadrature,
        typename base::Vector<SURFACEMESH::Node::dim>::Type& centroid )
    {
        const double volume = enclosedVolume( surfaceMesh, quadrature );

        centroid = enclosedVolumeMoment(      surfaceMesh, quadrature );
        centroid /= volume;

        return volume;
    }

    //--------------------------------------------------------------------------
    /** Compute the volume, the centroid and the inertia tensor
     *  \tparam SURFACEMESH Type of surface mesh
     *  \tparam QUADRATURE Surface quadrature
     */
    template<typename SURFACEMESH, typename QUADRATURE>
    double volumeCentroidAndInertia(
        const SURFACEMESH& surfaceMesh, 
        const QUADRATURE& quadrature,
        typename base::Vector<SURFACEMESH::Node::dim>::Type& centroid,
        typename base::Matrix<SURFACEMESH::Node::dim,
        SURFACEMESH::Node::dim>::Type& inertia )
    {
        static const unsigned dim = SURFACEMESH::Node::dim;
        
        const double volume = enclosedVolume( surfaceMesh, quadrature );

        centroid = enclosedVolumeMoment(      surfaceMesh, quadrature );
        centroid /= volume;

        inertia = enclosedVolumeSecondMoment( surfaceMesh, quadrature );
        // shift to centroid
        for ( unsigned d1 = 0; d1 < dim; d1++ ) {
            for ( unsigned d2 = 0; d2 < dim; d2++ ) {
                inertia(d1, d2) -= volume * centroid[d1] * centroid[d2];
            }
        }

        return volume;
    }

    //--------------------------------------------------------------------------
    /** Integrate over the normal-field of the surface, result should be
     *  zero for any closed surface.
     *  \tparam SURFACEMESH Type of surface mesh
     *  \tparam QUADRATURE Surface quadrature
     */
    template<typename SURFACEMESH, typename QUADRATURE>
    double isClosed( const SURFACEMESH& surfaceMesh,
                     const QUADRATURE&  quadrature )
    {
        // sum up the volume
        double result = 0.;

        typedef base::asmb::SurfaceFieldBinder<const SURFACEMESH>  SurfaceBinder;
        typedef typename SurfaceBinder::template TupleBinder<>::Type STB;
        SurfaceBinder surfaceBinder( surfaceMesh );

        // perform integration over surface
        base::asmb::simplyIntegrate<STB>(
            quadrature, result, surfaceBinder,
            boost::bind( &surf::Moments<typename STB::Tuple>::isClosed,
                         _1, _2, _3, _4 ) );
        
        return result;
    }

}

//------------------------------------------------------------------------------
/**  Geometry features of domains entire enclosed by a surface mesh.
 *   Using integration by parts, the enclosed volume or the centroid (which is
 *   the first moment of volume divided by volume) can be computed via surface
 *   integrals. This class provides the methods for these two quantities.
 *   \note Prerequisite is that the surface mesh is \em closed, otherwise the
 *         results are undefined.
 *   \tparam FIELDTUPLE Type of field tuple out of which only the geometry is
 *                      used
 */
template<typename FIELDTUPLE>
class surf::Moments
{
public:
    //! Template parameter
    typedef FIELDTUPLE FieldTuple;

    //! @name Extract element types from pointers
    //@{
    typedef typename FieldTuple::GeomElement  GeomElement;
    //@}

    //! Hack to allow boost::bind use this class
    typedef void result_type;

    static const unsigned dim = base::GeomTraits<GeomElement>::globalDim;

    //! @name Local and global coordinates
    //@{
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim   LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim  GlobalVecDim;
    typedef typename base::Matrix<dim,dim>::Type                  MatDimDim;
    //@}

    //! Global space dimension
    static const unsigned globalDim = base::GeomTraits<GeomElement>::globalDim;

    //--------------------------------------------------------------------------
    /** Computation of the volume enclosed by the mesh.
     *  Using the divergence theorem, one gets
     *  \f[
     *      V = \int_\Omega 1 d x
     *        = \int_\Omega \nabla \cdot r_1 d x
     *        = \int_\Omega x_1 n_1 d s
     *  \f]
     *  where \f$ r_1 \f$ is the vector \f$ (x_1,0,0) \f$ which holds only the
     *  1-component of the position vector \f$ x = (x_1, x_2, x_3) \f$. 
     *  Consequently, \f$ n_1 \f$ is the 1-component of the normal vector
     *  \f$ n \f$. 
     *  \param[in]  fieldTuple Tuple of field element pointers
     *  \param[in]  xi         Local surface evaluation coordinate
     *  \param[in]  weight     Integration weight
     *  \param[out] volume     Final sum equals the enclosed volume
     */
    static void enclosedVolume( const FieldTuple&  fieldTuple,
                                const LocalVecDim& xi,
                                const double       weight,
                                double& volume )
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();
        

        // Get surface normal and measure
        GlobalVecDim normal;
        const double detJ = 
            base::SurfaceNormal<GeomElement>()( geomEp, xi, normal );

        // get global coordinate
        const typename GeomElement::Node::VecDim x =
            base::Geometry<GeomElement>()( geomEp, xi );

        // integration kernel after integration by parts
        volume += detJ * weight * normal[0] *x[0];
    }

    //--------------------------------------------------------------------------
    /** First moment of volume of the domain enclosed by the mesh.
     *  The first moment of volume is defined as
     *  \f[
     *       S = \int_\Omega x d x
     *  \f]
     *  where both \f$ S \f$ and \f$ x \f$ are vectors. Using integration by
     *  parts, we get
     *  \f[
     *       S = \int_\Omega \frac{1}{2} \nabla \cdot X^2 d x
     *         = \int_\Gamma \frac{1}{2} X^2 \cdot n d s
     *         = \frac{1}{2} \int_\Gamma \sum_i n_i x_i^2 d s
     *  \f]
     *  where \f$ X^2 \f$ is a diagonal matrix with entries
     *  \f$ (x_1^2, x_2^2, x_3^2) \f$ along the diagonal. In order to compute
     *  the centroid of the enclosed volume, the vector \f$ S \f$ has to be
     *  divided component-wise by the enclosed volume \f$ V \f$,
     *  \f[
     *        C = \frac{1}{V} S
     *  \f]
     *  \param[in]  fieldTuple Tuple of field element pointers
     *  \param[in]  xi         Local surface evaluation coordinate
     *  \param[in]  weight     Integration weight
     *  \param[out] moment     Final sum equals the moment of volume
     */
    static void volumeMoment( const FieldTuple&  fieldTuple,
                            const LocalVecDim& xi,
                            const double       weight,
                            GlobalVecDim&      moment )
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();

        // Get surface normal and measure
        GlobalVecDim normal;
        const double detJ = 
            base::SurfaceNormal<GeomElement>()( geomEp, xi, normal );

        // get global coordinate
        const GlobalVecDim x = base::Geometry<GeomElement>()( geomEp, xi );

        //
        for ( unsigned d = 0; d < globalDim; d++ )
            moment[d] += 0.5 * x[d] * x[d] * normal[d] * detJ * weight;
    }

    //--------------------------------------------------------------------------
    /** First moment of volume of the domain enclosed by the mesh.
     *  The first moment of volume is defined as
     *  \f[
     *       I_{ij} = \int_\Omega x_i x_j d x
     *  \f]
     *  Using integration by
     *  parts, we get
     *  \f[
     *       I_{ii} = \int_\Omega \frac{1}{3} (x_i^3)_{,i} d x
     *              = \int_\Gamma \frac{1}{3} x_i^3 n_i d s \qquad
     *       I_{ij} = \int_\Omega \frac{1}{2} (x_i^2 x_j)_{,i} d x
     *              = \int_\Gamma \frac{1}{2} x_i^2 x_j n_i d s
     *  \f]
     *  \param[in]  fieldTuple Tuple of field element pointers
     *  \param[in]  xi         Local surface evaluation coordinate
     *  \param[in]  weight     Integration weight
     *  \param[out] moment     Final sum equals the moments of inertia
     */
    static void secondVolumeMoment( const FieldTuple&  fieldTuple,
                                    const LocalVecDim& xi,
                                    const double       weight,
                                    MatDimDim&         moment )
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();

        // Get surface normal and measure
        GlobalVecDim normal;
        const double detJ = 
            base::SurfaceNormal<GeomElement>()( geomEp, xi, normal );

        // get global coordinate
        const GlobalVecDim x = base::Geometry<GeomElement>()( geomEp, xi );

        //
        for ( unsigned d1 = 0; d1 < globalDim; d1++ ){
            for ( unsigned d2 = 0; d2 < globalDim; d2++ ){
                const double factor = (d1==d2 ? 1./3. : 1./2.);
                moment( d1, d2 ) +=
                    factor * x[d1] * x[d1] * x[d2] * normal[d1] * detJ * weight;
            }
        }
    }

    //--------------------------------------------------------------------------
    static void isClosed( const FieldTuple&  fieldTuple,
                          const LocalVecDim& xi,
                          const double       weight,
                          double&            result )
    {
        // Extract element pointer from tuple
        const GeomElement*  geomEp  = fieldTuple.geomElementPtr();

        // Get surface normal and measure
        GlobalVecDim normal;
        const double detJ = 
            base::SurfaceNormal<GeomElement>()( geomEp, xi, normal );


        for ( unsigned d = 0; d < globalDim; d++ )
            result += normal[d] * weight * detJ;
    }

};


#endif
