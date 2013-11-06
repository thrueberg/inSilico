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

//------------------------------------------------------------------------------
namespace surf{

    template<typename FIELDTUPLE>
    class Moments;

    //--------------------------------------------------------------------------
    /** Compute the volume enclosed by surface mesh.
     *  \tparam FTB Field-tuple binder (actually unsused as only geometry needed)
     *  \tparam QUADRATURE Surface quadrature
     *  \tparam FIELDBINDER Binder with surface mesh in front
     */
    template<typename FTB, typename QUADRATURE, typename FIELDBINDER>
    double enclosedVolume( const QUADRATURE&  quadrature,
                           const FIELDBINDER& fieldBinder )
    {
        // sum up the volume
        double result = 0.;

        // perform integration over surface
        base::asmb::simplyIntegrate<FTB>(
            quadrature, result, fieldBinder,
            boost::bind( &surf::Moments<typename FTB::Tuple>::enclosedVolume,
                         _1, _2, _3, _4 ) );
        
        return result;
    }

    //--------------------------------------------------------------------------
    /** Compute the volume moment as enclosed by a surface mesh.
     *  \tparam FTB Field-tuple binder (actually unsused as only geometry needed)
     *  \tparam QUADRATURE Surface quadrature
     *  \tparam FIELDBINDER Binder with surface mesh in front
     */
    template<typename FTB, typename QUADRATURE, typename FIELDBINDER>
    typename base::Vector<FIELDBINDER::Mesh::Node::dim>::Type
    enclosedVolumeMoment( const QUADRATURE&  quadrature,
                          const FIELDBINDER& fieldBinder )
    {
        static const unsigned dim = FIELDBINDER::Mesh::Node::dim;
        // 
        typename base::Vector<dim>::Type moment = base::constantVector<dim>(0.);
        
        base::asmb::simplyIntegrate<FTB>(
            quadrature, moment, fieldBinder,
            boost::bind( &surf::Moments<typename FTB::Tuple>::areaMoment,
                         _1, _2, _3, _4 ) );

        return moment;
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

    //! @name Local and global coordinates
    //@{
    typedef typename base::GeomTraits<GeomElement>::LocalVecDim   LocalVecDim;
    typedef typename base::GeomTraits<GeomElement>::GlobalVecDim  GlobalVecDim;
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
     *  The first moment of area is defined as
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
    static void areaMoment( const FieldTuple&  fieldTuple,
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

};


#endif
