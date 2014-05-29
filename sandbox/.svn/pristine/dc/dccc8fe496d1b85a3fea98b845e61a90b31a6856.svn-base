#ifndef surface_h
#define surface_h

#include <base/asmb/SimpleIntegrator.hpp>
#include <base/kernel/FieldIntegral.hpp> 

//------------------------------------------------------------------------------
//! Analytic level set function for a spherical domain
template<unsigned DIM>
bool spherical( const typename base::Vector<DIM>::Type& x,
                typename base::Vector<DIM>::Type& xClosest,
                const typename base::Vector<DIM>::Type& c,
                const double R )
{
    const typename base::Vector<DIM>::Type y = x - c;
    
    if ( y.norm() < coordTol ) {
        xClosest =    c;
        xClosest[0] = R;
    }
    else{
        xClosest = (R / y.norm()) * y + c;
    }

    if ( y.norm() <= R ) return true;
    return false;
}

//------------------------------------------------------------------------------
// Compute the enclosed volume, the centroid and principal directions
template<typename SUU, typename SURFFIELDBIND, typename SQUAD, typename SCQUAD>
void surfaceFeatures( const SURFFIELDBIND& surfaceFieldBinder,
                      const SURFFIELDBIND& boundaryFieldBinder,
                      const SQUAD&  surfaceQuadrature,
                      const SCQUAD& surfaceCutQuadrature,
                      std::ostream& out )
{
    static const unsigned dim = SURFFIELDBIND::Mesh::Node::dim;
    typedef typename base::Vector<dim>::Type     VecDim;
    typedef typename base::Matrix<dim,dim>::Type MatDimDim;
 
    // get shape features
    const double enclosed =
        surf::enclosedVolume<SUU>( surfaceQuadrature, surfaceFieldBinder ) +
        surf::enclosedVolume<SUU>( surfaceCutQuadrature, boundaryFieldBinder );

    const VecDim moment =
        surf::enclosedVolumeMoment<SUU>( surfaceQuadrature, surfaceFieldBinder ) +
        surf::enclosedVolumeMoment<SUU>( surfaceCutQuadrature, boundaryFieldBinder );

    const MatDimDim secondMoment =
        surf::enclosedVolumeSecondMoment<SUU>( surfaceQuadrature, surfaceFieldBinder ) +
        surf::enclosedVolumeSecondMoment<SUU>( surfaceCutQuadrature, boundaryFieldBinder );

    // process
    const VecDim centroid = moment / enclosed;
     
    MatDimDim inertia;
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            inertia(d1, d2) = secondMoment(d1,d2) - enclosed * centroid[d1] * centroid[d2];
        }
    }
         
    // compute principal values of inertia and the angle-axis rotation
    MatDimDim evec;
    const VecDim princVal = base::eigenPairs( inertia, evec );
    const std::pair<double,base::Vector<3>::Type> aa = base::angleAxis( evec );
     
    out << enclosed << "  " << centroid.transpose() << "  "
        << princVal.transpose() << "  "
        << aa.first << "  "  << (aa.second).transpose();
}

//------------------------------------------------------------------------------
// Compute the enclosed volume, the centroid and principal directions
template<typename FFI, typename FLUID, typename SQUAD, typename QUAD>
void bubbleFeatures( FFI& ffi, FLUID& fluid,
                     const SQUAD& surfaceQuadrature,
                     const QUAD&  quadrature, 
                     std::ostream& out )
{
    static const unsigned dim = FLUID::dim;
    typedef typename base::Vector<dim>::Type     VecDim;
    typedef typename base::Matrix<dim,dim>::Type MatDimDim;
 
    // get shape features
    const double enclosed =
        surf::enclosedVolume<typename FFI::SU1U1>( surfaceQuadrature, ffi.getBinder() );

    const VecDim moment =
        surf::enclosedVolumeMoment<typename FFI::SU1U1>( surfaceQuadrature, ffi.getBinder() );

    const MatDimDim secondMoment =
        surf::enclosedVolumeSecondMoment<typename FFI::SU1U1>( surfaceQuadrature,
                                                               ffi.getBinder() );

    // process
    const VecDim centroid = moment / enclosed;
     
    MatDimDim inertia;
    for ( unsigned d1 = 0; d1 < dim; d1++ ) {
        for ( unsigned d2 = 0; d2 < dim; d2++ ) {
            inertia(d1, d2) = secondMoment(d1,d2) - enclosed * centroid[d1] * centroid[d2];
        }
    }
         
    // compute principal values of inertia and the angle-axis rotation
    MatDimDim evec;
    const VecDim princVal = base::eigenPairs( inertia, evec );
    const std::pair<double,base::Vector<3>::Type> aa = base::angleAxis( evec );

    // average rise
    VecDim averageU = base::constantVector<dim>( 0. );
    base::asmb::simplyIntegrate<typename FLUID::UU>(
        quadrature, averageU,
        typename FLUID::FieldBinder( fluid.getMesh(), fluid.getVelocity(), fluid.getPressure() ),
        base::kernel::FieldIntegral<typename FLUID::UU::Tuple>() );
    averageU /= enclosed;
    
     
    out << enclosed << "  "
        << centroid.transpose() << "  "
        << averageU.transpose() << "  "
        << princVal.transpose() << "  "
        << aa.first << "  "  << (aa.second).transpose();
}

#endif
