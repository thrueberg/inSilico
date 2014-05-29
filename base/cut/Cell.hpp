//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Cell.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_cell_hpp
#define base_cut_cell_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <algorithm>
// base  includes
#include <base/shape.hpp>
#include <base/linearAlgebra.hpp>
#include <base/LagrangeShapeFun.hpp>
// base/cut includes
#include <base/cut/Marching.hpp>
#include <base/cut/MarchingUtils.hpp>
#include <base/cut/LevelSet.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<base::Shape SHAPE, unsigned DEGREE=1>
        class Cell;

        
    }
}

//------------------------------------------------------------------------------
/** Representation of the parametric geometry of a possibly cut cell.
 *  
 *  
 */
template<base::Shape SHAPE, unsigned DEGREE>
class base::cut::Cell
{
public:
    //! Template parameter: the shape of the cell
    static const base::Shape shape = SHAPE;
    //! Template parameter: Polynomial degree of the geometry representation
    static const unsigned geomDegree = DEGREE;

    //! @name Attributes
    //@{
    //! Local dimension 
    static const unsigned         dim = base::ShapeDim<shape>::value;
    //! Number of vertices of the shape
    static const unsigned numElementNodes =
        base::fe::LagrangeElement<SHAPE,DEGREE>::numTotalDoFs;
    //@}

    //! @name Simplex attributed deduced from geometry shape function
    //@{

    //! Volume
    static const base::Shape  volumeSimplex = base::SimplexShape<dim>::value;
    typedef base::LagrangeShapeFun<geomDegree, volumeSimplex>  VolumeFun;
    static const unsigned  nVolVert = VolumeFun::numFun;

    //! Surface
    static const base::Shape surfaceSimplex = base::SimplexShape<dim-1>::value;
    typedef base::LagrangeShapeFun<geomDegree,surfaceSimplex>  SurfaceFun;
    static const unsigned nSurfVert = SurfaceFun::numFun;

    //! Surface of surface (dummy for dim==2)
    //static const base::Shape surfSurfSimplex = base::SimplexShape<1>::value;
    //typedef base::LagrangeShapeFun<geomDegree,surfSurfSimplex> SurfSurfFun;
    //static const unsigned nSurfSurfVert = SurfSurfFun::numFun;
    //@}

    //! Linear lagrange function for vertices
    typedef base::LagrangeShapeFun<1,shape> LinearLagrange;


    //! @name Convenience typedefs
    //@{
    typedef typename base::Vector<dim,double>::Type             VecDim;
    typedef typename base::Vector<dim-1,double>::Type           VecLDim;
    
    typedef typename boost::array<unsigned,nSurfVert>     SurfIndexSimplex;
    typedef typename boost::array<unsigned,nVolVert >     VolIndexSimplex;
    //typedef typename boost::array<unsigned,nSurfSurfVert> SurfSurfIndexSimplex;
    typedef typename boost::array<VecDim,  nSurfVert>     SurfSimplex;
    typedef typename boost::array<VecDim,  nVolVert >     VolSimplex;
    //@}

    //! Default constructor initialises the tags
    Cell() // default state: in interior un-cut cell
        : volumeFun_( VolumeFun() ), surfaceFun_( SurfaceFun() ),
          isCut_( false ), isInside_( true ), isOutside_( false ) { }


    //! Necessary because shape functions cannot be copied
    base::cut::Cell<shape,geomDegree>&
    operator=(const base::cut::Cell<shape,geomDegree>& other )
    {
        isCut_     = other.isCut_;
        isInside_  = other.isInside_;
        isOutside_ = other.isOutside_;
        nodes_     = other.nodes_;
        surface_   = other.surface_;
        volumeIn_  = other.volumeIn_;
        volumeOut_ = other.volumeOut_;

        return *this;
    }
    
    //--------------------------------------------------------------------------
    //! Construct with signed distance function values of the vertices
    void create( const boost::array<double,numElementNodes>& signedDistances )
    {
        // Minimal value of signed distances
        const double distMin = *(std::min_element( signedDistances.begin(),
                                                   signedDistances.end() ) );
        // Maximal value of signed distances
        const double distMax = *(std::max_element( signedDistances.begin(),
                                                   signedDistances.end() ) );

        // Change of signs implies a cut cell
        if ( (distMin * distMax) <= 0. ) {
            isCut_     = true;
            isInside_  = false;
            isOutside_ = false;

            // get the vertices of the shape
            boost::array<VecDim,LinearLagrange::numFun> vertexCoordinates;
            LinearLagrange::supportPoints( vertexCoordinates );

            boost::array<unsigned,LinearLagrange::numFun> vertexIndices;
            
            // fill nodes with the vertices of this shape
            for ( unsigned v = 0; v < vertexCoordinates.size(); v++ ) {
                nodes_.push_back( vertexCoordinates[v] );
                vertexIndices[v] = v;
            }

            
            // create the internal structure of the cut cell
            std::map<base::cut::Edge,unsigned> uniqueNodes;
            base::cut::MarchingProxy<shape,geomDegree>::apply( signedDistances,
                                                               vertexIndices,
                                                               nodes_,
                                                               uniqueNodes,
                                                               surface_,
                                                               volumeIn_, volumeOut_,
                                                               false );

        }
        else {
            isCut_ = false;
            // decide location of cell which is not cut
            if ( distMax < 0. ) { isInside_ = false; isOutside_ = true;  }
            else                { isInside_ = true;  isOutside_ = false; }
        }

    }

    //--------------------------------------------------------------------------
    //! Construct with signed distance function values of the vertices
    void update( const boost::array<double,numElementNodes>& signedDistances )
    {

        if ( geomDegree > 1 ) {
            std::cerr << "(WW) The update of cut-cells does not work for "
                      << "higher-order geometry\n";

        }

        // Minimal value of signed distances
        const double distMin = *(std::min_element( signedDistances.begin(),
                                                   signedDistances.end() ) );
        // Maximal value of signed distances
        const double distMax = *(std::max_element( signedDistances.begin(),
                                                   signedDistances.end() ) );
        
        // Change of signs implies a cut cell
        if ( (distMin * distMax) <= 0. ) {

            // copy existing inside volume and surface to temporaries
            std::vector<VolIndexSimplex>  volumeInTmp = volumeIn_;
            std::vector<SurfIndexSimplex> surfaceTmp  = surface_;

            std::vector<VolIndexSimplex>  volumeOutTmp = volumeOut_;

            // delete the existing ones
            volumeIn_.clear();
            surface_.clear();

            volumeOut_.clear();

            // register already existing nodes
            std::map<base::cut::Edge,unsigned> uniqueNodes;

            // apply marching simplex to every volume simplex
            for ( std::size_t vs = 0; vs < volumeInTmp.size(); vs++ ) {
                // create an index vertex and an array of distance values
                boost::array<unsigned,dim+1>  volSimplex;
                boost::array<double,nVolVert> distances;
                
                for ( unsigned v = 0; v < nVolVert; v++ ) {
                    const unsigned vertexNum = volumeInTmp[vs][v];
                    
                    if ( v < dim+1 ) volSimplex[v] = vertexNum;
                    // interpolate distance to vertex node
                    distances[v]   =
                        detail_::InterpolatedDistance<shape,geomDegree>::evaluate(
                            nodes_[ vertexNum ], signedDistances );
                }

                base::cut::MarchingProxy<volumeSimplex,geomDegree>::apply(
                    distances,
                    volSimplex,
                    nodes_,
                    uniqueNodes,
                    surface_,
                    volumeIn_,
                    volumeOut_, true );
                
            }

#if 0 
            
            // apply marching simplex to every volume simplex
            for ( std::size_t vs = 0; vs < volumeOutTmp.size(); vs++ ) {
                // create an index vertex and an array of distance values
                boost::array<unsigned,dim+1>  volSimplex;
                boost::array<double,nVolVert> distances;
                
                for ( unsigned v = 0; v < nVolVert; v++ ) {
                    const unsigned vertexNum = volumeOutTmp[vs][v];
                    
                    if ( v < dim+1 ) volSimplex[v] = vertexNum;
                    // interpolate distance to vertex node
                    distances[v]   =
                        detail_::InterpolatedDistance<shape,geomDegree>::evaluate(
                            nodes_[ vertexNum ], signedDistances );
                }

                base::cut::MarchingProxy<volumeSimplex,geomDegree>::apply(
                    distances,
                    volSimplex,
                    nodes_,
                    uniqueNodes,
                    surface_,
                    volumeOut_,
                    volumeOut_, true );
                
            }

            // extract surfaces ...
            base::cut::GetInteriorSurface<shape,geomDegree>::apply(
                nodes_, surface_, volumeIn_, volumeOut_ );
#else

            // subtract cell boundary from surface
            static const unsigned surfSurfDim = dim > 1 ? dim - 2 : 0;
            static const base::Shape surfSurfSimplex = base::SimplexShape<surfSurfDim>::value;
            typedef base::LagrangeShapeFun<geomDegree,surfSurfSimplex> SurfSurfFun;
            static const unsigned nSurfSurfVert = SurfSurfFun::numFun;
            typedef typename boost::array<unsigned,nSurfSurfVert> SurfSurfIndexSimplex;
            
            // apply marching simplex to surface simplices
            for ( std::size_t ss = 0; ss < surfaceTmp.size(); ss++ ) {

                SurfIndexSimplex surfSimplex;
                boost::array<double,nSurfVert> distances;
                for ( unsigned v = 0; v < nSurfVert; v++ ) {
                    const unsigned vertexNum = surfaceTmp[ss][v];
                    surfSimplex[v] = vertexNum;
                    distances[v] =
                        detail_::InterpolatedDistance<shape,geomDegree>::evaluate(
                            nodes_[ vertexNum ], signedDistances );
                }
             
                std::vector<SurfIndexSimplex>     volOutDummy;
                std::vector<SurfSurfIndexSimplex> surfDummy;
                detail_::ApplyMarchingSimplex<dim,dim-1>::apply( surfSimplex,
                                                                 distances,
                                                                 nodes_,
                                                                 uniqueNodes,
                                                                 surfDummy,
                                                                 surface_,
                                                                 volOutDummy );
            }
            
#endif

        }
        else if ( distMax < 0. ) { // outside --> destroy old structure
            this -> destroy();
            isInside_  = false;
            isOutside_ = true;
            isCut_     = false;
            
        }

        return;
    }

    //! Clear all local data, set to default state
    void destroy()
    {
        isInside_  = true;
        isOutside_ = false;
        isCut_     = false;
        
        nodes_.clear();
        surface_.clear();
        volumeIn_.clear();
        volumeOut_.clear();
    }

    //-------------------------------------------------------------------------
    //! @name Predicates
    //@{
    //! Is this cell cut?
    bool isCut()     const { return isCut_; }
    //! Is this cell inside?
    bool isInside()  const { return isInside_; }
    //! Is this cell outside?
    bool isOutside() const { return isOutside_; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Query the sizes of the internal structure
    //@{
    std::size_t numNodes()             const { return nodes_.size(); }
    std::size_t numSurfaceElements()   const { return surface_.size(); }
    std::size_t numVolumeInElements()  const { return volumeIn_.size(); }
    std::size_t numVolumeOutElements() const { return volumeOut_.size(); }
    //@}

    //--------------------------------------------------------------------------
    //! @name Copy the internal structure
    //@{
    void getNodes( std::vector<VecDim>& nodes ) const
    {
        nodes = nodes_;
    }

    void getSurface( std::vector<SurfIndexSimplex> & surface ) const
    {
        surface = surface_;
    }

    void getVolumeIn( std::vector<VolIndexSimplex> & volumeIn ) const
    {
        volumeIn = volumeIn_;
    }

    void getVolumeOut( std::vector<VolIndexSimplex> & volumeOut ) const
    {
        volumeOut = volumeOut_;
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Perform coordinate maps (linear!)
    //@{

    /** Map from local simplex coordinate to cell-coordinate.
     *  The inside (outside) part of a cut cell composed of a collection of
     *  simplex elements. Given the number \f$ s \f$ of such a simplex and the
     *  simplex coordinate \f$ \eta \f$, this function returns the coordinate
     *  \f$ \xi \f$ within the cell
     *  \f[
     *      \xi^s = \sum \phi_J(\eta) \xi_J^s
     *  \f]
     *  based on the local nodes \f$ \xi_J^s \f$ of simplex \f$ s \f$ and the
     *  simple xshape functions \f$ \phi_J \f$.
     *  \param[in] eta    Simplex coordinate
     *  \param[in] s      Number of the local simplex
     *  \param[in] inside Flag if in- or outside volume is of interest
     *  \return           The cell-coordinate \f$ \xi \f$
     */
    VecDim mapVolumeCoordinate( const VecDim& eta,
                                const unsigned s,
                                const bool inside ) const
    {
        const std::vector<VolIndexSimplex> & volume
            = ( inside ? volumeIn_ : volumeOut_ );

        typename VolumeFun::FunArray phi;
        volumeFun_.fun( eta, phi );

        VecDim xi = base::constantVector<dim>( 0. );
        for ( unsigned f = 0; f < VolumeFun::numFun; f++ )
            xi += phi[f] * nodes_[ volume[s][f] ];

        return xi;
    }

    /** Compute the Jacobian of the map given in mapVolumeCoordinate
     *  \param[in] eta    Simplex coordinate
     *  \param[in] s      Number of the local simplex
     *  \param[in] inside Flag if in- or outside volume is of interest
     *  \return           Jacobian of coordinate map \f$ \xi(\eta) \f$
     */
    double volumeJacobian( const VecDim& eta,
                           const unsigned s,
                           const bool inside ) const
    {
        const std::vector<VolIndexSimplex> & volume
            = ( inside ? volumeIn_ : volumeOut_ );

        typename VolumeFun::GradArray dPhi;
        volumeFun_.gradient( eta, dPhi );

        typename base::Matrix<dim,dim>::Type J = base::constantMatrix<dim,dim>( 0. );
        for ( unsigned f = 0; f < VolumeFun::numFun; f++ )
            J += nodes_[ volume[s][f] ] * ( dPhi[f] ).transpose();

        return J.determinant();
    }

    /** Compute the Jacobian of the embedded surface simplex elements.
     *  Here, the map
     *  \f[
     *        \xi^s = \sum_K \psi_K(\eta) \xi^s_K
     *  \f]
     *  assigns a dim-1 coordinate \f$ \eta \f$ a cell coordinate
     *  \f$ \xi \f$ on the embedded surface for simplex \f$ s \f$.
     *  This function computes the Jacobian of this map.
     *  \param[in] eta  Local surface coordiante
     *  \param[in] s    Number of the surface simplex
     *  \return         Jacobian of the embedding.
     */
    double surfaceJacobian( const VecLDim& eta,
                            const unsigned s ) const
    {
        typename base::Matrix<dim,dim-1>::Type J;
        for ( unsigned d = 0; d < dim-1; d++ )
            J.col(d) = (nodes_[ surface_[s][d+1] ] -
                        nodes_[ surface_[s][0  ] ] );

        return
            std::sqrt( (J.transpose() * J).determinant() );
    }

    /** Return the approximate volume in- or outside of the domain.
     *  Calling volumeJacobian for all simplices inside (outside) at
     *  the simplex centroid and summing the result gives a measure of the
     *  paramtric volume inside (outside).
     *  \param[in] inside  in/outside flag
     *  \return            Approximate volume measure
     */
    double parameterVolume( const bool inside ) const
    {
        // size of the reference shape
        const double refSize = base::RefSize<shape>::apply();

        if (     inside and isInside_  ) return refSize;
        if ( not inside and isOutside_ ) return refSize;

        double result = 0.;

        // ask all simplices
        if ( isCut_ ) {
            const double simplexRefSize =
                base::RefSize<base::SimplexShape<dim>::value>::apply();
            
            const VecDim eta = base::detail_::SimplexCentroid<dim>::apply();
            const std::size_t numVolSimp
                = ( inside ? volumeIn_.size() : volumeOut_.size() );
            
            for ( unsigned s = 0; s < numVolSimp; s++ )
                result += volumeJacobian( eta, s, inside ) * simplexRefSize;
        }

        return result;
    }

    /** Return the approximate area of the embedded surface.
     *  Going over all surface simplex elements, call the surfaceJacobian
     *  function at the simplex centroid and sum up the result.
     *  \return Approximate embedded surface area
     */
    double surfaceArea( ) const
    {
        double result = 0.;

        if ( isCut_ ) {
            const double simplexRefSize =
                base::RefSize<base::SimplexShape<dim-1>::value>::apply();
            
            const VecLDim eta = base::detail_::SimplexCentroid<dim-1>::apply();
            
            for ( unsigned s = 0; s < surface_.size(); s++ )
                result += surfaceJacobian( eta, s ) * simplexRefSize;
        }

        return result;
    }
    //@}

    //--------------------------------------------------------------------------
    /** Discard the simplex elements with measure below a threshold.
     *  Due to the repeated cutting process, sometimes simplex elements appear
     *  which have a measure (volume or surface area size) that is too small
     *  for the following computations. This function deletes such elements from
     *  the lists.
     *  \param[in]  tol  The tolerance for this decision
     */
    void compress( const double tol )
    {
        // memorise the indices of the simplex elements to be deleted
        std::vector<std::size_t> volInSkip, volOutSkip, surfSkip;

        
        // check volume-in
        for ( std::size_t s = 0; s < volumeIn_.size(); s++ ) {
            typename base::cut::VecSimplex<dim,dim>::Type simplex;
            for ( unsigned d = 0; d < dim+1; d++ )
                simplex[d] = nodes_[ volumeIn_[s][d] ];
            
            const double orient =
                base::cut::simplexOrient<dim>( simplex );

            if ( orient < tol ) volInSkip.push_back( s );
        }

        // check volume-out
        for ( std::size_t s = 0; s < volumeOut_.size(); s++ ) {
            typename base::cut::VecSimplex<dim,dim>::Type simplex;
            for ( unsigned d = 0; d < dim+1; d++ )
                simplex[d] = nodes_[ volumeOut_[s][d] ];
            
            const double orient =
                base::cut::simplexOrient<dim>( simplex );

            if ( orient < tol ) volOutSkip.push_back( s );
        }

        // check surface
        for ( std::size_t s = 0; s < surface_.size(); s++ ) {
            typename base::cut::VecSimplex<dim,dim-1>::Type simplex;
            for ( unsigned d = 0; d < dim; d++ )
                simplex[d] = nodes_[ surface_[s][d] ];
            
            const double size =
                base::cut::surfSimplexSize<dim>( simplex );

            if ( size < tol ) surfSkip.push_back( s );
        }

        // Make sure to start deleting process from end of containers
        std::reverse( volInSkip.begin(),  volInSkip.end() );
        std::reverse( volOutSkip.begin(), volOutSkip.end() );
        std::reverse( surfSkip.begin(),   surfSkip.end() );

        // erase simplices, starting with the highest numbers
        for ( std::size_t s = 0; s < volInSkip.size();  s++ )
            volumeIn_.erase(  volumeIn_.begin() + volInSkip[s] );
        for ( std::size_t s = 0; s < volOutSkip.size(); s++ )
            volumeOut_.erase( volumeOut_.begin() + volOutSkip[s] );
        for ( std::size_t s = 0; s < surfSkip.size();   s++ )
            surface_.erase(   surface_.begin() + surfSkip[s] );

    }
    
private:
    //--------------------------------------------------------------------------
    //! @name Shape functions of the cut cell simplex structure
    //@{
    VolumeFun    volumeFun_;
    SurfaceFun   surfaceFun_;
    //@}
    
    //--------------------------------------------------------------------------
    //! @name Location predicates
    //@{
    bool isCut_;     //!< Qualifier, if cell is cut
    bool isInside_;  //!< Qualifier, if cell is inside
    bool isOutside_; //!< Qualifier, if cell is outside
    //@}
    
    //--------------------------------------------------------------------------
    //! @name Structure
    //@{
    std::vector<VecDim>           nodes_;     //!< Parametric coordinates
    std::vector<SurfIndexSimplex> surface_;   //!< Connectivity of surface simplices
    std::vector<VolIndexSimplex>  volumeIn_;  //!< Connectivity of inside volume
    std::vector<VolIndexSimplex>  volumeOut_; //!< Connectivity of outside volume
    //@}
};


#endif
