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
// base/cut includes
#include <base/cut/Marching.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<base::Shape SHAPE>
        class Cell;
    }
}

//------------------------------------------------------------------------------
/** Representation of the parametric geometry of a possibly cut cell.
 *  
 *  
 */
template<base::Shape SHAPE>
class base::cut::Cell
{
public:
    //! Template parameter: the shape of the cell
    static const base::Shape shape = SHAPE;

    //! @name Attributes
    //@{
    //! Local dimension 
    static const unsigned         dim = base::ShapeDim<shape>::value;
    //! Number of vertices of the shape
    static const unsigned numVertices = base::NumNFaces<shape,base::VERTEX>::value;
    //@}

    //! @name Convenience typedefs
    //@{
    typedef typename base::Vector<dim,double>::Type             VecDim;
    typedef typename base::Vector<dim-1,double>::Type           VecLDim;
    typedef typename base::cut::Simplex<dim-1,unsigned>::Type   SurfIndexSimplex;
    typedef typename base::cut::Simplex<dim,  unsigned>::Type   VolIndexSimplex;
    typedef typename base::cut::Simplex<dim-1,  VecDim>::Type   SurfSimplex;
    typedef typename base::cut::Simplex<dim,    VecDim>::Type   VolSimplex;
    //@}

    //! Default constructor initialises the tags
    Cell() // default state: in interior un-cut cell
        : isCut_( false ), isInside_( true ), isOutside_( false ) { }

    //--------------------------------------------------------------------------
    //! Construct with signed distance function values of the vertices
    void create( const boost::array<double,numVertices>& signedDistances )
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
            // create the internal structure of the cut cell
            base::cut::Marching<shape>::apply( signedDistances, nodes_, surface_,
                                               volumeIn_, volumeOut_ );
        }
        else {
            isCut_ = false;
            // decide location of cell which is not cut
            if ( distMax < 0. ) { isInside_ = false; isOutside_ = true;  }
            else                { isInside_ = true;  isOutside_ = false; }
        }

    }

    //! Clear all local data
    void destroy()
    {
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
    /** Perform a coordinate map (linear!).
     *
     */
    VecDim mapVolumeCoordinate( const VecDim& eta,
                                const unsigned s,
                                const bool inside ) const
    {
        const std::vector<VolIndexSimplex> & volume
            = ( inside ? volumeIn_ : volumeOut_ );

        VecDim xi0 = nodes_[ volume[s][0] ];
        VecDim xi  = xi0;
        for ( unsigned d = 0; d < dim; d++ )
            xi += ( nodes_[ volume[s][d+1] ] - xi0 ) * eta[d];
 
        return xi;
    }

    double volumeJacobian( const VecDim& eta,
                           const unsigned s,
                           const bool inside ) const
    {
        const std::vector<VolIndexSimplex> & volume
            = ( inside ? volumeIn_ : volumeOut_ );

        typename base::Matrix<dim,dim>::Type J;
        for ( unsigned d = 0; d < dim; d++ )
            J.col(d) = (nodes_[ volume[s][d+1] ] - nodes_[ volume[s][0] ] );

        return J.determinant();
    }

    double parameterArea( const bool inside ) const
    {
        const double refSize = base::RefSize<shape>::apply();

        if (     inside and isInside_  ) return refSize;
        if ( not inside and isOutside_ ) return refSize;

        double result = 0.;

        if ( isCut_ ) {
            const double simplexRefSize =
                base::RefSize<base::SimplexShape<dim>::value>::apply();
            
            const VecDim eta = base::constantVector<dim>( 0. );
            const std::size_t numVolSimp
                = ( inside ? volumeIn_.size() : volumeOut_.size() );
            
            for ( unsigned s = 0; s < numVolSimp; s++ )
                result += volumeJacobian( eta, s, inside ) * simplexRefSize;
        }

        return result;
    }

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
    
    double surfaceArea( ) const
    {
        double result = 0.;

        if ( isCut_ ) {
            const double simplexRefSize =
                base::RefSize<base::SimplexShape<dim-1>::value>::apply();
            
            const VecLDim eta = base::constantVector<dim-1>( 0. );
            
            for ( unsigned s = 0; s < surface_.size(); s++ )
                result += surfaceJacobian( eta, s ) * simplexRefSize;
        }

        return result;
    }

private:
    //------------------------------------------------------------------------------
    //! @name Location predicates
    //@{
    bool isCut_;     //!< Qualifier, if cell is cut
    bool isInside_;  //!< Qualifier, if cell is inside
    bool isOutside_; //!< Qualifier, if cell is outside
    //@}
    
    //------------------------------------------------------------------------------
    //! @name Structure
    //@{
    std::vector<VecDim>           nodes_;     //!< Parametric coordinates
    std::vector<SurfIndexSimplex> surface_;   //!< Connectivity of surface simplices
    std::vector<VolIndexSimplex>  volumeIn_;  //!< Connectivity of inside volume
    std::vector<VolIndexSimplex>  volumeOut_; //!< Connectivity of outside volume
    //@}
};


#endif
