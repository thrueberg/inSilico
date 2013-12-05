//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   LevelSet.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_levelset_hpp
#define base_cut_levelset_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<unsigned DIM>
        class LevelSet;
        
    }
}

//------------------------------------------------------------------------------
/** Representation of the level set datum at a given point.
 *  The level set method gives an implicit surface representation as the zero
 *  contour level of a function, most commonly a (signed) distance function.
 *  Moreover, boundary data on the original surface needs to be accessed.
 *  Therefore, the surface point \f$ x_i^* \f$ closest to the mesh (or grid)
 *  point \f$ x_i \f$ is needed and the index \f$ \tau_i^* \f$ of the surface
 *  element it lies on. Denoting by \f$ d_i = d(x_i) \f$ the (signed) distance
 *  of the point \f$ x_i \f$ of question, i.e. \f$ d_i = | x_i - x_i^* | \f$,
 *  this object holds all data of the map
 *  \f[
 *        x_i  \to \{  sign(d_i), x_i^*, \tau_i^* \}
 *  \f]
 *  \tparam DIM Spatial dimension of the problem
 */
template<unsigned DIM>
class base::cut::LevelSet
{
public:
    static const unsigned dim = DIM;
    typedef typename base::Vector<dim,double>::Type   VecDim;
    typedef typename base::Vector<dim-1,double>::Type LocalVecDim;

    //--------------------------------------------------------------------------
    //! Constructor with coordinate, invalidates other member data
    LevelSet( const VecDim& x = base::invalidVector<dim>() )
    : x_( x ),
      isInterior_( false ),
      closestPoint_( base::invalidVector<dim>() ),
      closestElement_( base::invalidInt ),
      localCoordinate_( base::invalidVector<dim-1>() )
    { }

    //--------------------------------------------------------------------------
    //! @name Mutators
    //@{
    void setInterior() { isInterior_ = true; }
    void setExterior() { isInterior_ = false; }

    void setClosestPoint(   const VecDim&     y ) { closestPoint_   = y; }
    void setClosestElement( const std::size_t e ) { closestElement_ = e; }
    void setClosestLocalCoordinate( const LocalVecDim& xi ) { localCoordinate_ = xi; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Accessors
    //@{
    VecDim getX()                const { return x_; }
    bool   isInterior()          const { return isInterior_; }
    VecDim getClosestPoint()     const { return closestPoint_; }
    double getUnsignedDistance() const { return base::norm(closestPoint_ - x_); }
    
    double getSignedDistance()  const
    {
        const double sign = ( isInterior_ ? +1.0 : -1.0 );
        return sign * (this -> getUnsignedDistance() );
    }
    
    std::size_t getClosestElement() const { return closestElement_; }

    LocalVecDim getClosestLocalCoordinate() const { return localCoordinate_; }
    //@}

    
private:
    //! Location of this datum
    VecDim x_;

    //! @name Variable distance data
    //@{
    bool        isInterior_;     //!< Flag, if the point x_ is interior
    VecDim      closestPoint_;   //!< Coordinates of the closest surface point 
    std::size_t closestElement_; //!< Index of the closest surface element
    LocalVecDim localCoordinate_;//!< Element coordinate of closest point
    //@}
};

#endif
