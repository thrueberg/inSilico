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
 */
template<unsigned DIM>
class base::cut::LevelSet
{
public:
    static const unsigned dim = DIM;
    typedef typename base::Vector<dim,double>::Type VecDim;

    //--------------------------------------------------------------------------
    //! Constructor with coordinate, invalidates other member data
    LevelSet( const VecDim& x )
    : x_( x ),
      isInterior_( false ),
      closestPoint_( base::invalidVector<DIM>() ),
      closestElement_( base::invalidInt )
    { }

    //--------------------------------------------------------------------------
    //! @name Mutators
    //@{
    void setInterior() { isInterior_ = true; }
    void setExterior() { isInterior_ = false; }

    void setClosestPoint(   const VecDim&     y ) { closestPoint_ = y; }
    void setClosestElement( const std::size_t e ) { closestElement_ = e; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Accessors
    //@{
    VecDim getX()                const { return x_; }
    bool   isInterior()          const { return isInterior_; }
    VecDim getClosestPoint()     const { return closestPoint_; }
    double getUnsignedDistance() const { return base::norm(closestPoint_ - x_); }
    
    double getSignedDistances()  const
    {
        const double sign = ( isInterior_ ? +1.0 : -1.0 );
        return sign * (this -> getUnsignedDistance() );
    }
    
    std::size_t getClosestElement() const { return closestElement_; }
    //@}

    
private:
    //! Location of this datum
    const VecDim x_;

    //! @name Variable distance data
    //@{
    bool        isInterior_;
    VecDim      closestPoint_;
    std::size_t closestElement_;
    //@}
};

#endif
