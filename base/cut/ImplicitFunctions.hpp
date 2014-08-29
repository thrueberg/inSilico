//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ImplicitFunctions.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_cut_implicitfunctions_h
#define base_cut_implicitfunctions_h

#include <boost/function.hpp>

#include <base/linearAlgebra.hpp>
#include <base/auxi/compareNumbers.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{
        
        template<unsigned DIM> class Sphere;
        template<unsigned DIM> class Cylinder;
        template<unsigned DIM> class Plane;

        namespace detail_{
            template<unsigned DIM>
            struct ImplicitFun
            {
                typedef boost::function<bool( const typename base::Vector<DIM,double>::Type&,
                                              typename base::Vector<DIM,double>::Type& )>
                Type;
            };
        }
    }
}


//------------------------------------------------------------------------------
/** Representation of a DIM-Sphere by means of an implicit function.
 *  The sphere is specified by its radius and the coordinates of its centre.
 *  A function call operator sets the closest point on the sphere surface and
 *  returns if the given point is inside or outside.
 *  \tparam DIM Spatial dimension 
 */
template<unsigned DIM>
class base::cut::Sphere
    : public detail_::ImplicitFun<DIM>::Type
{
public:
    //! Template parameter: dimension of the sphere
    static const unsigned dim = DIM;

    //! Convenience typedef for coordinate
    typedef typename base::Vector<dim,double>::Type VecDim;

    //! Construct with radius and centre
    Sphere( const double radius,
            const VecDim centre,
            const bool inOut = true )
        : radius_( radius ),
          centre_( centre ),
          inOut_(  inOut  ) { }

    //! Set closest point on surface and return in/outside flag
    bool operator()( const VecDim& x, VecDim& xClosest ) const
    {
        // Distance vector to sphere centre
        const VecDim y = x - centre_;
        const double dist = y.norm();

        // Check against centre of sphere
        if ( base::auxi::almostEqualNumbers( dist, 0. ) ) {
            xClosest = base::constantVector<dim>(0.);
            xClosest[0] = radius_;
        }
        else{
            xClosest = (radius_ / dist) * y;
        }

        xClosest += centre_;

        if ( dist <= radius_ ) return inOut_;
        return not inOut_;
    }
    

private:
    const double radius_; //!< Radius of the sphere
    const VecDim centre_; //!< Centre of the sphere
    const bool   inOut_;  //!< Sign flag
};

//------------------------------------------------------------------------------
/** Representation of a DIM-Cylinder by means of an implicit function.
 *  The cylinder is determined by its radius, a point on its main axis
 *  and the direction of that axis. A function call operator sets the
 *  closest point on the sphere surface and returns if the given point
 *  is inside or outside.
 *  \tparam DIM Spatial dimension 
 */
template<unsigned DIM>
class base::cut::Cylinder
    : public detail_::ImplicitFun<DIM>::Type
{
public:
    //! Template parameter: dimension of the sphere
    static const unsigned dim = DIM;

    //! Convenience typedef for coordinate
    typedef typename base::Vector<dim,double>::Type VecDim;

    //! Set radius, axis point and direction
    Cylinder( const double radius,
              const VecDim centre,
              const VecDim direction,
              const bool inOut = true )
        : radius_(    radius ),
          centre_(    centre ),
          direction_( direction / direction.norm() ),
          inOut_(     inOut ) { }

    //! Set closest point on cylinder surface and return in/out flag
    bool operator()( const VecDim& x, VecDim& xClosest ) const
    {
        // distance vector between point and cylinder axis point
        const VecDim y  = x - centre_;
        // coordinate along the cylinder axis
        const double xi = y.dot( direction_ );
        // parallel component of distance vector
        const VecDim yp = xi * direction_;
        // orthogonal component of distance vector
        const VecDim yo = y - yp;
        // distance to cylinder axis
        const double dist = yo.norm();
        
        // check for points on axis
        if ( base::auxi::almostEqualNumbers( dist, 0. ) ) {
            xClosest = base::constantVector<dim>(0.);
            xClosest[0] = radius_;
        }
        else{
            xClosest = (radius_ / dist) * yo;
        }

        xClosest += centre_ + yp;

        if ( dist <= radius_ ) return inOut_;
        return not inOut_;
    }
    

private:
    const double radius_;     //!< Radius of the cylinder
    const VecDim centre_;     //!< Point on axis of cylinder
    const VecDim direction_;  //!< Unit vector in direction of axis
    const bool   inOut_;      //!< Sign flag
};


//------------------------------------------------------------------------------
/** Representation of a plane.
 *  \tparam DIM Spatial dimension 
 */
template<unsigned DIM>
class base::cut::Plane
    : public detail_::ImplicitFun<DIM>::Type
{
public:
    //! Template parameter: dimension of the sphere
    static const unsigned dim = DIM;

    //! Convenience typedef for coordinate
    typedef typename base::Vector<dim,double>::Type VecDim;

    //! Construct with radius and centre
    Plane( const VecDim normal,
           const double distance, 
           const bool inOut = true )
        : normal_(   normal ),
          distance_( distance ),
          inOut_(    inOut  ) { }

    //! Set closest point on surface and return in/outside flag
    bool operator()( const VecDim& x, VecDim& xClosest ) const
    {
        // Distance vector to sphere centre
        const double distToPlane = x.dot( normal_ ) - distance_;
        xClosest = x - distToPlane * normal_;

        if ( distToPlane <= 0. ) return inOut_;
        return not inOut_;
    }
    

private:
    const VecDim normal_;   //!< Centre of the sphere
    const double distance_; //!< Radius of the sphere
    const bool   inOut_;    //!< Sign flag
};

#endif
