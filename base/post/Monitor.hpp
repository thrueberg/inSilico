//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Monitor.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_post_monitor_hpp
#define base_post_monitor_hpp

//------------------------------------------------------------------------------
// std includes
#include <string>
// base includes
#include <base/linearAlgebra.hpp>
// base/post includes
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace post{

        template<typename GEOMELEMENT, typename FIELDELEMENT>
        class Monitor;
    }
}


//------------------------------------------------------------------------------
/** Write out the solution value for a given elemeent and coordinate.
 *  For time-dependent problems, it is often interesting to 'watch' the solution
 *  or gradient at a specific location throughout the simulation process.
 *  This object provides a convenience interface by storing the pointers to the
 *  geometry and field elements and the value of the local evaluation
 *  coordinate. When called, it writes either the solution or the gradient to
 *  a provided stream object.
 *  \tparam GEOMELEMENT  Type of geometry element
 *  \tparam FIELDELEMENT Type of field element
 */
template<typename GEOMELEMENT, typename FIELDELEMENT>
class base::post::Monitor
{
public:
            
    //! @name Template parameter
    //@{
    typedef GEOMELEMENT  GeomElement;
    typedef FIELDELEMENT FieldElement;
    //@}

    //! Local coordinate type
    typedef typename FieldElement::FEFun::VecDim VecDim;

    //! Constructor with pointers and local evaluation coordinate
    Monitor( GeomElement*  geomElemPtr,
             FieldElement* fieldElemPtr,
             const VecDim  xi )
        : geomElemPtr_(  geomElemPtr ),
          fieldElemPtr_( fieldElemPtr ),
          xi_(           xi )
    { }

    //! Give geometry evaluation of the point of interest
    typename GeomElement::Node::VecDim location() const
    {
        return base::Geometry<GeomElement>()( geomElemPtr_, xi_ );
    }

    //! Write out the desired value
    void solution( std::ostream& out ) const
    {
        out << base::makeRow( base::post::evaluateField( geomElemPtr_,
                                                         fieldElemPtr_,
                                                         xi_ ) )
            << " ";
    }

    //! Write out the desired value
    void gradient( std::ostream& out ) const
    {
        out << base::makeRow( base::post::evaluateFieldGradient( geomElemPtr_,
                                                                 fieldElemPtr_,
                                                                 xi_ ) )
            << " ";
    }

private:
    const GeomElement*  geomElemPtr_;  //!< Pointer to geometry element
    const FieldElement* fieldElemPtr_; //!< Pointer to field element
    const VecDim        xi_;           //!< Local evaluation coordinate
};

#endif
