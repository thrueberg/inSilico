//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   base/cut/Quadrature.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_cut_quadrature_hpp
#define base_cut_quadrature_hpp

//------------------------------------------------------------------------------
// base  includes
#include <base/shape.hpp>
#include <base/meta.hpp>
#include <base/Quadrature.hpp>
// base/cut includes
#include <base/cut/Cell.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<unsigned DEGREE, base::Shape SHAPE,
                 typename CELL=base::cut::Cell<SHAPE> >
        class Quadrature;

        //----------------------------------------------------------------------
        namespace detail_{

            template<typename ELEMENT>
            struct AskVolumeID
            {
                static std::size_t apply( const ELEMENT* ep )
                {
                    return ep -> getID();
                }
            };

            template<typename ELEMENT>
            struct AskSurfaceID
            {
                static std::size_t apply( const ELEMENT* ep )
                {
                    return ep -> getSurfaceID();
                }
            };

            //! Hack
            template<typename ELEMENT>
            struct GetID
                : public base::IfElse<ELEMENT::dim==ELEMENT::Node::dim,
                                      AskVolumeID<ELEMENT>,
                                      AskSurfaceID<ELEMENT> >::Type
            { };

            template<typename ELEMENT>
            std::size_t getID( const ELEMENT* ep )
            {
                return GetID<ELEMENT>::apply( ep );
            }

        }
        
    }
}

//------------------------------------------------------------------------------
/** Quadrature rule for possibly cut elements.
 *  In case of immersed finite element methods, integration of subregions of
 *  an element have to be carried out. Such a subregion has been triangulated
 *  and the quadrature thus is a composite rule over the simplex elements the
 *  form the subregion. The simplex structure is stored in the \e CELL which
 *  corresponds to the considered element. If the considered element is full
 *  inside the domain a standard quadrature rule for that element will be
 *  chosen.
 *  \tparam DEGREE Polynomial degree to be integrated exactly by the rule
 *  \tparam SHAPE  Shape of the element
 *  \tparam CELL   Type of cell representing the cut-element structure
 */
template<unsigned DEGREE, base::Shape SHAPE, typename CELL>
class base::cut::Quadrature
{
public:
    //!@ Template parameter
    //@{
    static const unsigned    degree = DEGREE;
    static const base::Shape shape  = SHAPE;
    typedef      CELL                 Cell;
    //@}

    //! Shape of a volume simplex
    static const base::Shape simplexShape =
        base::SimplexShape<base::ShapeDim<shape>::value>::value;

    //!@ Quadratures to be used
    //@{
    typedef base::Quadrature<degree,shape>        StandardQuad;
    typedef base::Quadrature<degree,simplexShape> SimplexQuad;
    //@}

    //! for introspection
    typedef typename base::Vector<base::ShapeDim<shape>::value>::Type VecDim;

    //! Constructor with access to cut cells and flag for in/outside 
    Quadrature( std::vector<Cell>& cells,
                const bool inside = true )
        : cells_( cells ), inside_( inside )
    {
    }

    //--------------------------------------------------------------------------
    //! Main function to perform integration
    template<typename KERNEL>
    void apply( KERNEL& kernel,
                typename KERNEL::arg1_type& arg1,
                typename KERNEL::arg4_type& arg4 ) const
    {
        // get the ID of the element
        //const std::size_t elemID = arg1.geomElementPtr() -> getID();
        const std::size_t elemID = detail_::getID( arg1.geomElementPtr() );

        // For cut-elements go through stored simplex structure
        if ( cells_[elemID].isCut() ) {

            // number of simplices forming the cut-cell volume 
            const std::size_t numVolSimplices =
                ( inside_ ?
                  cells_[elemID].numVolumeInElements() :
                  cells_[elemID].numVolumeOutElements() );

            // go through simplices forming the volume of the cut element
            for ( unsigned s = 0; s < numVolSimplices; s++ ) {

                // go through the points of the simplex quadrature rule
                typename SimplexQuad::Iter sIter = simplexQuad_.begin();
                typename SimplexQuad::Iter sEnd  = simplexQuad_.end();
                for ( ; sIter != sEnd; ++sIter ) {

                    // quadrature point in reference element
                    const typename SimplexQuad::VecDim eta
                        = sIter -> second;

                    // apply map to sub-simplex of cut-element
                    const typename SimplexQuad::VecDim xi
                        = cells_[elemID].mapVolumeCoordinate( eta, s, inside_ );

                    // get modified weight for sub-simplex integration
                    const double weight = (sIter -> first)
                        * cells_[elemID].volumeJacobian( eta, s, inside_ );

                    // apply the quadrature rule
                    kernel( arg1, xi, weight, arg4 );

                }
            } // loop over simplices

        }
        else if ( ( cells_[elemID].isInside()  and     inside_ ) or
                  ( cells_[elemID].isOutside() and not inside_ ) ) {

            // apply standard quadrature
            standardQuad_.apply( kernel, arg1, arg4 );
        }

        return;
    }
        
    //! Change between inside or outside integration
    void flipInside() { inside_ = not inside_; }

private:
    const std::vector<Cell>&  cells_;  //!< Access to cut-cell structures
    bool                      inside_; //!< True for integration inside only

    StandardQuad  standardQuad_;       //!< Quadrature for non-cut cells
    SimplexQuad   simplexQuad_;        //!< Quadrature for cut cells
};

#endif
