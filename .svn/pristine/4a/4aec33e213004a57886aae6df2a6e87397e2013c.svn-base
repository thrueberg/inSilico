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
#include <base/Quadrature.hpp>
// base/cut includes
#include <base/cut/Cell.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace cut{

        template<unsigned DEGREE, base::Shape SHAPE,
                 typename CELL=base::cut::Cell<SHAPE> >
        class Quadrature;

    }
}

//------------------------------------------------------------------------------
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
    Quadrature( const std::vector<Cell>& cells,
                const bool inside = true )
        : cells_( cells ), inside_( inside )
    {
    }

    //! Main function to perform integration
    template<typename KERNEL>
    void apply( KERNEL& kernel,
                typename KERNEL::arg1_type& arg1,
                typename KERNEL::arg4_type& arg4 ) const
    {
        // get the ID of the element
        const std::size_t elemID = arg1.geomElementPtr() -> getID();

        // choose upon status
        if ( cells_[elemID].isCut() ) {

            // number of simplices forming the cut-cell volume 
            const std::size_t numVolSimplices =
                ( inside_ ?
                  cells_[elemID].numVolumeInElements() :
                  cells_[elemID].numVolumeOutElements() );

            for ( unsigned s = 0; s < numVolSimplices; s++ ) {

                typename SimplexQuad::Iter sIter = simplexQuad_.begin();
                typename SimplexQuad::Iter sEnd  = simplexQuad_.end();
                for ( ; sIter != sEnd; ++sIter ) {

                    const typename SimplexQuad::VecDim eta
                        = sIter -> second;

                    const double weight = (sIter -> first)
                        * cells_[elemID].volumeJacobian( eta, s, inside_ );

                    const typename SimplexQuad::VecDim xi
                        = cells_[elemID].mapVolumeCoordinate( eta, s, inside_ );

                    kernel( arg1, xi, weight, arg4 );

                }
            }

        }
        else if ( ( cells_[elemID].isInside()  and     inside_ ) or
                  ( cells_[elemID].isOutside() and not inside_ ) ) {
        
            typename StandardQuad::Iter qIter = standardQuad_.begin();
            typename StandardQuad::Iter qEnd  = standardQuad_.end();
            for ( ; qIter != qEnd; ++qIter )
                kernel( arg1, qIter -> second, qIter -> first, arg4 );

        }

    }


    void flipInside() { inside_ = not inside_; }

private:
    const std::vector<Cell>&  cells_;  //!< Access to cut-cell structures
    bool                      inside_; //!< True for integration inside only

    StandardQuad  standardQuad_;       //!< Quadrature for non-cut cells
    SimplexQuad   simplexQuad_;        //!< Quadrature for cut cells
};

#endif
