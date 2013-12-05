//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Node.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_mesh_node_hpp
#define base_mesh_node_hpp

//------------------------------------------------------------------------------
// boost includes
#include <boost/utility.hpp>
// base  includes
#include <base/linearAlgebra.hpp>
#include <base/numbers.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace mesh{
        template<unsigned DIM>
        class Node;
    }
}

//------------------------------------------------------------------------------
/**
 */
template<unsigned DIM>
class base::mesh::Node
    : public boost::noncopyable
{
public:
    //! Template parameter: spatial dimension
    static const unsigned dim = DIM;

    //! Coordinate type definition
    typedef typename base::Vector<dim>::Type VecDim;

    //--------------------------------------------------------------------------
    //! @name Constructors
    //@{

    /** Empty constructor invalidates the member data
     */
    Node()
        : id_( base::invalidInt ),
          x_(  base::invalidVector<dim>() )
    { }

    /** Basic constructor with ID and Coordinate
     *  \param[in] id  Global ID of this node
     *  \param[in] x   Coordinate of this node
     */
    Node( const std::size_t id,
          const VecDim & x )
        : id_( id ), x_( x )
    { }
    //@}

    //--------------------------------------------------------------------------
    //! @name Mutators
    //@{
    //! Set the global ID
    void setID( const std::size_t id ) { id_ = id; }
    //! Set the coordinate
    void setX(  const VecDim & x  ) { x_  = x;  }
    //! Set coordinates from iterator
    template<typename INPITER>
    void setX( INPITER iter )
    {
        for ( unsigned d = 0; d < dim; d ++ ) x_[d] = *iter++;
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Accessors
    //@{
    //! Return global ID of node
    std::size_t getID() const { return id_; }
    //! Return coordinate of node
    VecDim   getX()  const { return x_;  }
    //! Pass coordinates to an iterator
    template<typename OUTITER>
    void getX( OUTITER iter ) const
    {
        for ( unsigned d = 0; d < dim; d ++ ) *iter++ = x_[d];
    }
    //@}
    
private:
    std::size_t id_; //!< Global Node ID
    VecDim      x_;  //!< Coordinate of this node
};
//------------------------------------------------------------------------------
#endif
