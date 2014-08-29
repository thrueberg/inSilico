//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   AdjointBodyForce.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef tfm_adjointbodyforce_hpp
#define tfm_adjointbodyforce_hpp

//------------------------------------------------------------------------------
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace tfm{

    template<typename MESH>
    struct AlwaysPositive
        : public boost::function<bool( const typename MESH::Node::VecDim& ) >
    {
        bool operator()( const typename MESH::Node::VecDim& x )
        {
            return true;
        }
    };
    
    //--------------------------------------------------------------------------
    template<typename MESH, typename FIELD>
    class AdjointBodyForce
    {
    public:

        typedef typename base::Vector<FIELD::DegreeOfFreedom::size>::Type VecDoF;
        typedef typename base::Vector<MESH::Node::dim>::Type              VecDim;
        typedef boost::function<bool( const VecDim& )> FilterFun;
        
        AdjointBodyForce( MESH&  mesh,
                          FIELD& currentU,
                          FIELD& targetU,
                          FilterFun filterFun = boost::bind( tfm::AlwaysPositive<MESH>(), _1 ) )
            : mesh_( mesh ), 
              currentU_( currentU ),
              targetU_( targetU ),
              filterFun_( filterFun )
        { }


        VecDoF apply( const typename MESH::Element* ep,
                      const typename MESH::Element::GeomFun::VecDim& xi ) const
        {
            typename MESH::Node::VecDim x =
                base::Geometry<typename MESH::Element>()( ep, xi );

            if ( not filterFun_( x ) ) {
                return base::constantVector<FIELD::DegreeOfFreedom::size>( 0. );
            }
            
            const std::size_t elemID = ep -> getID();

            typename FIELD::Element* cuEPtr = currentU_.elementPtr( elemID );
            typename FIELD::Element* tuEptr = targetU_.elementPtr(  elemID );


            return
                base::post::evaluateField( ep, cuEPtr, xi ) -
                base::post::evaluateField( ep, tuEptr, xi );
        }
        
    private:
        MESH&  mesh_;
        FIELD& currentU_;
        FIELD& targetU_;
        FilterFun filterFun_;
    };
}

//------------------------------------------------------------------------------
#endif
