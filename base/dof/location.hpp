//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   location.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_location_hpp
#define base_dof_location_hpp

//------------------------------------------------------------------------------
// std includes
#include <vector>
#include <utility>
#include <iostream>
// boost array
#include <boost/array.hpp>
// base includes
#include <base/numbers.hpp>
#include <base/geometry.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename FIELD>
        void associateLocation( const FIELD& field,
                                std::vector<
                                std::pair<std::size_t,
                                typename FIELD::Element::FEFun::VecDim> >&
                                location );

        template<typename MESH>
        std::size_t findDoFWithLocation( const std::vector<
                                         std::pair<std::size_t,
                                         typename MESH::Element::GeomFun::VecDim> >&
                                         location,
                                         const MESH& mesh, 
                                         const typename MESH::Node::VecDim& x,
                                         const double tolerance );
    }
}

//------------------------------------------------------------------------------
/** Associate a location with every Degree of Freedom.
 *   By design, inSilico liberates the method from the standard approach of
 *   binding the field to the geometry (aka Nodal FEM). Although the chosen way
 *   is more flexible, it is sometimes necessary to have an association between
 *   the degrees of freedom and a geometric location. For instance, in case of
 *   immersed boundaries, the location of a degree of freedom with respect to
 *   that surface is of interest.  This function generates a pair of element
 *   for every degree of freedom of the given field by evluation of the
 *   geometry (via the provided mesh) at the support points associated with the
 *   inidividual degrees of freedom.
 *   \tparam FIELD Type of FE field
 *   \param[in]  field       The field for whose DoFs the location is sought
 *   \param[out] location    Outcome of this method
 */
template<typename FIELD>
void base::dof::associateLocation( const FIELD& field,
                                   std::vector< std::pair<std::size_t,
                                   typename FIELD::Element::FEFun::VecDim> >&
                                   location )
{
    typedef typename FIELD::Element::FEFun::VecDim LocalVecDim;

    // support points used for geometry association
    boost::array<LocalVecDim,FIELD::Element::FEFun::numFun> supportPoints;
    FIELD::Element::FEFun::supportPoints( supportPoints );

    // storage
    const std::size_t numDoFs = std::distance( field.doFsBegin(), field.doFsEnd() );
    location.resize( numDoFs );
    std::vector<bool> marker( numDoFs, false );

    // go through the elements of mesh and field simultaneously
    typename FIELD::ElementPtrConstIter fIter = field.elementsBegin();
    typename FIELD::ElementPtrConstIter fEnd  = field.elementsEnd();
    for ( ; fIter != fEnd; ++fIter ) {

        // go through the degrees of freedom of the field element
        typename FIELD::Element::DoFPtrConstIter dIter = (*fIter) -> doFsBegin();
        typename FIELD::Element::DoFPtrConstIter dEnd  = (*fIter) -> doFsEnd();
        for ( unsigned p = 0; dIter != dEnd; ++dIter, p++ ) {

            // Global ID of DoF object
            const std::size_t doFID = (*dIter) -> getID();

            // If DoF has not yet been handled
            if ( not marker[ doFID ] ) {

                // use current element for location
                const std::size_t elemID = (*fIter) -> getID();

                // local coordinate associated with this DoF
                const LocalVecDim xi = supportPoints[ p ];

                // store in provided vector
                location[ doFID ] = std::make_pair( elemID, xi );

                // mark as done
                marker[ doFID ] = true;

            } // condition of untouched DoF
            
        } // loop over DoFs

    } // loop over elements

    return;
}

//------------------------------------------------------------------------------
/** Find a degree of freedom close a point of interest.
 *  Based on the associated location of the degrees of freedom given by the
 *  function base::post::associateLocation , every location in this array is
 *  evaluated with the help of the mesh as a physcial coordinate. Once such a
 *  location is closer to the given point \a x to the given \a tolerance,
 *  the index of that degree of freedom in the array is returned.
 *  \tparam MESH   Type of mesh for geometry representation
 *  \param[in] location  Array of (element ID, coordinate) pairs of all DoFs
 *  \param[in] mesh      The used mesh
 *  \param[in] x         Point of interest
 *  \param[in] tolerance Coordinate comparison tolerance
 */
template<typename MESH>
std::size_t base::dof::findDoFWithLocation( const std::vector<
                                            std::pair<std::size_t,
                                            typename MESH::Element::GeomFun::VecDim> >&
                                            location,
                                            const MESH& mesh, 
                                            const typename MESH::Node::VecDim& x,
                                            const double tolerance )
{
    // go through the given array
    for ( std::size_t p = 0; p < location.size(); p++ ) {

        // extract values
        const typename MESH::Element* gep = mesh.elementPtr( location[p].first );
        const typename MESH::Element::GeomFun::VecDim xi  =  location[p].second;

        // check against point
        const double distance =
            (base::Geometry<typename MESH::Element>()( gep, xi ) - x ).norm();

        if ( distance < tolerance ) return p;
    }

    // this point shall not be reached
    std::cerr << "(WW) Cannot find any DoF closer than " << tolerance
              << " to the point (" << x.transpose() << ") \n\n";

    // return an invalid number
    return base::invalidInt;
}

#endif
