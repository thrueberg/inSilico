#ifndef advectfield_h
#define advectfield_h

#include <vector>
#include <utility>
#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
/** Copy the values of the field at the previous location to degree of freedom
 *  at current location. This function consists of two steps:
 *  1) Evaluate the field at the previous location of every DoF and store the
 *     value in a temporary vector.
 *  2) Copy the values from this temporary vector to the values of the field's
 *     DoFs.
 *  \tparam MESH  Type of geometry representation
 *  \tparam FIELD Type of field to advect
 */
template<typename MESH, typename FIELD>
void advectField( const MESH& mesh, FIELD& field,
                  const std::vector<std::pair<std::size_t,
                  typename MESH::Element::GeomFun::VecDim> >& previousDoFLocation,
                  const std::vector<double>& supports,
                  const double supportThreshold )
{
    static const unsigned nDoF = FIELD::DegreeOfFreedom::size;
    typedef typename base::Vector<nDoF>::Type       VecDoF;
    typedef typename MESH::Element::GeomFun::VecDim VecDim;

    // zero vector for initialisation
    const  VecDoF zero = base::constantVector<nDoF>( 0. );

    // initialise temporary storage
    const std::size_t numDoFs = previousDoFLocation.size();
    std::vector<VecDoF> tmp( numDoFs, zero );
    
    // Fill storage with solution at previous location
    typename FIELD::DoFPtrIter dIter = field.doFsBegin();
    typename FIELD::DoFPtrIter dEnd  = field.doFsEnd();
    for ( ; dIter != dEnd; ++dIter ) {

        const std::size_t doFID = (*dIter) -> getID();

        // consider only active dofs
        if ( supports[doFID] >= supportThreshold ) {
                    
            const std::size_t elemID = previousDoFLocation[doFID].first;
            const VecDim      prevXi = previousDoFLocation[doFID].second;

            // access to the affected element
            typename MESH::Element*  gOldElem = mesh.elementPtr(  elemID );
            typename FIELD::Element* dOldElem = field.elementPtr( elemID );

            // evaluation of the field at this dof
            const VecDoF old =
                base::post::evaluateField( gOldElem, dOldElem, prevXi );
                    
            // store the dof values in temporary storage
            tmp[ doFID ] = old;

        }
    }

    // copy advected data into dofs
    dIter = field.doFsBegin();
    for ( ; dIter != dEnd; ++dIter ) {
        // global ID of DoF 
        const std::size_t doFID = (*dIter) -> getID();

        // read these values from the tmp storage
        const VecDoF old = tmp[ doFID ];
                
        for ( unsigned d = 0; d < nDoF; d++ ) 
            (*dIter) -> setValue( d, old[ d] );

    } // finished all dofs
    
    return;
}

#endif
