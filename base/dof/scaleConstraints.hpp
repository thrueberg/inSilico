//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   scaleConstraints.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_dof_scaleConstraints_hpp
#define base_dof_scaleConstraints_hpp

//------------------------------------------------------------------------------
// base
#include <base/types.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace dof{

        template<typename DOFITER>
        void scaleConstraints( DOFITER first, DOFITER last, const double factor )
        {
            // size of individual dof object deduced from the iterator type
            static const unsigned size =
                base::TypeReduction<typename DOFITER::value_type>::Type::size;

            for ( DOFITER dIter = first; dIter != last; ++dIter ) {
                
                for ( unsigned d = 0; d < size; d ++ ) {
                    if ( not (*dIter) -> isActive( d ) ) {
                        const double oldValue = (*dIter) -> getConstraint( d );
                        //(*dIter) -> setValue( d, factor*oldValue );
                        (*dIter) -> constrainValue( d, factor * oldValue );
                    }
                }

            }
            return;
        }
                               

        
    }
}

#endif
