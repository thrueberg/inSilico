//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Constraint.hpp
//! @author Thomas Rueberg
//! @date   2012

#ifndef base_dof_constraint_hpp
#define base_dof_constraint_hpp

//------------------------------------------------------------------------------
// std   includes
#include <vector>
// boost includes
#include <boost/tuple/tuple.hpp>
// base  includes
#include <base/numbers.hpp>
// boost/io/includes
#include <base/io/Format.hpp>

//------------------------------------------------------------------------------
namespace base{

    namespace dof{
        
        template<typename DOF>
        class Constraint;
    }
}

//------------------------------------------------------------------------------
/** Representation of a linear constraint.
 *  A constrained degree of freedom component \f$ u^i_j \f$, where \f$ i \f$ is
 *  the identifier of the degree of freedom and \f$ j \f$ is the component of
 *  interest, obtains the relation
 *  \f[
 *       u^i_j = g^i_j + \sum_{k,l} c^{ik}_{jl} u^k_l
 *  \f]
 *  where \f$ g \f$ is the non-homogenous part and the matrix \f$ c \f$ gives
 *  the weights with which other degrees of freedom contribute to this one.
 *
 *  Internally, tuples of a pointer to a different degree of freedom \f$ k \f$
 *  together with the direction \f$ l \f$ and the value of \f$ c \f$ are
 *  stored in a dynamic list.
 *  \tparam DOF Type of degree contributing to this one
 */
template<typename DOF>
class base::dof::Constraint
{
public:
    //! Template parameter: type of DoF
    typedef DOF DegreeOfFreedom;

    //! Representation of one weighted DoF
    typedef boost::tuple<DegreeOfFreedom*,unsigned,base::number> WeightedDoF;

    //! @name Constructors
    //@{
    Constraint() : rhs_( 0. ) { }

    Constraint( const base::number given ) : rhs_( given ) { }

    Constraint( DegreeOfFreedom* doF,
                const unsigned dir, const base::number weight,
                const base::number rhs )
        : rhs_( rhs )
    {
        weightedDoFs_.push_back( boost::make_tuple( doF, dir, weight ) );
    }
    //@}

    void setValue( const base::number value ) { rhs_ = value; }
    
    //! Add another weighted dof
    void addWeightedDoF( DegreeOfFreedom* doF, const unsigned dir,
                         base::number weight )
    {
        ASSERT_MSG( doF -> isActive( dir ),
                    "DoF with ID " + x2s(doF->getID()) + " not active ");
        weightedDoFs_.push_back( boost::make_tuple( doF, dir, weight ) );
    }

    //! Evaluate the constraint equation
    base::number evaluate( const bool useRhsTerm = true ) const
    {
        base::number result = ( useRhsTerm ? rhs_ : 0. );
        for ( unsigned w = 0; w < weightedDoFs_.size(); w++ ) {
            const unsigned dir  = weightedDoFs_[w].template get<1>();
            const double weight = weightedDoFs_[w].template get<2>();
            result += weight *
                (weightedDoFs_[w].template get<0>()) -> getValue( dir );

            // guarantee that the master DoF is active
            ASSERT_MSG( (weightedDoFs_[w].template get<0>()) -> isActive( dir ),
                        "DoF with ID " +
                        x2s( (weightedDoFs_[w].template get<0>()) -> getID() ) +
                        " is not active" );

        }
        
        return result;
    }

    //! Get just the rhs value
    base::number getValue() const { return rhs_; }

    //! Copy the IDs of the weighted dofs
    void getWeightedDoFIDs( std::vector< std::pair<base::number,std::size_t> >&
                            weightedDoFIDs ) const
    {
        for ( unsigned w = 0; w < weightedDoFs_.size(); w++ ) {
            const unsigned     dir    = weightedDoFs_[w].template get<1>();
            const base::number weight = weightedDoFs_[w].template get<2>();
            const std::size_t id =
                weightedDoFs_[w].template get<0>() -> getIndex( dir );

            // guarantee that the master DoF is active
            ASSERT_MSG( (weightedDoFs_[w].template get<0>()) -> isActive( dir ),
                        "DoF with ID " +
                        x2s( (weightedDoFs_[w].template get<0>()) -> getID() ) +
                        " is not active" );

            weightedDoFIDs.push_back( std::make_pair( weight, id ) );
        }
        return;
    }

private:
    base::number             rhs_;          //!< RHS term of linear constraint
    std::vector<WeightedDoF> weightedDoFs_; //!< Weighted contributing DoFs
};

#endif

