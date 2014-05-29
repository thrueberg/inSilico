//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   ParticleTracer.hpp
//! @author Thomas Rueberg
//! @date   2014

#ifndef base_post_particletracer_hpp
#define base_post_particletracer_hpp

//------------------------------------------------------------------------------
// std  includes
#include <vector>
#include <ostream>
// boost includes
#include <boost/utility.hpp>
// base includes
#include <base/linearAlgebra.hpp>
#include <base/geometry.hpp>
// base/post includes
#include <base/post/evaluateField.hpp>
#include <base/post/findLocation.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace post{

        template<unsigned DIM>
        class ParticleTracer;

    }
}

//------------------------------------------------------------------------------
/** Store a number of nodes and keep track of their location.
 */
template<unsigned DIM>
class base::post::ParticleTracer
    : boost::noncopyable
{
public:
    //! Coordinate type
    typedef typename base::Vector<DIM>::Type VecDim;

    //! Constructor sets up the vector of vector structure
    ParticleTracer()
    {
        std::vector<VecDim> dummy;
        traces_.push_back( dummy );
    }

    //! Store a point for tracing
    void registerPoint( const VecDim& x )
    {
        traces_[0].push_back( x );
        valid_.push_back( true );
    }

    //! Return total number of registered points
    std::size_t numPoints() const
    {
        return valid_.size();
    }

    //! Number of tracing steps performed
    std::size_t numSteps() const
    {
        return traces_.size();
    }

    //--------------------------------------------------------------------------
    //! Find all latest point locations and evaluate their dislocation
    template<typename MESH,typename FIELD>
    void update( const MESH& mesh, const FIELD& field,
                 const double tolerance, const unsigned maxIter )
    {
        // find latest points in mesh
        const std::size_t nS = this -> numSteps();
        const std::size_t nP = this -> numPoints();

        std::vector<VecDim> latest = traces_[ nS-1 ];

        for ( std::size_t p = 0; p < nP; p++ ) {

            if ( valid_[p] ) {

                std::pair<std::size_t,
                          typename MESH::Element::GeomFun::VecDim> result;

                const bool found = 
                    base::post::findLocationInMesh( mesh, latest[p], tolerance,
                                                    maxIter, result );

                if ( found ) {
                    const VecDim increment =
                        base::post::evaluateField( mesh.elementPtr( result.first ),
                                                   field.elementPtr( result.first ),
                                                   result.second );
                    latest[p] += increment;
                }
                else {
                    valid_[p] = false;
                }
                
            }

        }

        traces_.push_back( latest );
    }

    //--------------------------------------------------------------------------
    void writeLatest( std::ostream& out ) const
    {
        const std::size_t nS = this -> numSteps();
        const std::size_t nP = this -> numPoints();
        for ( std::size_t p = 0; p < nP; p++ ) {
            for ( unsigned d = 0; d < DIM; d++ )
                out << traces_[nS-1][p][d] << "  ";
        }
        out << "\n";

        return;
    }

    //--------------------------------------------------------------------------
    void flush( std::ostream& out ) const
    {
        const std::size_t nS = this -> numSteps();
        const std::size_t nP = this -> numPoints();
        for ( std::size_t n = 0; n < nS; n++ ) {
            for ( std::size_t p = 0; p < nP; p++ ) {
                for ( unsigned d = 0; d < DIM; d++ )
                    out << traces_[n][p][d] << "  ";
            }
            out << "\n";
        }

        return;
    }

private:
    std::vector< std::vector<VecDim> > traces_;
    std::vector<bool>                  valid_;
};

#endif
