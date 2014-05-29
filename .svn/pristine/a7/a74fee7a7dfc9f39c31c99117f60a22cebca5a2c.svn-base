//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Parameters.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_nitsche_parameters_hpp
#define base_nitsche_parameters_hpp

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace base{
    namespace nitsche{

        template<typename DUMMY>
        class Parameters;

        class PrescribedParameters;

        class OuterBoundary;

        template<typename CELL>
        class ImmersedBoundary;

        template<typename CELL>
        class ImmersedInterface;
    }
}

//------------------------------------------------------------------------------
/**  
 *   The inverse estimate reads
 *   \f[
 *        \| \alpha \partial_n u \|^2_{f(e)} \leq C^2_e  a(u,u)
 *   \f]
 *   for the element \f$ \tau_e \f$ and its boundary face \f$ f(e) \f$. 
 */
template<typename DUMMY>
class base::nitsche::Parameters
{
    STATIC_ASSERT_MSG( sizeof(DUMMY) == 0,
                       "Documentation class only" );
};

//------------------------------------------------------------------------------
/**  Nitsche's method for boundary conditions on the mesh boundary.
 */
class base::nitsche::PrescribedParameters
{
public:
    PrescribedParameters( const double alpha, const double kappa )
        : alpha_( alpha ), kappa_( kappa ) { }

    template<typename ITER>
    double penaltyWeight( const ITER dummy ) const
    {
        return alpha_;
    }

    template<typename ITER>
    double energyWeight( const ITER dummy,
                         const bool inOut ) const
    {
        return kappa_;
    }

private:
    const double alpha_;
    const double kappa_;
};

//------------------------------------------------------------------------------
/**  Nitsche's method for boundary conditions on the mesh boundary.
 */
class base::nitsche::OuterBoundary
{
public:
    OuterBoundary( const double alpha )
        : alpha_( alpha ) { }

    template<typename ITER>
    double penaltyWeight( const ITER dummy ) const
    {
        return alpha_;
    }

    template<typename ITER>
    double energyWeight( const ITER dummy,
                         const bool inOut ) const
    {
        VERIFY_MSG( inOut, "Only one domain possible" );
        return 1.0;
    }

private:
    const double alpha_;
};

//------------------------------------------------------------------------------
/**  Nitsche's method for boundary conditions on an immersed boundary
 */
template<typename CELL>
class base::nitsche::ImmersedBoundary
{
public:
    ImmersedBoundary( const double alpha,
                      const std::vector<CELL>& cells )
        : alpha_( alpha ),
          cells_( cells ) { }

    template<typename ITER>
    double penaltyWeight( const ITER iter ) const
    {
        const bool inside = true;
        const std::size_t elemID = iter -> geomElementPtr() -> getID();
        
        const double sizeSurface = cells_[elemID].surfaceArea();
        const double sizeVolume  = cells_[elemID].parameterArea( inside );

        VERIFY_MSG( cells_[elemID].isCut(), "something wrong here" );

        return alpha_ * sizeSurface / sizeVolume;
    }

    template<typename ITER>
    double energyWeight( const ITER dummy, const bool inOut ) const
    {
        return (inOut ? 1.0 : -1.0);
    }
    
private:
    const double             alpha_;
    const std::vector<CELL>& cells_;
};

//------------------------------------------------------------------------------
/**  Nitsche's method for conditions on an immersed interface
 */
template<typename CELL>
class base::nitsche::ImmersedInterface
{
public:
    ImmersedInterface( const double alpha1, const double alpha2, 
                       const std::vector<CELL>& cells )
        : alpha1_( alpha1 ),
          alpha2_( alpha2 ),
          cells_(  cells ) { }

    template<typename ITER>
    double penaltyWeight( const ITER iter ) const
    {
        const std::size_t elemID = iter -> geomElementPtr() -> getID();

        const bool inside = true;
        const double sizeSurface = cells_[elemID].surfaceArea();
        const double sizeVolume1 = cells_[elemID].parameterArea( inside );
        const double sizeVolume2 = cells_[elemID].parameterArea( not inside );

        VERIFY_MSG( cells_[elemID].isCut(), "something wrong here" );

        return alpha1_ * alpha2_ * sizeSurface /
            (sizeVolume1 * alpha2_ + sizeVolume2 * alpha1_);
    }

    template<typename ITER>
    double energyWeight( const ITER iter,
                         const bool inOut ) const
    {
        const std::size_t elemID =  iter -> geomElementPtr() -> getID();

        const bool inside = true;
        const double sizeVolume1  = cells_[elemID].parameterArea( inside );
        const double sizeVolume2  = cells_[elemID].parameterArea( not inside );

        const double numer = (inOut ? sizeVolume1 * alpha2_ : sizeVolume2 * alpha1_ );
        const double denom = (sizeVolume1 * alpha2_ + sizeVolume2 * alpha1_);

        const double kappa = numer / denom;

        return kappa;
    }
    
private:
    const double             alpha1_;
    const double             alpha2_;
    const std::vector<CELL>& cells_;
};

#endif
