
#include <iostream>
#include <cmath>
#include <boost/array.hpp>
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
template<typename PHASE, unsigned N>
class Reactor
    : boost::noncopyable
{
public:
    typedef typename base::Vector<PHASE::DegreeOfFreedom::size>::Type VecDof;
    
    Reactor( const boost::array<PHASE*,   N>& phases, 
             const boost::array<unsigned, N>& stoichiometricA,
             const boost::array<unsigned, N>& stoichiometricB, 
             const double phiF,
             const double kF, const double kB )
        : phases_(           phases ),
          stoichiometricA_(  stoichiometricA ),
          stoichiometricB_(  stoichiometricB ),
          phiF_( phiF ),
          kF_( kF ), kB_( kB )
    {
        std::cout << "#  Chemical reaction \n"
                  << "# \n#";
        for ( unsigned n = 0; n < N; n++ ) {
            std::cout << stoichiometricA_[n] << " A_" << n << " ";
            if ( n < N-1 ) std::cout << "+  ";
        }

        std::cout << " <---> ";

        for ( unsigned n = 0; n < N; n++ ) {
            std::cout << stoichiometricB_[n] << " B_" << n << " ";
            if ( n < N-1 ) std::cout << "+  ";
        }

        std::cout << "\n# \n"
                  << "# Forward rate: " << kF_
                  << ", backward rate: " << kB_
                  << "\n#\n";
    }

    template<unsigned K, typename GEOMELEMENT>
    VecDof massProduction( const GEOMELEMENT* gEp,
                           const typename GEOMELEMENT::GeomFun::VecDim& xi ) const
    {
        boost::array<double,N> concentrations;
        for ( unsigned n = 0; n < N; n ++ )
            concentrations[n] = ( this -> concentration_( gEp, xi, phases_[n] ) )[0];

        const double factor = (static_cast<double>( stoichiometricB_[K] ) -
                               static_cast<double>( stoichiometricA_[K] ) );

        double a = 1.0;
        double b = 1.0;

        for ( unsigned n = 0; n < N; n++ ) {
            a *= std::pow( (phiF_ * concentrations[n]), stoichiometricA_[n] );
            b *= std::pow( (phiF_ * concentrations[n]), stoichiometricB_[n] );
        }


        VecDof result;
        result[0] = factor * (kF_ * a - kB_ * b );

        //if ( (gEp->getID()) == 0 ) {
        //    std::cout << " " << result << std::endl
        //              << " <- "
        //              << concentrations[0] << "  " << concentrations[1] << "  "
        //              << concentrations[2] << "  " << std::endl;
        //}

        return result;
    }
    

private:
    template<typename GEOMELEMENT>
    VecDof concentration_( const GEOMELEMENT* gEp,
                           const typename GEOMELEMENT::GeomFun::VecDim& xi,
                           const PHASE* phase) const
    {
        typename PHASE::Element* fieldEp = phase -> elementPtr( gEp -> getID() );

        return base::post::evaluateField( gEp, fieldEp, xi );
    }

private:
    const boost::array<PHASE*,   N> phases_;
    const boost::array<unsigned, N> stoichiometricA_;
    const boost::array<unsigned, N> stoichiometricB_;

    const double phiF_;
    const double kF_;
    const double kB_;
};
