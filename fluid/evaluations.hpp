//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   evaluations.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef fluid_evaluations_hpp
#define fluid_evaluations_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/linearAlgebra.hpp>
#include <base/post/evaluateField.hpp>

//------------------------------------------------------------------------------
namespace fluid{

    //--------------------------------------------------------------------------
    template<unsigned HIST,typename GEOMELEMENT, typename FIELDELEMENT>
    typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,double>::Type
    velocityHistory( const GEOMELEMENT*  geomEp,
                     const FIELDELEMENT* fieldEp,
                     const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return base::post::evaluateFieldHistory<HIST>( geomEp, fieldEp, xi );
    }

    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    typename base::Vector<FIELDELEMENT::DegreeOfFreedom::size,double>::Type
    velocity( const GEOMELEMENT*  geomEp,
              const FIELDELEMENT* fieldEp,
              const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return base::post::evaluateField( geomEp, fieldEp, xi );
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    typename base::Matrix<GEOMELEMENT::Node::dim,
                          FIELDELEMENT::DegreeOfFreedom::size,
                          double>::Type
    velocityGradientHistory( const GEOMELEMENT*  geomEp,
                             const FIELDELEMENT* fieldEp,
                             const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        // return the entry
        return base::post::evaluateFieldGradientHistory<HIST>( geomEp, fieldEp, xi );
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    double pressureHistory( const GEOMELEMENT*  geomEp,
                            const FIELDELEMENT* fieldEp,
                            const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        // return the entry
        return base::post::evaluateFieldHistory<HIST>( geomEp, fieldEp, xi )[0];
    }

    //--------------------------------------------------------------------------
    template<unsigned HIST, typename GEOMELEMENT, typename FIELDELEMENT>
    double velocityDivergenceHistory( const GEOMELEMENT*  geomEp,
                                      const FIELDELEMENT* fieldEp,
                                      const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        // Evaluate the displacement gradient
        const typename base::Matrix<GEOMELEMENT::Node::dim,
                                    FIELDELEMENT::DegreeOfFreedom::size,
                                    double>::Type
            gradU = base::post::evaluateFieldGradientHistory<HIST>( geomEp,
                                                                    fieldEp, xi );
        
        double divU = 0.;
        for ( unsigned d = 0; d < FIELDELEMENT::DegreeOfFreedom::size; d++ )
            divU += gradU(d,d);

        return divU;
    }

    //--------------------------------------------------------------------------
    template<typename GEOMELEMENT, typename FIELDELEMENT>
    double velocityDivergence( const GEOMELEMENT*  geomEp,
                               const FIELDELEMENT* fieldEp,
                               const typename FIELDELEMENT::FEFun::VecDim& xi )
    {
        return velocityDivergenceHistory<0>( geomEp, fieldEp, xi );
    }

    //--------------------------------------------------------------------------
    /** Computation of fluid stress.
     *  Based on a viscous, Newtonian fluid model the stress in the fluid
     *  reads
     *  \f[
     *         \sigma(u,p) = -p I + \mu (\nabla u + (\nabla u)^T)
     *  \f]
     *  based on the pressure \f$ p \f$, the viscosity \f$ \mu \f$ and the
     *  velocity gradient \f$ \nabla u \f$.
     */
    template<typename FIELDTUPLE>
    class Stress
        : public boost::function<
        typename base::Matrix<FIELDTUPLE::GeomElement::Node::dim,
                              FIELDTUPLE::GeomElement::Node::dim>::Type(
                                  const FIELDTUPLE&,
                                  const typename base::GeomTraits<typename
                                  FIELDTUPLE::GeomElement>::LocalVecDim& )>
    {
    public:
        //! Template parameter: a field tuple
        typedef FIELDTUPLE FieldTuple;

        //! @name Derive element pointer types
        //@{
        typedef typename FieldTuple::GeomElement              GeomElement;
        typedef typename FieldTuple::template Binder<1>::Type VelocityElementPtr;
        typedef typename FieldTuple::template Binder<2>::Type PressureElementPtr;
        //@}

        //! Spatial dimension
        static const unsigned dim = GeomElement::Node::dim;

        //! @name Result and coordinate types
        //@{ 
        typedef typename base::Matrix<dim,dim>::Type                ResultType;
        typedef typename base::GeomTraits<GeomElement>::LocalVecDim LocalVecDim;
        //@}

        //! Constructor with viscosity
        Stress( const double viscosity ) : viscosity_( viscosity ) { }

        //----------------------------------------------------------------------
        /** Function call operator computes Cauchy Stress
         *  For a Newtonian fluid with viscosity \f$ \mu \f$, the Cauchy Stress
         *  is given as
         *  \f[
         *      \sigma = - p I + \mu \nabla u + \mu \nabla u^T
         *  \f]
         *  based on the pressure \f$ p \f$ and the velocity \f$ u \f$.
         *  
         *  \param[in] fieldTuple  Tuple of element pointers
         *  \param[in] xi          Local evaluation coordinate
         *  \return                Cauchy stress matrix
         */
        ResultType operator()( const FieldTuple& fieldTuple,
                               const LocalVecDim& xi ) const
        {
            // extract field pointers
            const GeomElement*       geomEp  = fieldTuple.geomElementPtr();
            const VelocityElementPtr velocEp = fieldTuple.template get<1>();
            const PressureElementPtr pressEp = fieldTuple.template get<2>();

            // velocity gradient
            const ResultType
                gradU = base::post::evaluateFieldGradientHistory<0>( geomEp,
                                                                     velocEp, xi );

            // pressure field
            const double p =
                base::post::evaluateFieldHistory<0>( geomEp, pressEp, xi )[0];

            // compute resulting Cauchy stress
            ResultType sigma;
            for ( unsigned d1 = 0; d1 < dim; d1++ ) {
                for ( unsigned d2 = 0; d2 < dim; d2++ ) {

                    sigma(d1,d2) = viscosity_ * ( gradU(d1,d2) + gradU(d2,d1) ) +
                        (d1 == d2 ? -p : 0.0);
                }
            }

            return sigma;
        }
        
    private:
        const double viscosity_;
    };

    //--------------------------------------------------------------------------
    /** Computation of fluid traction on a surface.
     *  Using the fluid::Stress representation of the stress in the fluid,
     *  the traction forces a fluid applies on a surface become
     *  \f[
     *       t(u,p) = \sigma(u,p) n
     *  \f]
     *  using the outward unit normal vector \f$ n \f$.
     */
    template<typename FIELDTUPLE>
    class Traction
        : public boost::function<
        typename base::Vector<FIELDTUPLE::GeomElement::Node::dim>::Type(
            const FIELDTUPLE&,
            const typename base::GeomTraits<typename
            FIELDTUPLE::GeomElement>::LocalVecDim&,
            const typename base::GeomTraits<typename
            FIELDTUPLE::GeomElement>::GlobalVecDim&) >
    {
    public:
        //! Template parameter: a field tuple
        typedef FIELDTUPLE FieldTuple;

        //! @name Derive element pointer types
        //@{
        typedef typename FieldTuple::GeomElement              GeomElement;
        typedef typename FieldTuple::template Binder<1>::Type VelocityElementPtr;
        typedef typename FieldTuple::template Binder<2>::Type PressureElementPtr;
        //@}

        //! Spatial dimension
        static const unsigned dim = GeomElement::Node::dim;

        //! @name Result and coordinate types
        //@{ 
        typedef typename base::Matrix<dim,dim>::Type                MatrixType;
        typedef typename base::Vector<dim>::Type                    ResultType;
        typedef typename base::GeomTraits<GeomElement>::LocalVecDim LocalVecDim;
        typedef typename base::GeomTraits<GeomElement>::LocalVecDim GlobalVecDim;
        //@}

        //! Constructor with viscosity
        Traction( const double viscosity ) : viscosity_( viscosity ) { }

        //----------------------------------------------------------------------
        /** Function call operator computes Cauchy Stress
         *  
         *  \param[in] fieldTuple  Tuple of element pointers
         *  \param[in] xi          Local evaluaton coordinate
         *  \param[in] normal      Normal vector
         *  \return                Cauchy stress matrix
         */
        ResultType operator()( const FieldTuple& fieldTuple,
                               const LocalVecDim& xi,
                               const GlobalVecDim& normal ) const
        {
            const MatrixType sigma =
                fluid::Stress<FieldTuple>( viscosity_ )( fieldTuple, xi );
            ResultType traction;
            traction.noalias() = sigma * normal;
            
            return traction;
        }
        
    private:
        const double viscosity_;
    };

}
#endif
