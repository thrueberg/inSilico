//------------------------------------------------------------------------------
// <preamble>
// </preamble>
//------------------------------------------------------------------------------

//! @file   Basis.hpp
//! @author Thomas Rueberg
//! @date   2013

#ifndef base_fe_basis_hpp
#define base_fe_basis_hpp

//------------------------------------------------------------------------------
// base includes
#include <base/shape.hpp>
#include <base/funSpace.hpp>
#include <base/LagrangeShapeFun.hpp>
// base/sfun includes
#include <base/sfun/BSpline.hpp>
// base/fe includes
#include <base/fe/LagrangeElement.hpp>
#include <base/fe/BSplineCell.hpp>

//------------------------------------------------------------------------------
namespace base{
    namespace fe{
        
        //----------------------------------------------------------------------
        namespace detail_{

            //------------------------------------------------------------------
            //! Declare a default value for the across-element continuity
            template<unsigned DEGREE, base::FunSpace FESPACE>
            struct DefaultContinuity;


            //! By default, Lagrange shape functions are C^0 
            template<unsigned DEGREE>
            struct DefaultContinuity<DEGREE,base::LAGRANGE>
            {
                static const int value = 0;
            };

            //! By default, B-Splines are C^(degree-1) across elements
            template<unsigned DEGREE>
            struct DefaultContinuity<DEGREE,base::BSPLINE>
            {
                static const int value = DEGREE-1;
            };

            //! By default, Hermite polynomials are C^1 across elements
            template<unsigned DEGREE>
            struct DefaultContinuity<DEGREE,base::HERMITE>
            {
                static const int value = 1;
            };
        }

        //----------------------------------------------------------------------
        namespace detail_{

            template<base::Shape SHAPE, unsigned DEGREE, int CONTINUITY,
                     base::FunSpace FESPACE>
            struct FiniteElementTypeBinder;

            template<base::Shape SHAPE, unsigned DEGREE, int CONTINUITY>
            struct FiniteElementTypeBinder<SHAPE,DEGREE,CONTINUITY,
                                           base::LAGRANGE>
            {
                typedef base::fe::LagrangeElement<SHAPE,DEGREE> Type;
            };

            template<base::Shape SHAPE,unsigned DEGREE, int CONTINUITY>
            struct FiniteElementTypeBinder<SHAPE,DEGREE,CONTINUITY,
                                           base::BSPLINE>
            {
                static const unsigned dim = base::ShapeDim<SHAPE>::value;

                typedef base::fe::BSplineCell<dim,DEGREE,CONTINUITY> Type;
            };
            
            
            // Possibly: - HermiteElement
            //           - RaviartThomasElement

        }

        //----------------------------------------------------------------------
        template< base::Shape    SHAPE,
                  unsigned       DEGREE     = 1,
                  base::FunSpace FESPACE    = base::LAGRANGE,
                  int            CONTINUITY =
                  base::fe::detail_::DefaultContinuity<DEGREE,FESPACE>::value >
        struct Basis;
        
        //----------------------------------------------------------------------
        /** Finite element basis for isoparametric analysis.
         *  The fashionable term isoparamtric implies that the basis functions
         *  of the analysis are the same that are used for the geometry
         *  representation. Here, a mesh type is introspected for its geomtry
         *  representation and this types are used to define the analysis basis.
         *  \tparam MESH Type of mesh which determines the analysis
         */
        template<typename MESH>
        struct IsoparametricBasis
            : public Basis<MESH::Element::shape,
                           MESH::Element::GeomFun::degree,
                           MESH::Element::GeomFun::funSpace>
        {
            // empty, all is inherited
        };



        
    }
}

//------------------------------------------------------------------------------
/** Descriptor object for the chosen finite element basis.
 *  A Finite Element (Ciarlet, 1978) is given by the triplet
 *  (Brenner, Scott 2008)
 *  \f[
 *          F = {K, P, N}
 *  \f]
 *  with the following entities: \f$ K \f$ is the element domain, in our case
 *  represented by base::Shape (or an image thereof), \f$ P \f$ denotes the
 *  space of shape functions and \f$ N \f$ are the nodal variables. Commonly,
 *  \f$ P \f$ implies \f$ N \f$ as, for instance, the use of standard Lagrange
 *  shape functions is combined with evaluation of the functions at the
 *  so-called nodes. This nodal evaluations are the functionals which form the
 *  basis of \f$ N \f$ and we would get \f$ N_i(\phi_j ) = \delta_{ij} \f$,
 *  i.e. the nodal evaluation at node \f$ j \f$ of shape function \f$ i \f$
 *  gives one if and only if \f$ i = j \f$ and zero otherwise. Using other
 *  bases, e.g. a B-Spline basis, the dual space is not so obvious anymore.
 *  In any case, we assume that it is well-defined.
 *  
 *  In addition to this standard triplet, the inter-element continuity shall
 *  play a significan role and this global characteristic of the interpolation
 *  operator is included in this basis description.
 *
 *  \tparam SHAPE     Element domain shape (\f$ K \f$)
 *  \tparam DEGREE,
 *          FESPACE   Descriptors of the local shape functions on \f$ K \f$
 *  \tparam CONTINUITY Inter-element continuity of the global interpolation
 */
template< base::Shape SHAPE, unsigned DEGREE,
          base::FunSpace FESPACE, int CONTINUITY>
struct base::fe::Basis
{
    //! @name Template parameters
    //@{
    static const base::Shape    shape      = SHAPE;      //!< Geometric shape
    static const unsigned       degree     = DEGREE;     //!< Polynomial degree 
    static const base::FunSpace feSpace    = FESPACE;    //!< Classifier 
    static const int            continuity = CONTINUITY; //!< C^x of fe fun
    //@}

    //! Deduce type of finite elements
    typedef typename
    base::fe::detail_::FiniteElementTypeBinder<shape,degree,
                                               continuity,feSpace>::Type
    FiniteElement;

    //! Deduce type of shape functions used for this FE basis
    typedef typename FiniteElement::ShapeFun FEFun;
    
};

#endif
