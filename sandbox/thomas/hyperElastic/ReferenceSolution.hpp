#include <base/linearAlgebra.hpp>

//------------------------------------------------------------------------------
template<unsigned DIM>
struct ReferenceSolution
{
    static const unsigned dim = DIM;

    typedef typename base::Vector<dim>::Type     VecDim;
    typedef typename base::Vector<3>::Type       Vec3;
    typedef typename base::Matrix<3,3>::Type     Mat3x3;

    typedef
    boost::array<boost::array<boost::array<double,3>,3>,3> Array3x3x3;


    ReferenceSolution( const double factor )
        : factor_( factor ) { }

    // u[i] = u_i
    Vec3 u( const VecDim& x ) const
    {
        Vec3 result = base::constantVector<3>( 0. );
        result[0] = factor_ * x[0] * x[0] * x[1];
        //result[0] = factor_ * std::sin( M_PI * x[0] );
        return result;
    }

    // dU(i,J) = u_{i,J}
    Mat3x3 dU( const VecDim& x ) const
    {
        Mat3x3 result = base::constantMatrix<3,3>( 0. );
        result(0,0) = 2. * factor_ * x[0] * x[1];
        result(0,1) =      factor_ * x[0] * x[0];
        //result(0,0) = factor_ * M_PI * std::cos( M_PI * x[0] );
        return result;
    }

    // ddU[i][J][K] = u_{i,JK}
    Array3x3x3 ddU( const VecDim& x ) const
    {
        Array3x3x3 result;
        for ( unsigned i = 0; i < 3; i++ )
            for ( unsigned j = 0; j < 3; j++ )
                for ( unsigned k = 0; k < 3; k++ )
                    result[i][j][k] = 0.;

        result[0][0][0] = 2. * factor_ * x[1];
        result[0][0][1] = 2. * factor_ * x[0];
        result[0][1][0] = 2. * factor_ * x[0];
        //result[0][0][0] = -factor_ * M_PI * M_PI * std::sin( M_PI * x[0] );
        
        return result;
    }

    // true if on dirichlet boundary
    bool dirichletBoundaryFilter( const VecDim& x ) const
    {
        return true;
    }

private:
    const double factor_;
};

//------------------------------------------------------------------------------
template<typename REFSOL>
class BoundaryValueProblem
{
public:
    typedef REFSOL ReferenceSolution;
    
    static const unsigned dim = ReferenceSolution::dim;

    typedef typename ReferenceSolution::VecDim         VecDim;
    typedef typename ReferenceSolution::Vec3           Vec3;
    typedef typename ReferenceSolution::Mat3x3         Mat3x3;
    typedef typename ReferenceSolution::Array3x3x3     Array3x3x3;

    BoundaryValueProblem( const double lambda, const double mu,
                          const ReferenceSolution& refSol )
        : lambda_( lambda ), mu_( mu ), refSol_( refSol )
    { }

    //--------------------------------------------------------------------------
    //! Kronecker delta: \f$ \delta_{IJ} = (I==J) \f$
    double delta( const unsigned I, const unsigned J ) const
    {
        return (I==J ? 1. : 0.);
    }

    //--------------------------------------------------------------------------
    //! Deformation gradient \f$ F = I + \nabla_X u \f$
    Mat3x3 F( const VecDim& x ) const
    {
        Mat3x3 F = refSol_.dU( x );
        for ( unsigned d = 0; d < 3; d++ ) F(d,d) += 1.0;
        return F;
    }

    //--------------------------------------------------------------------------
    //! Gradient of deformation gradient: \f$ F_{iJ,K} = u_{i,JK} \f$
    Array3x3x3 gradF( const VecDim& x ) const
    {
        return refSol_.ddU( x );
    }

    //--------------------------------------------------------------------------
    /** Green-Lagrange tensor
     *  \f[
     *       E_{IJ} = \frac{1}{2}( F_{iI} F_{iJ} - \delta_{IJ} )
     *  \f]
     */
    Mat3x3 E( const VecDim& x ) const
    {
        const Mat3x3 F = this -> F( x );
        Mat3x3 C; C.noalias() = F.transpose() * F;
        Mat3x3 E;
        for ( unsigned I = 0; I < 3; I++ ) {
            for ( unsigned J = 0; J < 3; J++ ) {
                E(I,J) = 0.5 * (C(I,J) - delta(I,J));
            }
        }
        
        return E;
    }

    //--------------------------------------------------------------------------
    /** Divergence of the Green-Lagrange straint tensor
     *  \f[
     *      E_{JK,K} = \frac{1}{2}( F_{iJ,K} F_{iK} + F_{iJ} F_{iK,K} )
     *  \f]
     */
    Vec3 DivE( const VecDim& x ) const
    {
        const Array3x3x3 gradF = this -> gradF( x );
        const Mat3x3         F = this ->     F( x );

        Vec3 result;
        for ( unsigned J = 0; J < 3; J++ ) {
            result[J] = 0.;
            for ( unsigned i = 0; i < 3; i++ ) {
                for ( unsigned K = 0; K < 3; K++ ) {
                    result[J] += 0.5 *
                        (gradF[i][J][K] * F(i,K) + F(i,J) * gradF[i][K][K]);
                }
            }
        }

        return result;
    }

    //--------------------------------------------------------------------------
    /** Gradient of the trace of the  Green-Lagrange strain tensor
     *  \f[
     *        E_{JJ,K} = F_{iJ} F_{iJ,K}
     *  \f]
     */
    Vec3 GradTrE( const VecDim& x ) const
    {
        const Array3x3x3 gradF = this -> gradF( x );
        const Mat3x3         F = this ->     F( x );

        Vec3 result;
        for ( unsigned K = 0; K < 3; K++ ) {
            result[K] = 0.;
            for ( unsigned i = 0; i < 3; i++ ) {
                for ( unsigned J = 0; J < 3; J++ ) {
                    result[K] += F(i,J) * gradF[i][J][K];
                }
            }
        }
        return result;
    }

    //--------------------------------------------------------------------------
    /** Second Piola-Kirchhoff stress tensor due to St.Venant-Kirchhoff law:
     *  \f[
     *       S = \lambda tr(E) I + 2 \mu E
     *  \f]
     */
    Mat3x3 S( const VecDim& x ) const
    {
        const Mat3x3 E = this -> E( x );
        double trE = 0.;
        for ( unsigned d = 0; d < 3; d++ ) trE += E(d,d);

        Mat3x3 S = 2. * mu_ * E;
        for ( unsigned d = 0; d < 3; d++ )
            S(d,d) += lambda_ * trE;

        return S;
    }

    //--------------------------------------------------------------------------
    /** Divergence of the second Piola-Kirchhoff stress tensor
     *  \f[
     *      S_{JK,K} = \lambda [tr(E)]_{K} \delta_{JK} + 2 \mu E_{JK,K}
     *  \f]
     */
    Vec3 DivS( const VecDim& x ) const
    {
        const Vec3 divE    = this -> DivE( x );
        const Vec3 gradTrE = this -> GradTrE( x );

        const Vec3 result = lambda_ * gradTrE + 2. * mu_ * divE;
        return result;
    }
    
    //--------------------------------------------------------------------------
    //! First Piola-Kirchhoff stress tensor: \f$ P_{iJ} = F_{iK} S_{KJ} \f$
    Mat3x3 P( const VecDim& x ) const
    {
        return (this -> F(x)) * (this -> S(x));
    }

    //--------------------------------------------------------------------------
    /** Divergence of the first Piola-Kirchhoff stress tensor
     *  \f[
     *      P_{iK,K} = F_{iJ,K} S_{JK} + F_{iJ} S_{JK,K}
     *  \f]
     */
    Vec3 DivP( const VecDim& x ) const
    {
        const Vec3       divS  = this -> DivS( x );
        const Mat3x3     F     = this -> F(    x );
        const Mat3x3     S     = this -> S(    x );
        const Array3x3x3 gradF = this -> gradF( x );

        Vec3 result;
        for ( unsigned i = 0; i < 3; i++ ) {
            result[i] = 0.;

            for ( unsigned J = 0; J < 3; J++ ) {
                result[i] += F(i,J) * divS[J];
            
                for ( unsigned K = 0; K < 3; K++ )
                    result[i] += gradF[i][J][K] * S(J,K);
            }
        }

        return result;
    }

    //--------------------------------------------------------------------------
    //! First Piola-Kirchhoff traction vector \f$ T = P N \f$
    VecDim traction( const VecDim& x, const VecDim& N ) const
    {
        const Mat3x3 P = this -> P( x );

        VecDim result = base::constantVector<dim>( 0. );
        for ( unsigned i = 0; i < dim; i++ ) {
            result[i] = 0.;
            for ( unsigned J = 0; J < dim; J++ )
                result[i] += P(i,J) * N[J];
        }

        return result;
    }

    //--------------------------------------------------------------------------
    //! Body force \f$ f_{i} = - P_{iK,K} \f$
    VecDim bodyForce( const VecDim& x ) const
    {
        const Vec3 divP = this -> DivP( x );
        VecDim result;
        for ( unsigned d = 0; d < dim; d++ ) result[d] = -divP[d];
        return result;
    }

    //--------------------------------------------------------------------------
    //! Apply Dirichlet boundary conditions
    template<typename DOF>
    void dirichletBC( const VecDim& x, DOF* doFPtr ) const
    {
        const bool isDirichlet = refSol_.dirichletBoundaryFilter( x );

        if ( isDirichlet ) {
            const Vec3 u = refSol_.u( x );

            for ( unsigned d = 0; d < dim; d++ )
                if (doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, u[d] );
        }
    }

    //--------------------------------------------------------------------------
    //!
    VecDim solution( const VecDim& x ) const
    {
        VecDim result;
        const Vec3 u = refSol_.u( x );
        for ( unsigned d = 0; d < dim; d++ ) result[d] = u[d];
        return result;
    }

    //--------------------------------------------------------------------------
    //!
    typename base::Matrix<dim,dim>::Type solutionGradient( const VecDim& x ) const
    {
        typename base::Matrix<dim,dim>::Type result;
        const Mat3x3 gradU = refSol_.dU( x );
        for ( unsigned i = 0; i < dim; i++ )
            for ( unsigned j = 0; j < dim; j++ )
                result(i,j) = gradU( j, i );
        return result;
    }

private:
    const double             lambda_;
    const double             mu_;
    const ReferenceSolution& refSol_;
};
