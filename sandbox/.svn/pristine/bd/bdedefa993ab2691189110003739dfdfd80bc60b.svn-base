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

    void setFactor( const double factor ) { factor_ = factor; }

    double getFactor( ) const { return factor_; }

    // u[i] = u_i
    Vec3 u( const VecDim& X ) const
    {
        Vec3 result = base::constantVector<3>( 0. );

        //result[0] = factor_ * X[0] * X[0] * X[1];
        result[0] += factor_ * X[0] * X[0];
        return result;
    }

    Vec3 x( const VecDim& X ) const
    {
        Vec3 result = base::constantVector<3>( 0. );
        Vec3 u = (this -> u(X));
        for ( unsigned d = 0; d < dim; d++ ) result[d] = X[d] + u[d];
        return result;
    }

    // X(x) = x - u
    VecDim X( const VecDim& x ) const
    {
        VecDim result = x;
        
        //if (std::abs( x[1] < 1.e-10 ) ) return result;
        //result[0] = (std::sqrt(1. + 4.*factor_*x[0]*x[1]) - 1.) / (2. * factor_ * x[1]);

        result[0] = .5/factor_ * (std::sqrt(4.*factor_*x[0] +1.) -1.);
        return result;
    }

    // dU(i,J) = u_{i,J}
    Mat3x3 dU( const VecDim& X ) const
    {
        Mat3x3 result = base::constantMatrix<3,3>( 0. );
        //result(0,0) = 2. * factor_ * X[0] * X[1];
        //result(0,1) =      factor_ * X[0] * X[0];
        result(0,0) = 2. * factor_ * X[0];
        return result;
    }

    // ddU[i][J][K] = u_{i,JK}
    Array3x3x3 ddU( const VecDim& X ) const
    {
        Array3x3x3 result;
        for ( unsigned i = 0; i < 3; i++ )
            for ( unsigned j = 0; j < 3; j++ )
                for ( unsigned k = 0; k < 3; k++ )
                    result[i][j][k] = 0.;

        //result[0][0][0] = 2. * factor_ * X[1];
        //result[0][0][1] = 2. * factor_ * X[0];
        //result[0][1][0] = 2. * factor_ * X[0];
        result[0][0][0] = 2. * factor_;
        
        return result;
    }

    // true if on dirichlet boundary
    bool dirichletBoundaryFilter( const VecDim& X ) const
    {
        return true;
    }

private:
    double factor_;
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
                          const ReferenceSolution& refSol,
                          const bool lagrange )
        : lambda_( lambda ), mu_( mu ), refSol_( refSol ), lagrange_( lagrange )
    { }

private:
    //--------------------------------------------------------------------------
    //! Kronecker delta: \f$ \delta_{IJ} = (I==J) \f$
    double delta_( const unsigned I, const unsigned J ) const
    {
        return (I==J ? 1. : 0.);
    }

    //--------------------------------------------------------------------------
    //! Deformation gradient \f$ F = I + \nabla_X u \f$
    Mat3x3 F_( const VecDim& X ) const
    {
        Mat3x3 F = refSol_.dU( X );
        for ( unsigned d = 0; d < 3; d++ ) F(d,d) += 1.0;
        return F;
    }

    //--------------------------------------------------------------------------
    //! Gradient of deformation gradient: \f$ F_{iJ,K} = u_{i,JK} \f$
    Array3x3x3 gradF_( const VecDim& X ) const
    {
        return refSol_.ddU( X );
    }

    //--------------------------------------------------------------------------
    /** Green-Lagrange tensor
     *  \f[
     *       E_{IJ} = \frac{1}{2}( F_{iI} F_{iJ} - \delta_{IJ} )
     *  \f]
     */
    Mat3x3 E_( const VecDim& X ) const
    {
        const Mat3x3 F = this -> F_( X );
        Mat3x3 C; C.noalias() = F.transpose() * F;
        Mat3x3 E;
        for ( unsigned I = 0; I < 3; I++ ) {
            for ( unsigned J = 0; J < 3; J++ ) {
                E(I,J) = 0.5 * (C(I,J) - (this -> delta_(I,J) ) );
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
    Vec3 DivE_( const VecDim& X ) const
    {
        const Array3x3x3 gradF = this -> gradF_( X );
        const Mat3x3         F = this ->     F_( X );

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
    Vec3 GradTrE_( const VecDim& X ) const
    {
        const Array3x3x3 gradF = this -> gradF_( X );
        const Mat3x3         F = this ->     F_( X );

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
    Mat3x3 S_( const VecDim& X ) const
    {
        const Mat3x3 E = this -> E_( X );
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
    Vec3 DivS_( const VecDim& X ) const
    {
        const Vec3 divE    = this ->    DivE_( X );
        const Vec3 gradTrE = this -> GradTrE_( X );

        const Vec3 result = lambda_ * gradTrE + 2. * mu_ * divE;
        return result;
    }
    
    //--------------------------------------------------------------------------
    //! First Piola-Kirchhoff stress tensor: \f$ P_{iJ} = F_{iK} S_{KJ} \f$
    Mat3x3 P_( const VecDim& X ) const
    {
        return (this -> F_(X)) * (this -> S_(X));
    }

    //--------------------------------------------------------------------------
    /** Divergence of the first Piola-Kirchhoff stress tensor
     *  \f[
     *      P_{iK,K} = F_{iJ,K} S_{JK} + F_{iJ} S_{JK,K}
     *  \f]
     */
    Vec3 DivP_( const VecDim& X ) const
    {
        const Vec3       divS  = this ->  DivS_( X );
        const Mat3x3     F     = this ->     F_( X );
        const Mat3x3     S     = this ->     S_( X );
        const Array3x3x3 gradF = this -> gradF_( X );

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
    /** Compute lagrangian coordinate if necessary: \f$ X = x - u \f$
     */
    VecDim coord_( const VecDim& xX ) const
    {
        if ( lagrange_ ) return xX;

        return refSol_.X( xX );
    }

public:    

    //--------------------------------------------------------------------------
    //! First Piola-Kirchhoff traction vector \f$ T = P N \f$
    VecDim traction( const VecDim& xX, const VecDim& N ) const
    {
        const VecDim X = this -> coord_( xX );
        
        const Mat3x3 P = this -> P_( X );

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
    VecDim bodyForce( const VecDim& xX ) const
    {
        const VecDim X = this -> coord_( xX );
        
        const Vec3 divP = this -> DivP_( X );
        VecDim result;
        for ( unsigned d = 0; d < dim; d++ ) result[d] = -divP[d];

        //if ( false ) {
        //    const double alpha = refSol_.getFactor();
        //    std::cout << xX.transpose() << ":  "
        //              << -(lambda_ + 2 * mu_) *
        //        (12.*alpha*alpha*alpha*X[0]*X[0] + 12.*alpha*alpha*X[0] + 2.*alpha)
        //              << " vs. " << result[0] << std::endl;
        //}

        
        return result;
    }

    //--------------------------------------------------------------------------
    //! Apply Dirichlet boundary conditions
    template<typename DOF>
    void dirichletBC( const VecDim& xX, DOF* doFPtr ) const
    {
        const VecDim X = xX; //this -> coord_( xX );
        
        const bool isDirichlet = refSol_.dirichletBoundaryFilter( X );

        const Vec3 u = refSol_.u( X );

        if ( isDirichlet ) {

            for ( unsigned d = 0; d < dim; d++ )
                if (doFPtr -> isActive(d) ) doFPtr -> constrainValue( d, u[d] );
        }

    }

    //--------------------------------------------------------------------------
    //!
    VecDim solution( const VecDim& xX ) const
    {
        const VecDim X = this -> coord_( xX );
        
        VecDim result;
        const Vec3 x = refSol_.x( X );
        for ( unsigned d = 0; d < dim; d++ ) result[d] = x[d];

        return result;
    }

    //--------------------------------------------------------------------------
    //!
    typename base::Matrix<dim,dim>::Type solutionGradient( const VecDim& xX ) const
    {
        const VecDim X = this -> coord_( xX );
        
        typename base::Matrix<dim,dim>::Type result;
        const Mat3x3 gradU = refSol_.dU( X );
        for ( unsigned i = 0; i < dim; i++ )
            for ( unsigned j = 0; j < dim; j++ )
                result(i,j) = gradU( j, i );
        return result;
    }

private:
    const double             lambda_;
    const double             mu_;
    const ReferenceSolution& refSol_;
    const bool               lagrange_;
};
