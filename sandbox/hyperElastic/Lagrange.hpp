#include <mat/hypel/CompNeoHookean.hpp>

//------------------------------------------------------------------------------
class AnalyticLagrange1D
{
public:

    AnalyticLagrange1D( ) { }
    
    void setFactor( const double ubar )
    {
        alpha_ = ubar;
        beta_  = 3. * alpha_ * alpha_ + 3. * alpha_ + 1.;
    }

public:
    
    double x( const double X ) const
    {
        return ( std::pow(X+alpha_, 3.) - alpha_*alpha_*alpha_) / beta_;
    }

    double u( const double X ) const
    {
        return x(X) - X;
    }

    double F( const double X ) const
    {
        return 3./beta_ * (X+alpha_)*(X+alpha_);
    }

    double dF( const double X ) const
    {
        return 6./beta_ * (X+alpha_);
    }

private:
    double alpha_;
    double beta_;
};

//------------------------------------------------------------------------------
template<unsigned DIM>
class AnalyticLagrangeTensor
{
public:
    static const unsigned dim = DIM;
    typedef typename base::Vector<dim>::Type VecDim;
    
    AnalyticLagrangeTensor( const double lambda, const double mu,
                            boost::array<double,DIM> ubar )
        : lambda_( lambda ), mu_( mu )
    {
        for ( unsigned d = 0; d < DIM; d++ )
            oneDimensional_[d].setFactor( ubar[d] );
    }

public:
    mat::Tensor F( const VecDim& X ) const
    {
        mat::Tensor F = mat::Tensor::Identity();
        for ( unsigned d = 0; d < DIM; d++ )
            F(d,d) = oneDimensional_[d].F( X[d] );

        return F;
    }

    mat::Tensor P( const VecDim& X ) const
    {
        mat::Tensor FF = F( X );
        const double J = mat::determinant( FF );

        const double alpha = lambda_ * std::log( J ) - mu_;
        
        mat::Tensor P = alpha * FF.transpose() + mu_ * FF;

        return P;
    }
    
    VecDim force( const VecDim& X ) const
    {
        VecDim FF, DF;
        double J = 1.;
        for ( unsigned d = 0; d < dim; d++ ) {
            FF[d] = oneDimensional_[d].F(  X[d] );
            DF[d] = oneDimensional_[d].dF( X[d] );
            J *= FF[d];
        }

        VecDim result;
        for ( unsigned d = 0; d < dim; d++ )
            result[d] = -DF[d] / FF[d] / FF[d] *
                (lambda_ * (1. - std::log(J)) + mu_ * (1. + FF[d]*FF[d]));

        //result[d] = -DF[d] * mu_ * (1. + 1./FF[d]/FF[d]);
            

        return result;
    }

    VecDim solution( const typename base::Vector<DIM>::Type X ) const
    {
        VecDim result;
        for ( unsigned d = 0; d < dim; d++ )
            result[d] = oneDimensional_[d].u( X[d] );
        
        return result;
    }
    
private:
    const double lambda_;
    const double mu_;

    boost::array<AnalyticLagrange1D,dim> oneDimensional_;
};

