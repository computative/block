#include "blocker.h"
#include <cmath>

std::vector<double> operator-(std::vector<double> v, double scalar )
{
    for(double& vi : v) vi = vi - scalar;
    return v;
}

Blocker::Blocker(std::vector<double>& x)
{
    setData(x);
}

void Blocker::setData(std::vector<double>& x)
{
    unsigned long long n = x.size();
    // checks if size of input vector is a power of two
    if ( !isPowerOfTwo(n) )
    {
        // Error
        stdErr = mean = mse_stdErr = mse_mean = -1.0f;
        nK = 0;
        printf( "Error: Data size must be a power of two." );
    }
    else {
        // success. Proceeds to compute results
        mean = this->avg(x);
        mse_mean = 0 + estimate(x, n);
        stdErr = sqrt(mse_mean);
        mse_stdErr = (2*nK-1)*mse_mean*mse_mean/(nK*nK);
    }
}

// checks that integer n is a power of 2
bool Blocker::isPowerOfTwo(unsigned long long n)
{
    return (n>0 && ((n & (n-1)) == 0));
}

// the algorithm which computes the variance of the sample mean.
double Blocker::estimate(std::vector<double>& x, unsigned long long& n)
{
    // since n is a power of two, d is an integer.
    int d = log2(n);
    std::vector<double> X = x - mean;
    std::vector<double> s(d, 0.0f);
    std::vector<double> gamma(d, 0.0f);
    // computes covariance and variance and applies blocking transform
    for (int k = 0; k <= d-1; k++) {
        gamma[k] = this->gamma1(X, n);
        s[k] = this->var(X, n);
        this->transform(X, x, n);
    }
    // computes the test statistic Mj
    std::vector<double> M(d, 0.0f);
    for(int j = 0; j <= d-1; j++) {
        M[j] = pow(gamma[j]/s[j],2 )*pow(2,d-j);
    }
    M = this->cumsum(this->flip(M));
    // finds the first k such that Mk < quantile[k]
    int k;
    for(k = d-1; k >= 0; k--) {
        if(M[k] < quantile[k])
            break;
    }
    int K = d-(k+1);
    // returns answer
    nK = pow(2, d-K);
    return s[K]/nK;
}

// performs blocking transformation
void Blocker::transform(std::vector<double>& X, std::vector<double>& x, unsigned long long& n)
{
    for (unsigned long long i = 0; i < n/2; i++) {
        x[i] = 0.5*( x[2*i] + x[2*i+1] );
        X[i] = x[i] - mean;
    }
    n = n/2;
}

// reverses the order of the std::vector M
std::vector<double> Blocker::flip(std::vector<double> M)
{
    int _n = M.size();
    std::vector<double> _M(_n,0.0f);
    int i = _n-1;
    for (double& m : M) _M[i--] = m;
    return _M;
}

// performs cumulative sum on the elements of M
std::vector<double> Blocker::cumsum(std::vector<double> M)
{
    int _n = M.size();
    std::vector<double> _M(_n,0.0f);
    for (int i = 0; i < _n; i++) {
        for(int j = 0; j<=i; j++ ) _M[i] += M[j];
    }
    return _M;
}

// estimates gamma_h(1) for all h
double Blocker::gamma1(std::vector<double>& X, unsigned long long& n)
{
    double s = 0;
    for (unsigned long long i = 0; i < n-1; i++) s += X[i]*X[i+1];
    return s/n;
}

// estimates gamma_h(0) for all h
double Blocker::var(std::vector<double>& X, unsigned long long& n)
{
    double s = 0;
    for (unsigned long long i = 0; i < n; i++)  s += X[i]*X[i];
    return s/n;
}

// estimates mean of x
double Blocker::avg(std::vector<double>& x)
{
    double s = 0;
    for (double& xi : x)  s += xi;
    return s/double(x.size());
}
