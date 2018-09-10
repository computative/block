#include "blocker.h"

Blocker::Blocker(arma::vec& x)
{
    setData(x);
}

void Blocker::setData(arma::vec& x)
{
    unsigned long long n = arma::size(x).n_rows;
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
        mean = arma::mean(x);
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
double Blocker::estimate(arma::vec& x, unsigned long long& n)
{
    // since n is a power of two, d is an integer.
    int d = log2(n);
    arma::vec X = x - mean;
    arma::vec s = arma::zeros<arma::vec>(d);
    arma::vec gamma = arma::zeros<arma::vec>(d);
    // computes covariance and variance and applies blocking transform

    // Constructor accepts a vector
    // This vector must have a size which is a power of 2.
    for (int k = 0; k <= d-1; k++) {
        gamma[k] = this->gamma1(X, n);
        s[k] = this->var(X, n);
        this->transform(X, x, n);
    }
    // computes the test statistic Mj
    arma::vec M = arma::zeros<arma::vec>(d);
    for(int j = 0; j <= d-1; j++) {
        M[j] = pow(gamma[j]/s[j],2 )*pow(2,d-j);
    }
    M = arma::cumsum(arma::flipud(M));
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
void Blocker::transform(arma::vec& X, arma::vec& x, unsigned long long& n)
{
    for (unsigned long long i = 0; i < n/2; i++) {
        x[i] = 0.5*( x[2*i] + x[2*i+1] );
        X[i] = x[i] - mean;
    }
    n = n/2;
}

// estimates gamma_h(1) for all h
double Blocker::gamma1(arma::vec& X, unsigned long long& n)
{
    double s = 0;
    for (unsigned long long i = 0; i < n-1; i++) {
        s += X[i]*X[i+1];
    }
    return s/n;
}

// estimates gamma_h(0) for all h
double Blocker::var(arma::vec& X, unsigned long long& n)
{
    double s = 0;
    for (unsigned long long i = 0; i < n; i++) {
        s += X[i]*X[i];
    }
    return s/n;
}
