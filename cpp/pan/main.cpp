#include <iostream>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

vec block(vec X);
double estimate(vec x);
double gamma1(vec x);
double V(vec x);
void gather(int jobsize, int nthds, int AR, int samplesize);

int main(int, char *[])
{
    vec x;
    x.load("../../resources/data.txt");
    cout << estimate(x) << endl;
    return 0;
}

double estimate(vec x)
{
    vec q {6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181, 52.191395, 53.485772, 54.775540, 56.060909, 57.342073, 58.619215, 59.892500, 61.162087, 62.428121, 63.690740};
    int n = size(x).n_rows;
    int d = log2(n);
    double mu = mean(x);
    vec s = zeros<vec>(d);
    vec gamma = zeros<vec>(d);
    for (int k = 0; k <= d-1; k++) {
        vec X = x - mu;
        gamma[k] = gamma1(X);
        s[k] = V(X);
        x = block(x);
    }
    vec observ = zeros<vec>(d);
    for(int k = 0; k <= d-1; k++) {
        observ[k] = pow(gamma[k]/s[k],2 )*pow(2,d-k);
    }
    observ = cumsum(flipud(observ));
    int k;
    for(k = d-1; k >= 0; k--) {
        if(observ[k] < q[k])
            break;
    }
    int ind = d-1-k;
    return s[ind]/pow(2, d-ind);
}

vec block(vec x)
{
    int n = size(x).n_rows;
    vec y = zeros(n/2);
    for (int i = 0; i < n/2; i++) {
        y[i] = 0.5*( x[2*i] + x[2*i+1] );
    }
    return y;
}

double gamma1(vec X)
{
    double s = 0;
    int n = size(X).n_rows;
    for (int i = 0; i < n-1; i++) {
        s += X[i]*X[i+1];
    }
    return s/(n-1);
}

double V(vec X)
{
    double s = 0;
    int n = size(X).n_rows;
    for (int i = 0; i < n; i++) {
        s += X[i]*X[i];
    }
    return s/n;
}
