#include <iostream>
#include <fstream>
#include <armadillo>
#include <mpi.h>

using namespace std;
using namespace arma;

random_device rd; mt19937 gen(rd());
//normal_distribution<double> epsilon(0, 1); const double M_MU = 0.0;
gamma_distribution<double> epsilon(1.0,1.0); const double M_MU = 1.0;
uniform_real_distribution<double> param(2, 6 );
uniform_real_distribution<double> rand_float(0, 1);
uniform_int_distribution<int> rand_int(0, 100000);


vec block(vec X);
vec estimates(vec x);
double gamma1(vec x);
double V(vec x);
void gather(int jobsize, int nthds, int AR, int samplesize);


int main(int argc, char *argv[])
{

}



void gather(int jobsize, int nthds, int AR, int samplesize)
{
    double * result;
    result = new double[5];
    // samle inn data
    ofstream outfile;

    outfile.open("AR" + to_string(AR) + "Log2n" + to_string( (int) log2(samplesize)) + "_" + to_string(rand_int(gen)) + ".txt");
    for (int u = 0; u < jobsize; u++) {
        for (int k = 1; k<nthds; k++) {
            MPI_Recv(&result[0], 5, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            outfile << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << " " << result[4] << endl;
        }
    }
    outfile.close();
}

vec estimates(vec x)
{
    //vec q {3.841459,  5.991465,  7.814728,  9.487729, 11.070498, 12.591587, 14.067140, 15.507313, 16.918978, 18.307038, 19.675138, 21.026070, 22.362032, 23.684791, 24.995790, 26.296228, 27.587112, 28.869299, 30.143527, 31.410433, 32.670573, 33.924438, 35.172462, 36.415029, 37.652484, 38.885139, 40.113272, 41.337138, 42.556968, 43.772972}; // 95%
    vec q {6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181}; // 99%
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
    vec correction = zeros<vec>(d);
    double c = 0.85;
    for(int k = 0; k <= d-1; k++) {
        observ[k] = pow(gamma[k]/s[k],2 )*pow(2,d-k);
        correction[k] = pow(c,k)*gamma[k];
    }
    observ = cumsum(flipud(observ));
    int k;
    for(k = d-1; k >= 0; k--) {
        if(observ[k] < q[k])
            break;
    }
    int ind = d-1-k;
    double nk = pow(2, d-ind);
    vec results {(s[ind] + sum( correction(span(ind,d-1) ) ) )/nk, s[ind]/nk};
    return results;
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
