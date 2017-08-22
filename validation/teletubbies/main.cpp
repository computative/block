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

vec AR1(int n, vec phi);
vec ar1acf(vec h, double phi);
vec ar2acf(vec h, double phi1, double phi2);
vec AR2(int n, vec phi);
vec block(vec X);
vec estimates(vec x);
vec estimatesv2(vec x);
double exact(int n, vec phis);
double gamma1(vec x);
double V(vec x);
void gather(int jobsize, int nthds, int AR, int samplesize);
void generate(int AR, int jobsize, int samplesize);

int main(int argc, char *argv[])
{

    int AR = atoi(argv[1]);
    int samplesize = pow(2,atoi(argv[2]));
    int jobsize = atoi(argv[3]);
    int rank, nthds;


    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &nthds);

    if (rank == 0) {
        gather(jobsize, nthds, AR, samplesize);
    } else {
        generate(AR, jobsize, samplesize);
    }

    MPI_Finalize();
    return 0;
}

void generate(int AR, int jobsize, int samplesize) {
    double * result;
    result = new double[5];
    for (int l = 0; l<jobsize; l++) {
        double post = floor(log10(samplesize));
        double pre = pow(10, (log10(samplesize) - post) );
        vec phi, x;
        // AR1 eller AR2
        if (AR == 2) {
            double phi2 = -exp( - pow(10, param(gen)-post ) /pre );
            phi = {rand_float(gen)*sqrt(-4*phi2),phi2};
            x = AR2(samplesize,phi);
        } else {
            phi = { exp( - pow(10, param(gen)-post ) /pre ),0};
            x = AR1(samplesize,phi);
        }

        // beregne eksakt
        double target = exact(samplesize,phi);

        // blocking estimater
        vec blocking = estimatesv2(x);

        result[0] = target; result[1] = blocking[0];
        result[2] = blocking[1]; result[3] = phi[0]; result[4] = phi[1];

        MPI_Send(&result[0], 5, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
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

double exact(int n, vec phi)
{
    vec h = linspace(0,n-1,n);
    if (!phi[1]) {
        h = ar1acf(h,phi[0]);
    } else {
        h = ar2acf(h,phi[0],phi[1] );
    }
    double s = h[0];
    for (int i = 1; i < n; i++) {
        s += 2*(1 - (double)i/n)*h[i];
    }
    return s/n;
}

vec AR1(int n, vec phi)
{
    vec data = zeros<vec>(n);
    data[0] = epsilon(gen) - M_MU;
    for (int i = 1; i < n; i++)
        data[i] = phi[0]*data[i-1] + epsilon(gen) - M_MU;
    return data;
}

vec AR2(int n, vec phi)
{
    vec data = zeros<vec>(n);
    data[0] = epsilon(gen) - M_MU;
    data[1] = phi[0]*data[0] + epsilon(gen) - M_MU;
    for (int i = 2; i < n; i++)
        data[i] = phi[0]*data[i-1] + phi[1]*data[i-2] + epsilon(gen)-M_MU;
    return data;
}

vec estimatesv2(vec x)
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
/*
vec estimates(vec x)
{
    int n = size(x).n_rows;
    int d = log2(n);
    double mu = mean(x);
    double S; // old blocking method
    double s = V(x - mu); // new method
    for (int k = 0; k <= d-1; k++) {
        int nk = size(x).n_rows;
        vec X = x - mu;
        double gammai = gamma1(X);
        double Vx = V(X);
        S = Vx/nk;
        x = block(x);
        if ( abs(gammai) < 1.96*Vx/sqrt(nk) )
            break;
        s += pow(2,k)*gammai;
    }
    vec results {s/(n-1), S};
    return results;
}
*/

vec ar2acf(vec h, double phi1, double phi2)
{
    int n = size(h).n_rows;

    // arg and norm of complex num
    double x = -0.5*phi1/phi2; double y = -0.5*sqrt(-phi1*phi1 - 4*phi2)/phi2;
    double z = 1/sqrt(-phi2);
    double theta = atan2(y,x);

    // compute acf
    double b = atan( (phi1/(z*(1-phi2)) - cos(theta) )/sin(theta) ) ; double a = 1.0/cos(b);
    double sigma2 = (1 - phi2)/( (1 + phi2)*( (1-phi2)*(1-phi2) - phi1*phi1 ) );
    for (int i = 0; i<n; i++)
        h[i] != 0.0 ? h[i] = a*pow(z,-h[i])*cos(h[i]*theta + b) : h[i] = 1.0;
    return sigma2*h;
}

vec ar1acf(vec h, double phi)
{
    int n = size(h).n_rows;
    for (int i = 0; i<n; i++) {
        h[i] = pow(phi,abs(h[i]) );
    }
    return h/(1 - phi*phi);
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
