#ifndef BLOCKER_H
#define BLOCKER_H

#include <iostream>
#include <vector>


class Blocker
{
public:
    unsigned long long nK;
    // variables that hold the results
    double stdErr, mean, mse_stdErr, mse_mean;
    // a vector of magic numbers
    std::vector<double> quantile {
        3.841459,  5.991465,  7.814728,  9.487729,  11.070498, 12.591587, 14.067140, 15.507313,
        16.918978, 18.307038, 19.675138, 21.026070, 22.362032, 23.684791, 24.995790, 26.296228,
        27.587112, 28.869299, 30.143527, 31.410433, 32.670573, 33.924438, 35.172462, 36.415029,
        37.652484, 38.885139, 40.113272, 41.337138, 42.556968, 43.772972, 44.985343, 46.194260,
        47.399884, 48.602367, 49.801850, 50.998460, 52.192320, 53.383541, 54.572228, 55.758479,
        56.942387, 58.124038, 59.303512, 60.480887, 61.656233, 62.829620, 64.001112, 65.170769,
        66.338649, 67.504807, 68.669294, 69.832160, 70.993453, 72.153216, 73.311493, 74.468324,
        75.623748, 76.777803, 77.930524, 79.081944, 80.232098, 81.381015, 82.528727, 83.675261};

    Blocker(std::vector<double>& x);

    void setData(std::vector<double>& x);

    // checks that integer n is a power of 2
    bool isPowerOfTwo(unsigned long long n);

private:

    // the algorithm which computes the variance of the sample mean.
    double estimate(std::vector<double>& x, unsigned long long& n);

    // performs blocking transformation
    void transform(std::vector<double>& X, std::vector<double>& x, unsigned long long& n);

    // estimates gamma_h(1) for all h
    double gamma1(std::vector<double>& X, unsigned long long& n);

    // estimates gamma_h(0) for all h
    double var(std::vector<double>& X, unsigned long long& n);

    std::vector<double> cumsum(std::vector<double> M);

    std::vector<double> flip(std::vector<double> M);

    double avg(std::vector<double>& x);
};

#endif // BLOCKER_H
