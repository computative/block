#include "blocker.h"
#include <fstream>
#include <string>

int main(int, char *[])
{
    std::ifstream infile("../resources/test/data1.txt");
    std::vector<double> x;
    std::string line;
    while(std::getline(infile,line)) x.push_back(strtod(line.c_str(),nullptr) );

    // Constructor accepts a vector
    // This vector must have a size which is a power of 2.
    Blocker block(x);

    // the public variables mean, mse_mean, stdErr and mse_stdErr are output
    printf("Expected value = %g (with mean sq. err. = %g) \n", block.mean, block.mse_mean);
    printf("Standard error = %g (with mean sq. err. = %g) \n", block.stdErr, block.mse_stdErr);

    return 0;
}
