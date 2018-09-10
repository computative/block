#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include <time.h>
#include <thread>
#include <unistd.h>
#include <cmath>
#include <memory>
#include <jobmanager.h>
#include <file_info.h>
#include <buffer.h>
#include <armadillo>

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


std::ostream& operator << (std::ostream& stream, File_info & file )
{
   return stream << file.name << ", " << file.id << ", " << file.row_count;
}

std::ostream& operator << (std::ostream& stream, std::vector<Job_info> & file )
{
   return stream << file[0].path;
}

std::vector<File_info> jobManager::files;
unsigned long long jobManager::d;
unsigned long long jobManager::N;
unsigned int jobManager::num_threads;
unsigned long long jobManager::row_count;
unsigned int jobManager::col_count;
std::vector<std::vector<Job_info>> jobManager::job_table;

void parsing(
        std::vector<Job_info> job_list,
        std::vector<double>::iterator start,
        std::vector<double>::iterator end
        )
{
    jobManager manager;
    unsigned int D = std::round(std::log2(jobManager::N/jobManager::num_threads));
    unsigned int b = D +1 -jobManager::d;
    unsigned int s = std::pow(2, b);
    Buffer buffer( s ); // must be a power of 2
    buffer.set_drain(D-jobManager::d, start, end); // perform D-d blocking transf to take 2^16 into
    manager.open_job(job_list); // seems to work.
    std::string line;
    int j = 0;
    while( manager.get_line(line) )
    {
        buffer.parse(line);
        if(buffer.is_full)
        {
            buffer.transform_and_flush();
        }
    }
}

int main(int argc, char *argv[])
{
    std::string dir = "/home/marius/Dokumenter/master/blocking/peter.pan/resources/test/";
    std::string file = "/home/marius/Dokumenter/master/blocking/peter.pan/resources/data.txt";
    unsigned int num_cols = 1;

    // husk å trimme av whitespace og annet som kan komme inn

    if(std::thread::hardware_concurrency() != 0)
    {
        jobManager::num_threads = std::thread::hardware_concurrency();
    }
    else {
        jobManager::num_threads = 1;
    }
    std::cout << "Using " << jobManager::num_threads << " threads. Use keyword --threads to override" << std::endl;

    try
    {
        jobManager::set_files(dir, num_cols); // sets num of cols, determines if file or path. listfiles og count rows
    }
    catch(std::string e)
    {
        std::cout << e << std::endl; return 1;
    }
    catch(...)
    {
        std::cout << "Error: There were problems accessing the supplied path." << std::endl; return 1;
    }

    std::thread thread[jobManager::num_threads-1];

    std::vector<double> bucket;
    jobManager::d = 8;
    unsigned int s = std::pow(2,jobManager::d)*jobManager::num_threads;
    bucket.reserve( s );
    for (unsigned int j = 0; j < s; j++) bucket.emplace_back(0.0f);

    // start parsing job
    /*
            */
    for (unsigned int i = 0; i < jobManager::num_threads-1; i++) thread[i] = std::thread (parsing,
            (jobManager::job_table)[i],
            bucket.begin() + i*s/jobManager::num_threads,
            bucket.begin() + (i+1)*s/jobManager::num_threads
            );
    parsing(
            (jobManager::job_table)[jobManager::num_threads-1],
            (bucket.begin() + (jobManager::num_threads-1)*s/jobManager::num_threads),
            (bucket.begin() + s)
            );
    for (unsigned int i = 0; i < jobManager::num_threads-1; i++) thread[i].join();
    /*
    */
    std::cout << "yolo" << std::endl;
    // routine to perform blocking
    return 0;
}
