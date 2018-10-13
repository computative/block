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
#include <blocker.h>


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

    jobManager::d = 8;
    unsigned int s = jobManager::num_threads*std::pow(2,jobManager::d);
    //bucket.reserve( s );
    std::vector<double> bucket(s,0.0f);
    //for (unsigned int j = 0; j < s; j++) bucket.emplace_back(0.0f);

    // start parsing job

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
    Blocker blocker(bucket);
    printf("Expected value = %g (with mean square error %g)\n", blocker.mean, blocker.mse_mean );
    printf("Standard error = %g (with mean square error %g)\n", blocker.stdErr, blocker.mse_stdErr );
    return 0;
}
