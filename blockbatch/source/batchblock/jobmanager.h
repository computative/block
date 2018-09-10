#ifndef JOBMANAGER_H
#define JOBMANAGER_H

#include <iostream>
#include <file_info.h>
#include <job_info.h>
#include <vector>
#include <fstream>
#include <cmath>

class jobManager {

public:
    static unsigned long long row_count;
    static unsigned int col_count;
    static std::vector<File_info> files;
    static std::vector<std::vector<Job_info>> job_table;
    static unsigned long long d;
    static unsigned long long N;
    static unsigned int num_threads;

    unsigned long long counter;
    std::ifstream filestream;
    std::vector<Job_info> job_info;

    static void set_files(std::string& path, unsigned int& col_count);
    static void set_thread_count();
    bool get_line(std::string& line);
    void open_job(std::vector<Job_info>& info);

private:
    static void set_row_count();
    void seek_line( const unsigned long long pos );
};

#endif // JOBMANAGER_H
