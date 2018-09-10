#include "jobmanager.h"
#include <filesystem>

void jobManager::set_files(std::string &path, unsigned int &col_count )
{
    namespace fs = std::filesystem;
    // set number of columns in data
    jobManager::col_count = col_count;
    // check if path is directory or file
    if( fs::is_directory((fs::path) path) )
    {
        for (auto & p : fs::directory_iterator(path))
        {
            if ( fs::is_regular_file(((fs::path) p)) )
            {
                //this->files.push_back( {p.path().string() , 0, 0, 0} );
                files.push_back( {p.path().string() , 0, 0, 0} );
            }
        }
    }
    else if ( fs::is_regular_file(((fs::path) path)) )
    {
        //this->files.push_back( {path, 0, 0, 0} );
        files.push_back( {path, 0, 0, 0} );
    }
    else
    {
        throw (std::string) "Error: A supplied path is not a directory and not a regular file.";
    }
    // count rows in data and tuncate if necessary
    set_row_count();
    set_thread_count();
}

void jobManager::set_thread_count()
{
    unsigned int thread_count = jobManager::num_threads;
    unsigned int num_chunks = std::pow(2, (int)std::ceil( std::log2(thread_count) ) );
    int leftover_chunks = num_chunks % thread_count;
    unsigned long long chunk_size = jobManager::row_count/num_chunks;
    std::vector<Job_info> _files;
    for(File_info& file : jobManager::files)
        _files.push_back({file.name,0,file.row_count});

    for (unsigned int i = 0; i < thread_count; i++)
    {
        unsigned long long rest = 0;
        std::vector<Job_info> job_list;
        leftover_chunks-- > 0 ? rest = 2*chunk_size : rest = chunk_size;
        for (Job_info& file : _files)
        {
            unsigned int diff = file.pos_end - file.pos_start;
            if (diff == 0)
            {
                continue;
            }
            else
            {
                if(rest > diff)
                {
                    job_list.push_back({file.path, file.pos_start, file.pos_end});
                    file.pos_start = file.pos_end;
                    rest -= diff;
                }
                else if(rest <= diff)
                {
                    job_list.push_back({file.path, file.pos_start, file.pos_start + rest});
                    file.pos_start = file.pos_start + rest;
                    break;
                }
            }
        }
        jobManager::job_table.push_back(job_list);
    }
}

void jobManager::seek_line(const unsigned long long pos)
{
    std::string line;
    for(unsigned long long i = 0; i < pos; i++) std::getline( filestream, line );
}

void jobManager::open_job(std::vector<Job_info> &info)
{
    job_info = info;
    filestream.open(job_info[0].path);
    seek_line(job_info[0].pos_start);
    counter = job_info[0].pos_start;
}

bool jobManager::get_line( std::string& line )
{
    if( counter == (job_info[0]).pos_end )
    {
        job_info.erase( job_info.begin() );
        filestream.close();
        if (job_info.empty()) {
            return false;
        }
        else
        {
            filestream.open( job_info[0].path );
            seek_line( job_info[0].pos_start );
            std::getline(filestream, line);
            counter = job_info[0].pos_start+1;
            return true;
        }
    }
    std::getline(filestream, line);
    counter++;
    return true;
}

void jobManager::set_row_count()
{
    if( !(jobManager::files.size() > 0) )
    {
        throw (std::string) "Error: It was impossible to process the supplied path.";
    }
    std::ifstream infile;
    std::string line;
    unsigned int n, k = 0;
    // count total data and write to jobManager::files
    for (File_info& file : files)
    {
        k = 0;
        std::cout << "Opening " << file.name << std::endl;
        infile.open(file.name);
        while( getline(infile, line) )
            k++;
        infile.close();
        file.row_count = k;
        //std::cout << k << " lines found" << std::endl;
        n += k;
    }

    // truncate row_count to base 2
    jobManager::row_count = std::pow( 2, std::floor(std::log2(double(n)/double(jobManager::num_threads)) ) );
    jobManager::N = jobManager::row_count;
    k = 0; bool TRUC_FLAG = false;
    for (File_info& file : files)
    {
        if ( k + file.row_count > jobManager::row_count )
        {
            file.row_count = jobManager::row_count - k;
            TRUC_FLAG = true;
        }
        k += file.row_count;
    }
    // if data was truncated, print warning.
    if (TRUC_FLAG)
    {
        std::cout << "Notice: observation count must be a power of 2. Some data is ignored." << std::endl;
    }

}
