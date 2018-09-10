#include "buffer.h"

Buffer::Buffer(unsigned int size)
    : size(size)
{
    load.reserve(size);
    for (unsigned int i = 0; i < size; i++) load.emplace_back( 0.0f );
    load.shrink_to_fit();
    it = load.begin();
}

void Buffer::parse(std::string& line)
{
    *it = std::strtod( line.c_str(), nullptr ) ;
    if (++it == load.end())
        is_full = true;
}

void Buffer::set_drain(const int k, std::vector<double>::iterator& start, std::vector<double>::iterator& end)
{
    this->k = k;
    this->drain_start = start;
    this->drain_end = end;
}

void Buffer::transform_and_flush()
{
    int nk = size;
    // apply k blocking transformations to load and return data
    for(int j = 0; j < k; j++)
    {
        nk = nk/2;
        for(int i = 0; i < nk; i++)
            load[i] = 0.5*( load[2*i] + load[2*i-1] );
    }
    int i = 0;
    for(std::vector<double>::iterator itr = drain_start; itr!= drain_end; itr++)
    {
        *itr = load[i];
        i++;
    }
    // reset pointer to beginning of buffer
    it = load.begin();
    is_full = false;
}
