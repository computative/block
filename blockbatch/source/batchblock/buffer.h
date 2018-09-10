#ifndef BUFFER_H
#define BUFFER_H

#include <vector>
#include <string>

class Buffer
{
private:
    std::vector<double>::iterator it;
    std::vector<double>::iterator drain_start;
    std::vector<double>::iterator drain_end;
    int k;
public:
    bool is_full;
    unsigned int size;
    std::vector<double> load;
    int files;

    Buffer(unsigned int size);
    void parse(std::string& line);
    void set_drain(const int k, std::vector<double>::iterator& start, std::vector<double>::iterator& end);
    void transform_and_flush();
};

#endif // BUFFER_H
