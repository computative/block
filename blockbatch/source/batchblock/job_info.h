#ifndef JOBINFRO_H
#define JOBINFRO_H

#include <string>

struct Job_info {
    std::string path;
    unsigned long long pos_start;
    unsigned long long pos_end;
};

#endif // JOBINFRO_H
