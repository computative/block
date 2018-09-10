#ifndef FILEINFO_H
#define FILEINFO_H

#include <string>

struct File_info {
    std::string name;
    unsigned int id;
    unsigned int col_count;
    unsigned int row_count;
};

#endif // FILEINFO_H
