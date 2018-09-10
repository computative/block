#include <iostream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

int main()
{
  std::string inpath = "/home/marius/Dokumenter/master/blocking/peter.pan/resources/test";
  for (auto & p : fs::directory_iterator(inpath))
      std::cout << p << std::endl;
    return 0;
}
