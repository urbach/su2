// h5-dump.cpp
// reading .h5 file and dumping its content in the working directory

#include <fstream>
#include <highfive/H5Easy.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "boost/filesystem.hpp"
#include "boost/range/iterator_range.hpp"
#include <boost/algorithm/string.hpp>


int main(int argc, char **argv) {
  // check that argc == 3

  H5Easy::File file(argv[1], H5Easy::File::ReadWrite);
  std::vector<std::string> path = H5Easy::load<std::vector<std::string>>(file, "path");

  const int n = path.size();

  for (std::string p : path) {
    std::ofstream ofs(p);

    if (!boost::filesystem::exists(p)) {
      std::vector<std::string> strs;
      boost::split(strs, p, boost::is_any_of("/"));
      boost::filesystem::create_directory(strs[0]);
    }

    const std::string result = H5Easy::load<std::string>(file, p);
    ofs << result;
  }

  return 0;
}