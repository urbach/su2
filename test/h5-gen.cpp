// h5-gen.cpp
// generate .h5 file from directory path

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "boost/filesystem.hpp"
#include "boost/range/iterator_range.hpp"

#include <highfive/H5Easy.hpp>

std::vector<std::string> ls_recursive(const std::string &path) {
  boost::filesystem::path dir = path;
  boost::filesystem::recursive_directory_iterator it(dir), end;

  std::vector<std::string> files;
  for (auto &entry : boost::make_iterator_range(it, end)) {
    if (boost::filesystem::is_regular(entry)) {
      files.push_back(entry.path().native());
    }
  }

  return files;
}

int main(int argc, char **argv) {
  // check that argc == 3

  // list all files in directory argv[1]

  std::string h5_file = argv[2];
  std::vector<std::string> path; // vector of paths

  if (boost::filesystem::exists(h5_file)) {
    std::cerr << "Error. File already exists: " << h5_file << "\n";
    std::abort();
  }

  H5Easy::File file2(h5_file, H5Easy::File::Overwrite);

  const std::vector<std::string> files = ls_recursive(argv[1]);
  for (std::string f : files) {
    path.push_back(f);

    std::string cnt;
    boost::filesystem::load_string_file(f, cnt);
    H5Easy::dump(file2, f, cnt); // dumping vector of paths

  }

  H5Easy::dump(file2, "path", path); // dumping vector of paths
 
  return 0;
}