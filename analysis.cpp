#include "analysis.hpp"

#include "boost/filesystem.hpp"
#include <fstream>
#include <iostream>
#include <string>

namespace gsinfmax {
std::ofstream get_logfile(std::string fname, std::string directory) {
    boost::filesystem::path rootPath(directory);
    boost::system::error_code returnedError;

    boost::filesystem::create_directories(rootPath, returnedError);

    if (returnedError) {
        std::cout << "Did not successfully create directory " << directory
                  << std::endl;
    }

    return std::ofstream(directory + "/" + fname);
}
}
