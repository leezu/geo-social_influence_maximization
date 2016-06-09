#include "analysis.hpp"

#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
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

double get_average_node_degree(const network &g) {
    double node_degree{0};

    BGL_FORALL_VERTICES(user, g, network) { node_degree += degree(user, g); }

    return node_degree / num_vertices(g);
}
double write_node_degrees(const network &g, const std::string fname) {
    auto f = get_logfile(fname);
    BGL_FORALL_VERTICES(user, g, network) {
        f << user << "\t" << degree(user, g) << "\n";
    }
}
}
