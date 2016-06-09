#pragma once

#include "Graph_reader.hpp"

#include <fstream>
#include <string>

namespace gsinfmax {
std::ofstream get_logfile(std::string fname, std::string directory = "./log");

double get_average_node_degree(const network &g);
double write_node_degrees(const network &g,
                          const std::string fname = "node_degrees");
}
