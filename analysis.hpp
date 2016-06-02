#pragma once

#include <fstream>
#include <string>

namespace gsinfmax {
std::ofstream get_logfile(std::string fname, std::string directory = "./log");
}
