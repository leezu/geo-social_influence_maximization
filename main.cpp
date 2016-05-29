#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <fstream>
#include <string>

#include "docopt.h"

#include "Graph_reader.hpp"

static const char USAGE[] =
R"(Geo-social influence maximization

    Usage:
      gsinfmax gowalla <edges> <locations>
      gsinfmax gowalla_austin_dallas <edges> <locations> <events>
      gsinfmax (-h | --help)
      gsinfmax --version

    Options:
      -h --help     Show this screen.
      --version     Show version.
)";

using namespace gsinfmax;

network get_network(auto args) {
	if (args.at("gowalla").asBool()) {
		return reader::gowalla::read_network(args["<edges>"].asString(),
				args["<locations>"].asString());
	} else if (args.at("gowalla_austin_dallas").asBool()) {
		return reader::gowalla_austin_dallas::read_network(args["<edges>"].asString(),
				args["<locations>"].asString(),
				args["<events>"].asString());
	} else {
		throw std::runtime_error("Couldn't get network for given argument");
	}
}

int main(int argc, char* argv[]) {
	std::map<std::string, docopt::value> args
		= docopt::docopt(USAGE,
				{ argv + 1, argv + argc },
				true,						// show help if requested
				"geo-social influence maximzation 0.1");	// version string

	network g = get_network(args);

	return 0;
}
