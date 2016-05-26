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
      gsinfmax (-h | --help)
      gsinfmax --version

    Options:
      -h --help     Show this screen.
      --version     Show version.
)";

// Define our graph
// We use setS to enforce our graph not to become a multigraph
using graph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, vertex_property, edge_property>;

int main(int argc, char* argv[]) {
	std::map<std::string, docopt::value> args
		= docopt::docopt(USAGE,
				{ argv + 1, argv + argc },
				true,						// show help if requested
				"geo-social influence maximzation 0.1");	// version string


	// Gowalla dataset
	if (args.find("gowalla") != args.end()) {
		auto reader = Gowalla_reader();

		reader.read_edges(args["<edges>"].asString());
		reader.read_locations(args["<locations>"].asString());
	}

	return 0;
}
