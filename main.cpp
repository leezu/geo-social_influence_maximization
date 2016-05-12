#include <boost/graph/adjacency_list.hpp>
#include <boost/spirit/home/x3.hpp>
#include <iostream>

#include "docopt.h"

static const char USAGE[] =
R"(Geo-social influence maximization

    Usage:
      gsinfmax (-h | --help)
      gsinfmax --version

    Options:
      -h --help     Show this screen.
      --version     Show version.
)";

// Define classes for vertex and edge properties
struct vertex_property {
	double longitude;
	double latitude;
};

struct edge_property {
	double weight;
};

// Define our graph
// We use setS to enforce our graph not to become a multigraph
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, vertex_property, edge_property > graph;
//Some typedefs for simplicity
typedef boost::graph_traits<graph>::vertex_descriptor vertex_descr_t;
typedef boost::graph_traits<graph>::edge_descriptor edge_descr_t;

int main(int argc, char* argv[]) {
	graph g;

	return 0;
}
