#include <boost/graph/adjacency_list.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include "docopt.h"

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

// Parse gowalla dataset
void parse_gowalla(auto g, std::string edge_path, std::string location_path) {
	namespace x3 = boost::spirit::x3;
	namespace fusion = boost::fusion;

	std::unordered_map<int, vertex_descr_t>  vertexes;

	// Lambda function that adds vertex to graph if not already added
	auto add_vertex = [&](auto& ctx){
		// Return if the vertex is already known
		if (vertexes.find(x3::_attr(ctx)) != vertexes.end())	{
			return false;
		}

		// Otherwise add vertex to graph
		auto v = boost::add_vertex(g);
		// TODO: Add the location to the vertex

		// And add vertex descriptor to map
		vertexes[x3::_attr(ctx)] = v;
	};

	// Lambda function that adds edge to graph
	auto add_edge = [&](auto& ctx){
		// _attr(ctx) returns a boost fusion tuple
		auto attr = x3::_attr(ctx);

		// Add edge from the vertexes returned from context
		boost::add_edge(vertexes[fusion::at_c<0>(attr)],
				vertexes[fusion::at_c<1>(attr)], g);
	};

	std::ifstream edge_file(edge_path);
	if (!edge_file) {
		std::cerr << "Couldn't open " << edge_path << " for reading";
	}
	// By default whitespace is skipped, so we disable that
	edge_file >> std::noskipws;
	// TODO: open location file

	// Parse the gowalla edge file
	boost::spirit::istream_iterator edge_file_iterator_first(edge_file), eof;

	x3::phrase_parse(edge_file_iterator_first, eof,
			// Begin grammar
			(
			 *((x3::int_[add_vertex] >> x3::int_[add_vertex])[add_edge])
			),
			// End grammar
			x3::space
			);

	// Fail if we couldn't parse the whole file
	if (edge_file_iterator_first != eof) {
		std::cerr << "Couldn't parse whole gowalla file" << std::endl;
	}

	std::cout << "Parsed " << boost::num_vertices(g) << " vertices" << std::endl;
	std::cout << "Parsed " << boost::num_edges(g) << " edges" << std::endl;

	return;
}

int main(int argc, char* argv[]) {
	std::map<std::string, docopt::value> args
		= docopt::docopt(USAGE,
				{ argv + 1, argv + argc },
				true,						// show help if requested
				"geo-social influence maximzation 0.1");	// version string

	graph g;

	// Gowalla dataset
	if (args.find("gowalla") != args.end()) {
		parse_gowalla(g, args["<edges>"].asString(), args["<locations>"].asString());
	}

	return 0;
}
