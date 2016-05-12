#include <boost/graph/adjacency_list.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

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

	std::unordered_map<int, vertex_descr_t>  vertices; // Store vertex descriptors
	std::unordered_set<int>  location_already_added; // Store vertices, for which we already added location

	// Lambda function that adds vertex to graph if not already added
	auto add_vertex = [&](auto& ctx){
		// Return if the vertex is already known
		if (vertices.find(x3::_attr(ctx)) != vertices.end())	{
			return false;
		}

		// Otherwise add vertex to graph
		auto v = boost::add_vertex(g);

		// And add vertex descriptor to map
		vertices[x3::_attr(ctx)] = v;
	};

	// Lambda function that adds edge to graph
	auto add_edge = [&](auto& ctx){
		// _attr(ctx) returns a boost fusion tuple
		auto attr = x3::_attr(ctx);

		// Add edge from the vertices returned from context
		boost::add_edge(vertices[fusion::at_c<0>(attr)],
				vertices[fusion::at_c<1>(attr)], g);
	};

	// Lambda function that adds locations to vertices in the graph
	auto add_location = [&](auto& ctx){
		// _attr(ctx) returns a boost fusion tuple
		auto attr = x3::_attr(ctx);
		auto vertex_id = fusion::at_c<0>(attr);

		if (location_already_added.find(vertex_id) != location_already_added.end())	{
			// Exit, as we already stored the location for this vertex
			return true;
		}
		location_already_added.insert(vertex_id);

		// Test if vertex is in our graph
		// We are parsing locations from a different file than the graph,
		// so there might be inconsistencies
		if (vertices.find(vertex_id) == vertices.end())	{
			std::cerr << "Tried to add location to vertex " << vertex_id << ", but this vertex is not in our graph" << std::endl;
			return false;
		}

		auto vertex = vertices[vertex_id];

		// Add location to the vertex
		g[vertex].latitude = fusion::at_c<2>(attr);
		g[vertex].longitude = fusion::at_c<3>(attr);

		return true;
	};

	std::ifstream edge_file(edge_path);
	std::ifstream location_file(location_path);
	if (!edge_file) {
		std::cerr << "Couldn't open " << edge_path << " for reading";
	}
	if (!location_file) {
		std::cerr << "Couldn't open " << location_path << " for reading";
	}
	// By default whitespace is skipped, so we disable that
	edge_file >> std::noskipws;
	location_file >> std::noskipws;

	// Parse the gowalla edge file
	boost::spirit::istream_iterator file_iterator(edge_file), eof;

	x3::phrase_parse(file_iterator, eof,
			// Begin grammar
			(
			 *((x3::int_[add_vertex] >> x3::int_[add_vertex])[add_edge])
			),
			// End grammar
			x3::space
			);

	// Fail if we couldn't parse the whole edges file
	if (file_iterator != eof) {
		std::cerr << "Couldn't parse whole edges file" << std::endl;
	}

	// Parse the gowalla location file
	file_iterator = boost::spirit::istream_iterator(location_file);

	x3::phrase_parse(file_iterator, eof,
			// Begin grammar
			(
			 // vertex_id 	time of checkin 	  latitude 	longitude 		      location id
			 *((x3::int_ >> x3::lexeme[*x3::graph] >> x3::double_ >> x3::double_)[add_location] >> x3::int_ >> x3::eol)
			),
			// End grammar
			x3::blank
			);

	// Fail if we couldn't parse the whole location file
	if (file_iterator != eof) {
		std::cerr << "Couldn't parse whole location file" << std::endl;
	}

	std::cout << "Parsed " << boost::num_vertices(g) << " vertices" << std::endl;
	std::cout << "Parsed " << boost::num_edges(g) << " edges" << std::endl;
	std::cout << "Added location to " << location_already_added.size() << " vertices" << std::endl;

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
