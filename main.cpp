#include <boost/graph/adjacency_list.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include "docopt.h"

namespace x3 = boost::spirit::x3;
namespace fusion = boost::fusion;

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
using vertex_descriptor = graph::vertex_descriptor;
using edge_descriptor = graph::edge_descriptor;

class Reader {
	protected:
		graph g;
};

class Gowalla_reader : public Reader {
	public:
		bool read_edges(std::string fname) {
			// Lambda function that adds edge to graph
			auto add_edge = [&](auto& ctx){
				vertex_descriptor source, target;
				auto tup = std::tie(source, target);

				// Add edge from the vertices returned from context
				fusion::copy(x3::_attr(ctx), tup);
				auto aer = boost::add_edge(source, target, g);

				edge_property edge {};
				edge.weight = distribution(generator);

				g[aer.first] = edge;
			};

			// Parse the gowalla edge file
			boost::iostreams::mapped_file_source mm(fname);

			auto f = mm.begin(), l = mm.end();

			x3::parse(f, l, *((x3::int_ >> '\t' >> x3::int_ >> x3::eol)[add_edge]));

			std::cout << "Parsed " << boost::num_vertices(g) << " vertices" << std::endl;
			std::cout << "Parsed " << boost::num_edges(g) << " edges" << std::endl;

			// Fail if we couldn't parse the whole edges file
			return f == l;
		}

		bool read_locations(std::string fname) {
			boost::dynamic_bitset<> location_already_added(num_vertices(g));

			auto add_location = [&](auto& ctx){
				// _attr(ctx) returns a boost fusion tuple
				auto attr = x3::_attr(ctx);
				auto vertex_id = fusion::at_c<0>(attr);

				if (location_already_added.test(vertex_id)) {
					return true; // Exit, as we already stored the location for this vertex
				}
				location_already_added.set(vertex_id);

				// Test if vertex is in our graph
				// We are parsing locations from a different file than the graph,
				// so there might be inconsistencies
				if (graph::null_vertex() == vertex_id) {
					std::cerr << "Tried to add location to vertex " << vertex_id << ", but this vertex is not in our graph" << std::endl;
					return false;
				}

				// Add location to the vertex
				auto& props = g[vertex_id];
				props.latitude = fusion::at_c<1>(attr);
				props.longitude = fusion::at_c<2>(attr);

				return true;
			};

			boost::iostreams::mapped_file_source mm(fname);

			auto f = mm.begin(), l = mm.end();

			x3::parse(f, l,
					// [vertex_id]          [time of checkin]              [latitude]             [longitude]                          [location] id
					*((x3::int_ >> '\t' >> x3::omit[*x3::graph] >> '\t' >> x3::double_ >> '\t' >> x3::double_)[add_location] >> '\t' >> x3::int_ >> x3::eol)
				 );

			std::cout << "Added location to " << location_already_added.count() << " vertices" << std::endl;

			// Fail if we couldn't parse the whole location file
			return f == l;
		}

	private:
		// Random number generator for the edge weights
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution = std::uniform_real_distribution<double>(0.0, 1.0);
};

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
