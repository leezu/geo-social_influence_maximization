#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <random>

namespace gsinfmax {
	// Define classes for vertex and edge properties
	class user {
		public:
			double longitude;
			double latitude;

			std::vector<double> importances;

			user(): longitude(0), latitude(0) {};
	};

	class friendship {
		public:
			double weight;

			friendship(): weight(0) {};
	};

	using network = boost::adjacency_list<
		boost::setS, // This enforces the graph not to become a multigraph
		boost::vecS,
		boost::undirectedS,
		user,
		friendship,
		boost::property<boost::graph_name_t, int> // Store the number of parties as graph_name_t
			>;
	using vertex_descriptor = network::vertex_descriptor;
	using edge_descriptor = network::edge_descriptor;

	void set_number_of_parties(const int number_of_parties, network& g);
	int get_number_of_parties(const network& g);

	namespace reader {
		namespace gowalla_austin_dallas {
			network read_network(std::string edge_file, std::string location_file, std::string events_file);
		}

		namespace gowalla {
			network read_network(std::string edge_file, std::string location_file);

		}
	}
}
