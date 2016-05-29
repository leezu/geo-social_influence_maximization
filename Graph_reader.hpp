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

	// Define our graph
	// We use setS to enforce our graph not to become a multigraph
	using network = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, user, friendship>;
	using vertex_descriptor = network::vertex_descriptor;
	using edge_descriptor = network::edge_descriptor;


	namespace reader {
		namespace gowalla_austin_dallas {
			network read_network(std::string edge_file, std::string location_file, std::string events_file);
		}

		namespace gowalla {
			network read_network(std::string edge_file, std::string location_file);

		}
	}
}
