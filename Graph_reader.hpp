#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <random>

// Define classes for vertex and edge properties
class vertex_property {
	public:
		double longitude;
		double latitude;

		vertex_property(): longitude(0), latitude(0) {};
};

class edge_property {
	public:
		double weight;

		edge_property(): weight(0) {};
};

class Graph_reader {
	protected:
		boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, vertex_property, edge_property> g;

	public:
		auto get_graph() {
			return g;
		}
};

class Gowalla_reader : public Graph_reader {
	public:
		bool read_edges(std::string fname);
		bool read_locations(std::string fname);

	private:
		// Random number generator for the edge weights
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution = std::uniform_real_distribution<double>(0.0, 1.0);
};
