#include "Graph_reader.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <string>

namespace x3 = boost::spirit::x3;
namespace fusion = boost::fusion;


// Define our graph
// We use setS to enforce our graph not to become a multigraph
using graph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, vertex_property, edge_property>;
using vertex_descriptor = graph::vertex_descriptor;
using edge_descriptor = graph::edge_descriptor;

bool Gowalla_austin_dallas_reader::read_edges(std::string fname) {
	// Lambda function that adds edge to graph
	auto add_edge = [&](auto& ctx){
		auto attr = x3::_attr(ctx);
		auto vertex_id = fusion::at_c<0>(attr);
		auto neighbors = fusion::at_c<1>(attr);

		for (auto i : neighbors) {
			auto aer = boost::add_edge(vertex_id, i, g);

			edge_property edge {};
			edge.weight = distribution(generator);

			g[aer.first] = edge;
		}
	};

	// Parse the gowalla edge file
	boost::iostreams::mapped_file_source mm(fname);

	auto f = mm.begin(), l = mm.end();

	//		[vertex_id]		[number_of_neigbhors]		[neighbor]		[1]
	x3::parse(f, l, *((x3::int_ >> ' ' >> x3::omit[x3::int_] >> ' ' >> ((x3::int_ >> ' ' >> x3::omit[x3::int_]) % ' '))[add_edge] >> x3::eol));

	std::cout << "Parsed " << boost::num_vertices(g) << " vertices" << std::endl;
	std::cout << "Parsed " << boost::num_edges(g) << " edges" << std::endl;

	// Fail if we couldn't parse the whole edges file
	return f == l;
}

bool Gowalla_austin_dallas_reader::read_locations(std::string fname) {
	auto add_location = [&](auto& ctx){
		// _attr(ctx) returns a boost fusion tuple
		auto attr = x3::_attr(ctx);
		auto vertex_id = fusion::at_c<0>(attr);

		// Add location to the vertex
		auto& props = g[vertex_id];
		props.latitude = fusion::at_c<1>(attr);
		props.longitude = fusion::at_c<2>(attr);

		return true;
	};

	boost::iostreams::mapped_file_source mm(fname);

	auto f = mm.begin(), l = mm.end();

	x3::parse(f, l,
			// [vertex_id]        [latitude]            [longitude]
			*((x3::int_ >> ' ' >> x3::double_ >> ' ' >> x3::double_)[add_location] >> x3::eol)
		 );

	// Fail if we couldn't parse the whole location file
	return f == l;
}

// Declare X3 rules to parse event location
BOOST_FUSION_ADAPT_STRUCT(
		event_location,
		(double, longitude)
		(double, latitude)
		)

x3::rule<class event, event_location> const event = "event";
auto const event_def = x3::double_ >> x3::double_;
BOOST_SPIRIT_DEFINE(event);

bool Gowalla_austin_dallas_reader::read_events(std::string fname) {
	// Parse the event locations
	boost::iostreams::mapped_file_source mm(fname);
	auto f = mm.begin(), l = mm.end();

	this->events = std::vector<event_location>();
	x3::phrase_parse(f, l, *event, x3::space, events);

	// Fail if we couldn't parse the whole edges file
	return f == l;
}

bool Gowalla_reader::read_edges(std::string fname) {
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

bool Gowalla_reader::read_locations(std::string fname) {
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
