#include "Graph_reader.hpp"
#include "misc.hpp"

#include <boost/dynamic_bitset.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/home/x3.hpp>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace gsinfmax {
void set_number_of_colors(const int number_of_colors, network &g) {
    get_property(g, boost::graph_name) = number_of_colors;
};

int get_number_of_colors(const network &g) {
    return get_property(g, boost::graph_name);
}

namespace reader {
namespace x3 = boost::spirit::x3;
namespace fusion = boost::fusion;

struct event_location {
    double longitude;
    double latitude;
};

// Random number generator for the edge weights
std::ranlux24_base generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);

void generic_reader::set_random_weights() { random_weights = true; }
void generic_reader::set_no_random_weights() { random_weights = false; }
void generic_reader::set_weights(double weight) { this->weight = weight; }
void generic_reader::limit_colors(int limit) { this->color_limit = limit; }

/**
 * Return vertex_descriptor for id (potentially create new vertex_descriptor).
 *
 * If id has no associated vertex yet, create a new vertex.
 * Otherwise simply return the associated vertex.
 */
vertex_descriptor generic_reader::map_vertex(int id, network &g) {
    auto match = mapped_vertices.find(id);

    if (match != mapped_vertices.end()) {
        return match->second; // vertex already known
    } else {                  // Otherwise add vertex
        mapped_vertices[id] = boost::add_vertex(g);
        return mapped_vertices[id];
    }
};

/**
 * Return vertex_descriptor for id (don't create new vertex_descriptor).
 *
 * If id has no associated vertex yet, return null_vertex().
 * Otherwise simply return the associated vertex.
 */
vertex_descriptor generic_reader::mapped_vertex(int id) {
    auto mapped = mapped_vertices.find(id);

    return mapped == mapped_vertices.end() ? network::null_vertex()
                                           : mapped->second;
}

/**
 * Construct network by reading edges of the gowalla austin dallas dataset.
 */
network gowalla_austin_dallas::read_edges(std::string edge_file) {
    network g;

    // Lambda function that adds edge to graph
    auto add_edge = [&](auto &ctx) {
        auto attr = x3::_attr(ctx);
        auto vertex_id = fusion::at_c<0>(attr);
        auto neighbors = fusion::at_c<1>(attr);

        for (auto i : neighbors) {
            auto aer = boost::add_edge(this->map_vertex(vertex_id, g),
                                       this->map_vertex(i, g), g);

            friendship edge{};
            if (random_weights == true) {
                edge.weight = distribution(generator);
            } else {
                edge.weight = weight;
            }

            g[aer.first] = edge;
        }
    };

    // Parse the gowalla edge file
    boost::iostreams::mapped_file_source mm(edge_file);
    auto f = mm.begin(), l = mm.end();

    //               [vertex_id]        [number_of_neigbhors]
    auto n_friends = x3::int_ >> ' ' >> x3::omit[x3::int_] >> ' ' >>
                     // [neighbor]                 [1]
                     ((x3::int_ >> ' ' >> x3::omit[x3::int_]) % ' ') >> x3::eol;
    //               [vertex_id]   [0 neigbhors]
    auto zero_friends = x3::int_ >> " 0" >> x3::eol;

    x3::parse(f, l, *n_friends[add_edge]);

    // We do not need to add any edges for the users without friends
    // The users will be added to the graph when adding their location attribute
    x3::parse(f, l, *zero_friends);

    if (f != l) {
        throw std::runtime_error(
            "Could only parse " + std::to_string(std::distance(mm.begin(), f)) +
            " of " + std::to_string(std::distance(mm.begin(), l)) +
            " bytes of the edges file");
    }

    return g;
}

/**
 * Read user locations from location file and add them to the network.
 */
void gowalla_austin_dallas::read_locations(std::string location_file,
                                           network &g) {
    auto add_location = [&](auto &ctx) {
        // _attr(ctx) returns a boost fusion tuple
        auto attr = x3::_attr(ctx);
        auto vertex_id = fusion::at_c<0>(attr);

        // Add location to the user
        if (this->mapped_vertex(vertex_id) != network::null_vertex() ||
            drop_users_without_friends == false) {

            auto &props = g[this->map_vertex(vertex_id, g)];
            props.latitude = fusion::at_c<1>(attr);
            props.longitude = fusion::at_c<2>(attr);

            return true;
        } else {
            return false;
        }
    };

    boost::iostreams::mapped_file_source mm(location_file);

    auto f = mm.begin(), l = mm.end();

    x3::parse(f, l,
              // [vertex_id]        [latitude]            [longitude]
              *((x3::int_ >> ' ' >> x3::double_ >> ' ' >>
                 x3::double_)[add_location] >>
                x3::eol));

    if (f != l) {
        throw std::runtime_error(
            "Could only parse " + std::to_string(std::distance(mm.begin(), f)) +
            " of " + std::to_string(std::distance(mm.begin(), l)) +
            " bytes of the locations file");
    }

    return;
}

// Declare X3 rules to parse event location
x3::rule<class event, event_location> const event = "event";
auto const event_def = x3::double_ >> x3::double_;
BOOST_SPIRIT_DEFINE(event);

/**
 * Read event locations from file and return vector of event_location.
 */
auto generic_reader::read_events(std::string events_file) {
    // Parse the event locations
    boost::iostreams::mapped_file_source mm(events_file);
    auto f = mm.begin(), l = mm.end();

    std::vector<event_location> events;
    x3::phrase_parse(f, l, *event, x3::space, events);

    if (f != l) {
        throw std::runtime_error(
            "Could only parse " + std::to_string(std::distance(mm.begin(), f)) +
            " of " + std::to_string(std::distance(mm.begin(), l)) +
            " bytes of the events file");
    }

    if (color_limit > 0) {
        events.resize(color_limit);
    }

    return events;
}

void generic_reader::compute_user_importances(const auto events, auto &g) {
    set_number_of_colors(events.size(), g);

    BGL_FORALL_VERTICES(user, g, network) {
        g[user].importances = Eigen::ArrayXd(events.size());

        int i{0}; // Loop index
        for (auto event : events) {
            double d_double =
                misc::great_circle_length(g[user].latitude, g[user].longitude,
                                          event.latitude, event.longitude);

            // Round to next power of two
            int d = (int)d_double;
            d = std::max(1, d);
            --d;
            d |= d >> 1;
            d |= d >> 2;
            d |= d >> 4;
            d |= d >> 8;
            d |= d >> 16;
            d += 1;

            g[user].importances[i] = 1.0 / d;

            i++;
        }
    }

    return;
};

/**
 * Create network from synthetic dataset.
 */
network synthetic::read_network(std::string edge_file,
                                std::string location_file,
                                std::string events_file) {
    auto g = read_edges(edge_file);
    read_locations(location_file, g);

    auto events = read_events(events_file);

    compute_user_importances(events, g);

    return g;
}

/**
 * Construct network by reading edges of adjacency list file.
 */
network synthetic::read_edges(std::string edge_file) {
    network g;

    // Lambda function that adds edge to graph
    auto add_edge = [&](auto &ctx) {
        auto attr = x3::_attr(ctx);
        auto vertex_id = fusion::at_c<0>(attr);
        auto neighbors = fusion::at_c<1>(attr);

        for (auto i : neighbors) {
            auto aer = boost::add_edge(this->map_vertex(vertex_id, g),
                                       this->map_vertex(i, g), g);

            friendship edge{};
            if (random_weights == true) {
                edge.weight = distribution(generator);
            } else {
                edge.weight = weight;
            }

            g[aer.first] = edge;
        }
    };

    // Parse the edge file
    boost::iostreams::mapped_file_source mm(edge_file);
    auto f = mm.begin(), l = mm.end();

    //               [vertex_id]        [number_of_neigbhors]
    auto n_friends = x3::int_ >> ' ' >> x3::omit[x3::int_] >> ' ' >>
                     // [neighbor]
                     (x3::int_ % ' ') >> x3::eol;

    x3::parse(f, l, *n_friends[add_edge]);

    if (f != l) {
        throw std::runtime_error(
            "Could only parse " + std::to_string(std::distance(mm.begin(), f)) +
            " of " + std::to_string(std::distance(mm.begin(), l)) +
            " bytes of the edges file");
    }

    return g;
}

/**
 * Read user locations from location file and add them to the network.
 */
void synthetic::read_locations(std::string location_file, network &g) {
    auto add_location = [&](auto &ctx) {
        // _attr(ctx) returns a boost fusion tuple
        auto attr = x3::_attr(ctx);
        auto vertex_id = fusion::at_c<0>(attr);

        // Add location to the user
        if (this->mapped_vertex(vertex_id) != network::null_vertex()) {
            auto &props = g[this->map_vertex(vertex_id, g)];
            props.latitude = fusion::at_c<1>(attr);
            props.longitude = fusion::at_c<2>(attr);

            return true;
        } else {
            return false;
        }
    };

    boost::iostreams::mapped_file_source mm(location_file);

    auto f = mm.begin(), l = mm.end();

    x3::parse(f, l,
              // [vertex_id]        [latitude]            [longitude]
              *(x3::int_ >> ' ' >> (x3::int_ >> ' ' >> x3::double_ >> ' ' >>
                                    x3::double_)[add_location] >>
                ' ' >> x3::double_ >> ' ' >> x3::double_ >> x3::eol));

    if (f != l) {
        throw std::runtime_error(
            "Could only parse " + std::to_string(std::distance(mm.begin(), f)) +
            " of " + std::to_string(std::distance(mm.begin(), l)) +
            " bytes of the locations file");
    }

    return;
}

/**
 * Create network from gowalla austin dallas dataset.
 */
network gowalla_austin_dallas::read_network(std::string edge_file,
                                            std::string location_file,
                                            std::string events_file) {
    auto g = read_edges(edge_file);
    read_locations(location_file, g);

    auto events = read_events(events_file);

    compute_user_importances(events, g);

    return g;
}
void gowalla_austin_dallas::set_drop_users_without_friends() {
    drop_users_without_friends = true;
}
void gowalla_austin_dallas::set_no_drop_users_without_friends() {
    drop_users_without_friends = false;
}

/**
 * Construct network by reading edges of the gowalla dataset.
 */
network gowalla::read_edges(std::string edge_file) {
    network g;

    // Lambda function that adds edge to graph
    auto add_edge = [&](auto &ctx) {
        vertex_descriptor source, target;
        auto tup = std::tie(source, target);

        // Add edge from the vertices returned from context
        fusion::copy(x3::_attr(ctx), tup);
        auto aer = boost::add_edge(this->map_vertex(source, g),
                                   this->map_vertex(target, g), g);

        friendship edge{};
        if (random_weights == true) {
            edge.weight = distribution(generator);
        } else {
            edge.weight = weight;
        }

        g[aer.first] = edge;
    };

    // Parse the gowalla edge file
    boost::iostreams::mapped_file_source mm(edge_file);
    auto f = mm.begin(), l = mm.end();

    x3::parse(f, l, *((x3::int_ >> '\t' >> x3::int_ >> x3::eol)[add_edge]));

    if (f != l) {
        throw std::runtime_error(
            "Could only parse " + std::to_string(std::distance(mm.begin(), f)) +
            " of " + std::to_string(std::distance(mm.begin(), l)) +
            " bytes of the edges file");
    }

    return g;
}

void gowalla::read_locations(std::string location_file, network &g) {
    boost::dynamic_bitset<> location_already_added(num_vertices(g));

    auto add_location = [&](auto &ctx) {
        // _attr(ctx) returns a boost fusion tuple
        auto attr = x3::_attr(ctx);
        auto vertex_id = fusion::at_c<0>(attr);

        if (location_already_added.test(vertex_id)) {
            return true; // Exit, as we already stored the location for this
                         // vertex
        }
        location_already_added.set(vertex_id);

        // Test if vertex is in our graph
        // We are parsing locations from a different file than the graph,
        // so there might be inconsistencies
        if (network::null_vertex() == this->mapped_vertex(vertex_id)) {
            std::cerr << "Tried to add location to vertex " << vertex_id
                      << ", but this vertex is not in our graph" << std::endl;
            return false;
        }

        // Add location to the vertex
        auto &props = g[this->map_vertex(vertex_id, g)];
        props.latitude = fusion::at_c<1>(attr);
        props.longitude = fusion::at_c<2>(attr);

        return true;
    };

    boost::iostreams::mapped_file_source mm(location_file);

    auto f = mm.begin(), l = mm.end();

    x3::parse(f, l,
              // [vertex_id]       [time of checkin]
              *((x3::int_ >> '\t' >> x3::omit[*x3::graph] >> '\t' >>
                 // [latitude]             [longitude]
                 x3::double_ >> '\t' >> x3::double_)[add_location] >>
                // [location] id
                '\t' >> x3::int_ >> x3::eol));

    std::cout << "Added location to " << location_already_added.count()
              << " vertices" << std::endl;

    if (f != l) {
        throw std::runtime_error(
            "Could only parse " + std::to_string(std::distance(mm.begin(), f)) +
            " of " + std::to_string(std::distance(mm.begin(), l)) +
            " bytes of the locations file");
    }

    return;
}

/**
 * Create network from gowalla.
 *
 * Users will be assigned the location of their first checkin.
 * If a user has never checked in anywhere, his location is set to (0, 0)
 */
network gowalla::read_network(std::string edge_file,
                              std::string location_file) {
    auto g = read_edges(edge_file);
    read_locations(location_file, g);

    return g;
}
}
}

// Register event_location struct with boost fusion
BOOST_FUSION_ADAPT_STRUCT(gsinfmax::reader::event_location,
                          (double, latitude)(double, longitude))
