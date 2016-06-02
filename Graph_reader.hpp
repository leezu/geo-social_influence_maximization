#pragma once

#include <Eigen/Core>
#include <boost/graph/adjacency_list.hpp>
#include <random>

namespace gsinfmax {
// Define classes for vertex and edge properties
class user {
  public:
    double longitude;
    double latitude;

    Eigen::ArrayXd importances; // X-dimensional vector

    user() : longitude(0), latitude(0){};
};

class friendship {
  public:
    double weight;

    friendship() : weight(0){};
};

using network = boost::adjacency_list<
    boost::setS, // This enforces the graph not to become a multigraph
    boost::vecS, boost::undirectedS, user, friendship,
    boost::property<boost::graph_name_t,
                    int> // Store the number of colors as graph_name_t
    >;
using vertex_descriptor = network::vertex_descriptor;
using edge_descriptor = network::edge_descriptor;

void set_number_of_colors(const int number_of_colors, network &g);
int get_number_of_colors(const network &g);

namespace reader {
class generic_reader {
  public:
    void set_random_weights();
    void set_no_random_weights();
    void set_weights(double weight);

  protected:
    bool random_weights{false};
    double weight{0.1};
};

class gowalla_austin_dallas : public generic_reader {
  public:
    network read_network(std::string edge_file, std::string location_file,
                         std::string events_file);

  private:
    network read_edges(std::string edge_file);
    void read_locations(std::string location_file, network &g);
    auto read_events(std::string events_file);
    void compute_user_importances(const auto events, auto &g);
};

class gowalla : public generic_reader {
  public:
    network read_network(std::string edge_file, std::string location_file);

  private:
    network read_edges(std::string edge_file);
    void read_locations(std::string location_file, network &g);
};
}
}
