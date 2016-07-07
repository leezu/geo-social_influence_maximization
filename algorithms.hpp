#pragma once

#include "Graph_reader.hpp"

#include <unordered_map>
#include <unordered_set>

namespace gsinfmax {
namespace algorithms {
struct user_color {
    vertex_descriptor user;
    int color;

    user_color(vertex_descriptor v, int c) : user(v), color(c){};
    user_color(std::pair<vertex_descriptor, int> u_c)
        : user(u_c.first), color(u_c.second){};
};

struct user_distance {
    vertex_descriptor user;
    int distance;

    user_distance(vertex_descriptor v, int d) : user(v), distance(d){};
    user_distance(std::pair<vertex_descriptor, int> u_d)
        : user(u_d.first), distance(u_d.second){};
};

class influence_maximization_algorithm {
  protected:
    influence_maximization_algorithm(network g) : g(g){};

    using color = int;

    const network g;
    std::ranlux24_base generator{std::random_device{}()};

    std::vector<user_color> random_propagation(
        const std::unordered_map<vertex_descriptor, color> &s,
        const std::unordered_map<vertex_descriptor, color> &ignore =
            std::unordered_map<vertex_descriptor, color>());
    std::vector<user_distance>
    reverse_random_propagation(vertex_descriptor user);
    std::vector<vertex_descriptor>
    random_neighbors(const vertex_descriptor user);
    std::vector<vertex_descriptor>
    reverse_random_neighbors(const vertex_descriptor user);
};
}
}
