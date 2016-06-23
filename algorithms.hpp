#pragma once

#include "Graph_reader.hpp"

#include <unordered_map>
#include <unordered_set>

namespace gsinfmax {
namespace algorithms {
class influence_maximization_algorithm {
  protected:
    influence_maximization_algorithm(network g) : g(g){};

    using color = int;

    network g;
    std::ranlux24_base generator;

    std::unordered_map<vertex_descriptor, color> random_propagation(
        const std::unordered_map<vertex_descriptor, color> &s,
        const std::unordered_map<vertex_descriptor, color> &ignore =
            std::unordered_map<vertex_descriptor, color>());
    std::unordered_map<vertex_descriptor, int>
    reverse_random_propagation(vertex_descriptor user);
    std::vector<vertex_descriptor>
    random_neighbors(const vertex_descriptor user);
    std::vector<vertex_descriptor>
    reverse_random_neighbors(const vertex_descriptor user);
};
}
}
