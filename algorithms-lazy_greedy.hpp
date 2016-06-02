#pragma once

#include "Graph_reader.hpp"

#include <Eigen/Core>
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace gsinfmax {
namespace algorithms {
class lazy_greedy {
  public:
    std::unordered_map<vertex_descriptor, int>
    maximize_influence(std::vector<unsigned int> budgets);
    std::unordered_map<vertex_descriptor, int>
    maximize_influence_baseline(std::vector<unsigned int> budgets);

    lazy_greedy(network g) : g(g){};
    lazy_greedy(network g, int mc_sim) : g(g), number_of_mc_sim(mc_sim){};

  private:
    int number_of_mc_sim{10'000};
    int number_of_parties;
    std::default_random_engine generator;
    network g;

    Eigen::ArrayXd marginal_influence_gain(
        const vertex_descriptor u,
        std::unordered_map<vertex_descriptor, int> s =
            std::unordered_map<vertex_descriptor, int>(),
        const std::unordered_map<vertex_descriptor, int> ignore =
            std::unordered_map<vertex_descriptor, int>());
    std::unordered_map<vertex_descriptor, int> random_propagation(
        const std::unordered_map<vertex_descriptor, int> &s,
        const std::unordered_map<vertex_descriptor, int> &ignore =
            std::unordered_map<vertex_descriptor, int>());
    std::unordered_set<vertex_descriptor>
    random_neighbors(const vertex_descriptor user);
    Eigen::ArrayXd importance_of_user_set(const auto &set);
};
}
}
