#pragma once

#include "Graph_reader.hpp"

#include <Eigen/Core>
#include <fstream>
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

    void enable_generate_statistics();
    void disable_generate_statistics();

    lazy_greedy(network g) : g(g){};
    lazy_greedy(network g, int mc_sim) : g(g), number_of_mc_sim(mc_sim){};

  private:
    int number_of_mc_sim{10'000};
    int number_of_parties;
    std::default_random_engine generator;
    network g;

    bool statusline_printed{false};
    bool generate_statistics{false};
    int current_iteration{
        0}; ///< Stores iteration to facilitate generating statistics
    std::ofstream logfile;

    Eigen::ArrayXd marginal_influence_gain(
        const vertex_descriptor u, std::unordered_map<vertex_descriptor, int> s,
        Eigen::ArrayXd sigma_s,
        const std::unordered_map<vertex_descriptor, int> ignore =
            std::unordered_map<vertex_descriptor, int>());
    std::unordered_map<vertex_descriptor, int> random_propagation(
        const std::unordered_map<vertex_descriptor, int> &s,
        const std::unordered_map<vertex_descriptor, int> &ignore =
            std::unordered_map<vertex_descriptor, int>());
    std::unordered_set<vertex_descriptor>
    random_neighbors(const vertex_descriptor user);
    Eigen::ArrayXd global_importance_of_user_set(const auto &set);
    Eigen::ArrayXd importance_of_user_set(const auto &set);

    void
    update_statusline_seeds(const int current_seedset_size,
                            const std::vector<unsigned int> &budgets,
                            const int number_of_mgs_updates_current_iteration,
                            const int total_number_of_mgs_updates);
    Eigen::ArrayXd
    influence(const std::unordered_map<vertex_descriptor, int> &seedset,
              const std::unordered_map<vertex_descriptor, int> &ignore =
                  std::unordered_map<vertex_descriptor, int>());
};
}
}
