#pragma once

#include "Graph_reader.hpp"
#include "algorithms.hpp"

#include <Eigen/Core>
#include <fstream>
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace gsinfmax {
namespace algorithms {
class lazy_greedy : public influence_maximization_algorithm {
  public:
    std::unordered_map<vertex_descriptor, int>
    maximize_influence(std::vector<unsigned int> budgets);
    std::unordered_map<vertex_descriptor, int>
    maximize_influence_baseline(std::vector<unsigned int> budgets);

    void enable_generate_statistics();
    void disable_generate_statistics();

    void enable_statusline();
    void disable_statusline();

    lazy_greedy(network g) : influence_maximization_algorithm(g){};
    lazy_greedy(network g, int mc_sim)
        : influence_maximization_algorithm(g), number_of_mc_sim(mc_sim){};

  private:
    class logger {
      public:
        void sigma_set(const Eigen::ArrayXd &sigma_set, const int &iteration,
                       const std::string &fname);

      private:
        std::unordered_map<std::string, std::ofstream> files;
        std::unordered_set<std::string> initialized;
    };

    int number_of_mc_sim{10'000};
    int number_of_parties;

    bool statusline_enabled{true};
    bool statusline_printed{false};
    bool generate_statistics{false};
    logger log;

    Eigen::ArrayXd marginal_influence_gain(
        const vertex_descriptor u, std::unordered_map<vertex_descriptor, int> s,
        Eigen::ArrayXd sigma_s,
        const std::unordered_map<vertex_descriptor, int> ignore =
            std::unordered_map<vertex_descriptor, int>());
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
