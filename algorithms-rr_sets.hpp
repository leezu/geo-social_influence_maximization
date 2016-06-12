#pragma once

#include "Graph_reader.hpp"
#include "algorithms.hpp"

namespace gsinfmax {
namespace algorithms {
class rr_sets : public influence_maximization_algorithm {
  public:
    std::unordered_map<vertex_descriptor, int>
    maximize_influence(std::vector<unsigned int> budgets);

    rr_sets(network g);

  private:
    using hypergraph = std::vector<std::vector<vertex_descriptor>>;

    std::vector<double> color_sumimportances;
    std::map<color, std::discrete_distribution<>> color_distributions;
    std::map<color, hypergraph> color_hypergraphs;
    std::map<color, hypergraph> color_user_rrsets;

    double
    estimate_influence(std::unordered_map<vertex_descriptor, color> seedset);
    double
    estimate_influence(std::unordered_map<vertex_descriptor, color> seedset,
                       color c);
    double
    estimate_influence(std::unordered_map<vertex_descriptor, color> seedset,
                       color c, int rrset_id);

    std::vector<double> estimate_lower_bound(double epsilon,
                                             std::vector<unsigned int> budgets);
    void theta_based_on_lower_bound(double epsilon,
                                    std::vector<double> color_OPT_prime,
                                    std::vector<unsigned int> budgets);
    bool assure_theta_i(const std::vector<int> &color_i,
                        const std::vector<unsigned int> &budgets,
                        double epsilon_prime);

    std::unordered_map<vertex_descriptor, color>
    build_seedset(std::vector<unsigned int> budgets);
    std::vector<std::vector<int>> get_color_user_rrsetdegrees();
    std::pair<vertex_descriptor, color> get_best_user_color(
        const std::vector<unsigned int> budgets,
        const std::vector<std::vector<int>> color_user_rrsetsdegrees);

    bool build_hypergraph_r_color(const int r, const color c);
    void add_hypergraph_edge(const vertex_descriptor user, const color c);
    vertex_descriptor random_user(const color c);

    class Math {
      public:
        static double logcnk(int n, int k);
    };
};
}
}
