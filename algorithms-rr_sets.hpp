#pragma once

#include "Graph_reader.hpp"
#include "algorithms.hpp"

namespace gsinfmax {
namespace algorithms {
class ris : public influence_maximization_algorithm {
  public:
    std::pair<std::vector<double>, std::unordered_map<vertex_descriptor, int>>
    maximize_influence(std::vector<unsigned int> budgets, double epsilon = 0.1,
                       double delta = 0.1);
    ris(network g);

  private:
    using user_rr_set_ids = std::vector<std::vector<int>>;
    using color_user_rr_set_ids = std::vector<user_rr_set_ids>;
    using rr_set = std::vector<user_distance>;
    using rr_sets = std::vector<std::vector<user_distance>>;
    using color_rr_sets = std::vector<rr_sets>;

    std::vector<double> color_sumimportances;
    std::vector<std::discrete_distribution<>> color_distributions;

    double estimate_influence(
        const std::unordered_map<vertex_descriptor, color> &seedset,
        const color_rr_sets &color_rs, const color_user_rr_set_ids &color_ids);
    double estimate_influence(
        const std::unordered_map<vertex_descriptor, color> &seedset, color c,
        const rr_sets &rs, const user_rr_set_ids &ids);
    double estimate_influence(
        const std::unordered_map<vertex_descriptor, color> &seedset, color c,
        const rr_set &r);
    double coverage(const std::unordered_map<vertex_descriptor, color> &seedset,
                    const color_rr_sets &color_rs,
                    const color_user_rr_set_ids &color_ids);
    double coverage(const std::unordered_map<vertex_descriptor, color> &seedset,
                    color c, const rr_sets &rs, const user_rr_set_ids &ids);

    bool stopping_criterion(std::vector<unsigned int> budgets, double epsilon,
                            double delta, const color_rr_sets &color_rs);

    std::pair<std::vector<rr_sets>, std::vector<user_rr_set_ids>>
    generate_rrsets(const int n);
    std::pair<rr_sets, user_rr_set_ids> generate_rrsets(const int n,
                                                        const color c);

    std::pair<color_rr_sets, color_user_rr_set_ids>
    merge_rr_sets(const color_rr_sets &color_rs,
                  const color_rr_sets &tmp_color_rs);

    std::unordered_map<vertex_descriptor, color>
    build_seedset(const std::vector<unsigned int> &budgets,
                  const color_rr_sets &color_r,
                  const color_user_rr_set_ids &color_ids);
    std::vector<std::vector<double>>
    get_color_user_rrsetdegrees(const color_rr_sets &color_rs,
                                const color_user_rr_set_ids &color_ids);
    std::pair<vertex_descriptor, color> get_best_user_color(
        const std::vector<unsigned int> &budgets,
        const std::vector<std::vector<double>> &color_user_rrsetsdegrees,
        const color_rr_sets &color_rs);

    vertex_descriptor random_user(const color c);

    class Math {
      public:
        static double logcnk(int n, int k);
    };
};
}
}
