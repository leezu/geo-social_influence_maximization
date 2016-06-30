#include "algorithms-rr_sets.hpp"

#include "Graph_reader.hpp"

#include <algorithm>
#include <boost/graph/iteration_macros.hpp>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <tuple>

namespace gsinfmax {
namespace algorithms {
namespace {
using user_rr_set_ids = std::vector<std::vector<int>>;
using color_user_rr_set_ids = std::vector<user_rr_set_ids>;
using rr_set = std::vector<user_distance>;
using rr_sets = std::vector<std::vector<user_distance>>;
using color_rr_sets = std::vector<rr_sets>;
using color = int;
}
std::pair<std::vector<double>, std::unordered_map<vertex_descriptor, int>>
ris::maximize_influence(std::vector<unsigned int> budgets, double epsilon,
                        double delta) {
    double Lambda = 2.0 * (2 * (std::exp(1)) - 2) * (1 + epsilon) *
                    (1 + epsilon) *
                    std::log(2.0 / delta * 1 / (epsilon * epsilon));

    color_rr_sets color_rs;
    color_user_rr_set_ids color_ids;
    std::tie(color_rs, color_ids) = generate_rrsets(Lambda);

    auto seedset = build_seedset(budgets, color_rs, color_ids);
    auto cov = estimate_influence(seedset, color_rs, color_ids);

    while (!stopping_criterion(budgets, epsilon, delta, color_rs)) {
        color_rr_sets tmp_color_rs;
        color_user_rr_set_ids tmp_color_ids;
        std::tie(tmp_color_rs, tmp_color_ids) =
            generate_rrsets(color_rs[0].size());

        auto tmp_cov = estimate_influence(seedset, tmp_color_rs, tmp_color_ids);

        std::cout << "influence " << cov << " influence(tmp) " << tmp_cov
                  << std::endl;

        double epsilon_1 = cov / tmp_cov - 1;
        if (epsilon_1 <= epsilon) {
            double epsilon_2 = (epsilon - epsilon_1) / (2 * (1 + epsilon_1));
            double epsilon_3 = (epsilon - epsilon_1) / (2 * (1 - epsilon_1));

            double delta_1 = std::exp(
                (-1.0) * (coverage(seedset, color_rs, color_ids) * epsilon_3 *
                          epsilon_3 / (2.0 * 2.0 * (std::exp(1) - 2.0) *
                                       (1 + epsilon_1) * (1 + epsilon_2))));
            double delta_2 = std::exp(
                (-1.0) * (coverage(seedset, tmp_color_rs, tmp_color_ids) *
                          epsilon_2 * epsilon_2 /
                          (2.0 * 2.0 * (std::exp(1) - 2.0) * (1 + epsilon_2))));

            std::cout << "coverage " << coverage(seedset, color_rs, color_ids)
                      << " coverage(tmp) "
                      << coverage(seedset, tmp_color_rs, tmp_color_ids)
                      << std::endl;

            std::cout << "delta_1 + delta_2 <= delta: " << delta_1 << " + "
                      << delta_2 << " <= " << delta << std::endl;

            if (delta_1 + delta_2 <= delta) {
                std::cout << "Returning seedset early" << std::endl;

                std::vector<double> influences(color_rs.size());
                for (color c{0}; c < color_rs.size(); c++) {
                    influences[c] = estimate_influence(seedset, c, color_rs[c],
                                                       color_ids[c]);
                }

                return {influences, seedset};
            }
        }

        std::tie(color_rs, color_ids) = merge_rr_sets(color_rs, tmp_color_rs);
        seedset = build_seedset(budgets, color_rs, color_ids);
        cov = estimate_influence(seedset, color_rs, color_ids);
    }

    std::cout << "Returning seedset late" << std::endl;

    std::vector<double> influences(color_rs.size());
    for (color c{0}; c < color_rs.size(); c++) {
        influences[c] =
            estimate_influence(seedset, c, color_rs[c], color_ids[c]);
    }

    return {influences, seedset};
}
/**
 * Checks if the number of the RR sets surpasses the maximum number established
 * by the older paper.
 *
 * |Rc| > (8+2ε)Wc*(ln(2/δ)+ln(n over k))/ε².
 */
bool ris::stopping_criterion(std::vector<unsigned int> budgets, double epsilon,
                             double delta, const color_rr_sets &color_rs) {
    int number_of_users = num_vertices(g);
    int budgets_sum = std::accumulate(budgets.begin(), budgets.end(), 0);
    double constant =
        (8 + 2 * epsilon) *
        (std::log(2.0 / delta) + Math::logcnk(number_of_users, budgets_sum)) /
        (epsilon * epsilon);
    for (color c{0}; c < color_rs.size(); c++) {
        if (color_rs[c].size() <= constant) {
            return false;
        }
    }

    return true;
}

std::pair<color_rr_sets, color_user_rr_set_ids>
ris::merge_rr_sets(const color_rr_sets &color_rs,
                   const color_rr_sets &tmp_color_rs) {
    color_rr_sets new_color_rs(color_rs.size());
    color_user_rr_set_ids new_color_ids(color_rs.size());

    int number_of_users = num_vertices(g);
    for (color c{0}; c < color_rs.size(); c++) {
        new_color_rs[c].reserve(color_rs[c].size() + tmp_color_rs[c].size());

        std::copy(color_rs[c].begin(), color_rs[c].end(),
                  std::back_inserter(new_color_rs[c]));

        std::copy(tmp_color_rs[c].begin(), tmp_color_rs[c].end(),
                  std::back_inserter(new_color_rs[c]));

        new_color_ids[c].resize(num_vertices(g));

        for (int i{0}; i < new_color_rs[c].size(); i++) {
            for (user_distance &u_d : new_color_rs[c][i]) {
                new_color_ids[c][u_d.user].push_back(i);
            }
        }
    }

    return {new_color_rs, new_color_ids};
}

/**
 * Calculate the rrset coverages.
 */
double
ris::coverage(const std::unordered_map<vertex_descriptor, color> &seedset,
              const color_rr_sets &color_rs,
              const color_user_rr_set_ids &color_ids) {
    double influence{0};

    for (color c{0}; c < color_rs.size(); c++) {
        influence += coverage(seedset, c, color_rs[c], color_ids[c]);
    }

    return influence;
}

double
ris::coverage(const std::unordered_map<vertex_descriptor, color> &seedset,
              color c, const rr_sets &rs, const user_rr_set_ids &ids) {
    double rrset_sum{0};
    std::unordered_set<int> covered_rrsets;

    for (auto user : seedset) {
        if (user.second == c) {
            covered_rrsets.insert(ids[user.first].begin(),
                                  ids[user.first].end());
        }
    }

    for (auto rrset_id : covered_rrsets) {
        rrset_sum += estimate_influence(seedset, c, rs[rrset_id]);
    }
    return rrset_sum;
}

/**
 * Estimate the influence of and seedset according to its rrset coverage.
 *
 * The relation to the rrset coverage is only guaranteed to hold in the if all
 * nodes in the seedsets are of the same color.
 */
double ris::estimate_influence(
    const std::unordered_map<vertex_descriptor, color> &seedset,
    const color_rr_sets &color_rs, const color_user_rr_set_ids &color_ids) {
    double influence{0};

    for (color c{0}; c < color_rs.size(); c++) {
        influence += estimate_influence(seedset, c, color_rs[c], color_ids[c]);
    }

    return influence;
}

double ris::estimate_influence(
    const std::unordered_map<vertex_descriptor, color> &seedset, color c,
    const rr_sets &rs, const user_rr_set_ids &ids) {
    double rrset_sum = coverage(seedset, c, rs, ids);

    return color_sumimportances[c] * rrset_sum / rs.size();
}

/**
 * Returns the probability color c has to influence rrset_id, given the seedset
 * (which also contains the seeds of the other colors)
 *
 * Essentially this computes |S_c ∩ RR_i| / | S ∩ RR_i|.
 */
double ris::estimate_influence(
    const std::unordered_map<vertex_descriptor, color> &seedset, color c,
    const rr_set &r) {
    double c_size{0};
    double all_size{0};
    int closest_seed_distance{std::numeric_limits<int>::max()};

    for (const user_distance &u_d : r) {
        auto seed_it = seedset.find(u_d.user);
        if (seed_it != seedset.end()) {
            if (u_d.distance < closest_seed_distance) {
                closest_seed_distance = u_d.distance;
                c_size = 0;
                all_size = 0;
            }
            if (u_d.distance == closest_seed_distance) {
                if ((*seed_it).second == c) {
                    c_size++;
                }
                all_size++;
            }
        }
    }

    // If all_size == 0 also c_size is 0.
    return (all_size == 0) ? 0 : c_size / all_size;
}

/**
 * Extract the number of rr_sets each user is in for each color
 */
std::vector<std::vector<double>>
ris::get_color_user_rrsetdegrees(const color_rr_sets &color_rs,
                                 const color_user_rr_set_ids &color_ids) {
    std::vector<std::vector<double>> color_user_rrsetsdegrees(color_ids.size());
    for (color c{0}; c < color_ids.size(); c++) {
        color_user_rrsetsdegrees[c].resize(color_ids[c].size());
        for (int i{0}; i < color_ids[c].size(); i++) {
            color_user_rrsetsdegrees[c][i] = color_ids[c][i].size() *
                                             color_sumimportances[c] /
                                             color_rs[c].size();
        }
    }

    return color_user_rrsetsdegrees;
}

/**
 * Return the user with highest rrsetdegree over all colors.
 *
 * Given that the color still has a budget left.
 */
std::pair<vertex_descriptor, color> ris::get_best_user_color(
    const std::vector<unsigned int> &budgets,
    const std::vector<std::vector<double>> &color_user_rrsetsdegrees,
    const color_rr_sets &color_rs) {
    color c{-1};
    double current_max = std::numeric_limits<double>::lowest();
    std::vector<double>::const_iterator current_iterator;
    for (color i{0}; i < color_user_rrsetsdegrees.size(); i++) {
        auto t = std::max_element(color_user_rrsetsdegrees[i].begin(),
                                  color_user_rrsetsdegrees[i].end());

        // We have to normalize this with Wc and θ
        auto value = color_sumimportances[i] * (*t) / color_rs[i].size();

        if (budgets[i] > 0 && value > current_max) {
            current_max = value;
            current_iterator = t;
            c = i;
        }
    }
    assert(c != -1);

    return {current_iterator - color_user_rrsetsdegrees[c].begin(), c};
}

std::unordered_map<vertex_descriptor, int>
ris::build_seedset(const std::vector<unsigned int> &const_budgets,
                   const color_rr_sets &color_rs,
                   const color_user_rr_set_ids &color_ids) {
    std::unordered_map<vertex_descriptor, color> seedset;
    std::vector<unsigned int> budgets(const_budgets.begin(),
                                      const_budgets.end());

    auto color_user_rrsetsdegrees =
        get_color_user_rrsetdegrees(color_rs, color_ids);

    // We store the rrset_ids which we have processed already
    std::vector<std::vector<bool>> color_seen_rrsets(color_rs.size());
    for (color c{0}; c < color_rs.size(); c++) {
        color_seen_rrsets[c] = std::vector<bool>(color_rs[c].size());
    }

    while (std::accumulate(budgets.begin(), budgets.end(), 0) > 0) {
        vertex_descriptor seed;
        color c;
        std::tie(seed, c) =
            get_best_user_color(budgets, color_user_rrsetsdegrees, color_rs);

        seedset.insert({seed, c});
        budgets[c]--;

        // If we add a seed we have to adjust the rrsets for all colors
        for (color i{0}; i < color_user_rrsetsdegrees.size(); i++) {
            // A user can only be seed once, by setting its degree to 0 we
            // prevent it from being selected again
            color_user_rrsetsdegrees[i][seed] = 0;

            // For all rrsets, that the seed is part of
            for (int rrset_id : color_ids[i][seed]) {
                // If we have not processed the rrset before
                if (!color_seen_rrsets[i][rrset_id]) {
                    color_seen_rrsets[i][rrset_id] = true;

                    // Decrease the degree of all users, that are in those
                    // rrsets.
                    // Because the rrsets, that our seed is part of are not
                    // interesting anymore,
                    // we therefore "ignore" them by decreasing the degree.
                    for (const user_distance &u_d : color_rs[i][rrset_id]) {
                        color_user_rrsetsdegrees[i][u_d.user] -=
                            color_sumimportances[i] / color_rs[c].size();
                    }
                }
            }
        }
    }

    return seedset;
}

ris::ris(network g) : influence_maximization_algorithm(g) {
    int number_of_colors = get_number_of_colors(g);
    int number_of_users = num_vertices(g);

    for (color c{0}; c < number_of_colors; c++) {
        std::vector<double> tmp_importances;
        tmp_importances.reserve(number_of_users);

        BGL_FORALL_VERTICES(user, g, network) {
            tmp_importances.push_back(g[user].importances[c]);
        }

        color_distributions.push_back(std::discrete_distribution<>(
            tmp_importances.begin(), tmp_importances.end()));
        color_sumimportances.push_back(std::accumulate(
            tmp_importances.begin(), tmp_importances.end(), 0.0));
    }
}

/**
 * Generates n RR sets for each color.
 *
 * Returns the RR sets and a mapping from users to the RR sets they are part of.
 */
std::pair<std::vector<rr_sets>, std::vector<user_rr_set_ids>>
ris::generate_rrsets(const int n) {
    auto C = get_number_of_colors(g);

    std::vector<rr_sets> color_r(C);
    std::vector<user_rr_set_ids> color_ids(C);

    for (int c{0}; c < C; c++) {
        auto r_ids = generate_rrsets(n, c);
        color_r[c] = r_ids.first;
        color_ids[c] = r_ids.second;
    }

    return {color_r, color_ids};
}

/**
 * Generates n RR sets for color c.
 *
 * Returns the RR sets and a mapping from users to the RR sets they are part of.
 */
std::pair<rr_sets, user_rr_set_ids> ris::generate_rrsets(const int n,
                                                         const color c) {
    rr_sets r(n);
    user_rr_set_ids ids(num_vertices(g));

    for (int i{0}; i < n; i++) {
        r[i] = reverse_random_propagation(random_user(c));
        for (user_distance &u_d : r[i]) {
            ids[u_d.user].push_back(i);
        }
    }

    return {r, ids};
}

vertex_descriptor ris::random_user(const color c) {
    return color_distributions[c](generator);
}

double ris::Math::logcnk(int n, int k) {
    double ans = 0;
    for (int i = n - k + 1; i <= n; i++) {
        ans += std::log(i);
    }
    for (int i = 1; i <= k; i++) {
        ans -= std::log(i);
    }
    return ans;
}
}
}
