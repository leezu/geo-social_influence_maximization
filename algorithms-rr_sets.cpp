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

namespace gsinfmax {
namespace algorithms {
namespace {
using vector_vector = std::vector<std::vector<vertex_descriptor>>;
using rr_sets = std::vector<std::vector<user_distance>>;
using color = int;
}
std::pair<std::vector<double>, std::unordered_map<vertex_descriptor, int>>
ris::dssa(std::vector<unsigned int> budgets, double epsilon, double delta) {
    double Lambda = 2.0 * (2 * (std::exp(1)) - 2) * (1 + epsilon) *
                    (1 + epsilon) *
                    std::log(2.0 / delta * 1 / (epsilon * epsilon));
    int number_of_colors = get_number_of_colors(g);

    for (color c{0}; c < number_of_colors; c++) {
        build_hypergraph_r_color(Lambda, c);
    }
    std::cout << Lambda << std::endl;

    auto seedset = build_seedset(budgets);
    auto cov = estimate_influence(seedset);

    while (dssa_double(budgets, epsilon, delta)) {
        std::cout << "DSSA doubled" << std::endl;
        auto tmp_cov = estimate_influence(seedset, true);
        std::cout << "influence " << cov << " influence(tmp) " << tmp_cov
                  << std::endl;
        double epsilon_1 = cov / tmp_cov - 1;
        if (epsilon_1 <= epsilon) {
            double epsilon_2 = (epsilon - epsilon_1) / (2 * (1 + epsilon_1));
            double epsilon_3 = (epsilon - epsilon_1) / (2 * (1 - epsilon_1));

            double delta_1 =
                std::exp((-1.0) * (coverage(seedset) * epsilon_3 * epsilon_3 /
                                   (2.0 * 2.0 * (std::exp(1) - 2.0) *
                                    (1 + epsilon_1) * (1 + epsilon_2))));
            double delta_2 = std::exp(
                (-1.0) * (coverage(seedset, true) * epsilon_2 * epsilon_2 /
                          (2.0 * 2.0 * (std::exp(1) - 2.0) * (1 + epsilon_2))));

            std::cout << "coverage " << coverage(seedset) << " coverage(tmp) "
                      << coverage(seedset, true) << std::endl;

            std::cout << "delta_1 + delta_2 <= delta: " << delta_1 << " + "
                      << delta_2 << " <= " << delta << std::endl;

            if (delta_1 + delta_2 <= delta) {
                std::cout << "Returning seedset early" << std::endl;

                std::vector<double> influences;
                for (color c{0}; c < number_of_colors; c++) {
                    influences.push_back(estimate_influence(seedset, c));
                }

                return {influences, seedset};
            }
        }

        merge_RR_sets();
        seedset = build_seedset(budgets);
        cov = estimate_influence(seedset);
    }

    std::cout << "Returning seedset late" << std::endl;

    std::vector<double> influences;
    for (color c{0}; c < number_of_colors; c++) {
        influences.push_back(estimate_influence(seedset, c));
    }

    return {influences, seedset};
}
/**
 * Doubles the amount of RR-sets for each color that does not yet exceed the
 * maximimum.
 * Returns true, if at least one RR-set was doubled.
 *
 * It does not double the amount, if |Rc| > (8+2ε)Wc*(ln(2/δ)+ln(n over k))/ε².
 */
bool ris::dssa_double(std::vector<unsigned int> budgets, double epsilon,
                      double delta) {
    int number_of_users = num_vertices(g);
    bool doubled{false};

    for (color c{0}; c < color_user_rrsets.size(); c++) {
        double constant = (8 + 2 * epsilon) *
                          (std::log(2.0 / delta) +
                           Math::logcnk(number_of_users, budgets[c])) /
                          (epsilon * epsilon);
        if (color_hypergraphs[c].size() <= constant) {
            doubled = true;
            build_hypergraph_r_color(color_hypergraphs[c].size(), c, true);
        }
    }

    return doubled;
}

void ris::merge_RR_sets() {
    std::cout << "Merging RR sets" << std::endl;
    int number_of_users = num_vertices(g);
    for (color c{0}; c < color_hypergraphs.size(); c++) {
        int previous_size = color_hypergraphs[c].size();

        std::copy(tmp_color_hypergraphs[c].begin(),
                  tmp_color_hypergraphs[c].end(),
                  std::back_inserter(color_hypergraphs[c]));

        for (int i{previous_size}; i < color_hypergraphs[c].size(); i++) {
            for (user_distance &u_d : color_hypergraphs[c][i]) {
                color_user_rrsets[c][u_d.user].push_back(i);
            }
        }

        tmp_color_hypergraphs[c].clear();
        tmp_color_user_rrsets.erase(c);
        tmp_color_user_rrsets.insert({c, vector_vector(number_of_users)});
    }
}

std::unordered_map<vertex_descriptor, int>
ris::maximize_influence(std::vector<unsigned int> budgets) {
    double epsilon = 0.1;

    auto opt_lower_bound = estimate_lower_bound(epsilon, budgets);
    theta_based_on_lower_bound(epsilon, opt_lower_bound, budgets);

    auto seedset = build_seedset(budgets);

    return seedset;
}

std::vector<double>
ris::estimate_lower_bound(double epsilon, std::vector<unsigned int> budgets) {
    auto number_of_colors = get_number_of_colors(g);
    double epsilon_prime = epsilon * std::sqrt(2);

    // Initialize i for each color to 1
    std::vector<int> color_i(number_of_colors, 1);
    std::vector<double> color_OPT_prime(number_of_colors, -1);

    // Assure that there exist theta_i RR-sets for each color based on i in
    // color_i
    // While new edges had to be added for one of the colors, we need to check
    // if any of the color_i needs to be increased
    while (assure_theta_i(color_i, budgets, epsilon_prime)) {
        auto seedset = build_seedset(budgets);
        // Test, for which colors theta_i shall be increased
        // Increase for all those (in the while loop check)
        for (color c{0}; c < number_of_colors; c++) {
            double ept =
                estimate_influence(seedset, c) / color_sumimportances[c];
            color_OPT_prime[c] =
                ept * color_sumimportances[c] / (1 + epsilon_prime);

            // TODO: Shouldn't "1 / std::pow(2.0, x)" be multiplied with
            // (1+epsilon_prime) ?
            if (!(ept > 1 / std::pow(2.0, color_i[c]))) {
                color_i[c]++;
            }
        }
    }

    return color_OPT_prime;
}

/**
 * Assure that there are theta_i RR-sets for each color, where theta_i is
 * calculated for each color based on color_i.
 *
 * Returns true, if new RR-sets were created.
 */
bool ris::assure_theta_i(const std::vector<int> &color_i,
                         const std::vector<unsigned int> &budgets,
                         double epsilon_prime) {
    int number_of_users = num_vertices(g);
    bool edges_added{false};

    for (color c{0}; c < color_i.size(); c++) {
        // i does not have to be bigger than log_2(W) - 1
        auto i = color_i[c];
        if (i >= std::log2(color_sumimportances[c])) {
            i = std::log2(color_sumimportances[c]);
        }
        // Calculate theta_i according to equation 9: λ'/(n/2^i)
        int theta_i = (2 + 2 / 3 * epsilon_prime) *
                      (std::log(color_sumimportances[c]) +
                       Math::logcnk(number_of_users, budgets[c]) +
                       std::log(std::log2(color_sumimportances[c]))) *
                      std::pow(2.0, i) / (epsilon_prime * epsilon_prime);

        if (build_hypergraph_r_color(theta_i, c)) {
            edges_added = true;
        };
    }

    return edges_added;
}

/**
 * Ensure that there are at least R=λ* / OPT_prime RR-sets (according to
 * Equation 6, TangShiXiao2015)
 */
void ris::theta_based_on_lower_bound(double epsilon,
                                     std::vector<double> color_OPT_prime,
                                     std::vector<unsigned int> budgets) {
    int number_of_colors = get_number_of_colors(g);
    int number_of_users = num_vertices(g);
    double e = std::exp(1);

    for (color c{0}; c < number_of_colors; c++) {
        double alpha =
            std::sqrt(std::log(color_sumimportances[c]) + std::log(2));
        double beta = std::sqrt(
            (1 - 1 / e) * (Math::logcnk(number_of_users, budgets[c]) +
                           std::log(color_sumimportances[c]) + std::log(2)));

        int R = 2.0 * color_sumimportances[c] *
                std::sqrt((1 - 1 / e) * alpha + beta) / color_OPT_prime[c] /
                epsilon / epsilon;

        build_hypergraph_r_color(R, c);
    }
}

/**
 * Calculate the rrset coverages.
 */
double ris::coverage(std::unordered_map<vertex_descriptor, color> seedset,
                     bool tmp) {
    double influence{0};
    int number_of_colors = get_number_of_colors(g);

    for (color c{0}; c < number_of_colors; c++) {
        influence += coverage(seedset, c, tmp);
    }

    return influence;
}

double ris::coverage(std::unordered_map<vertex_descriptor, color> seedset,
                     color c, bool tmp) {
    double rrset_sum{0};
    std::unordered_set<int> covered_rrsets;

    for (auto user : seedset) {
        if (user.second == c) {
            if (!tmp) {
                covered_rrsets.insert(color_user_rrsets[c][user.first].begin(),
                                      color_user_rrsets[c][user.first].end());
            } else {
                covered_rrsets.insert(
                    tmp_color_user_rrsets[c][user.first].begin(),
                    tmp_color_user_rrsets[c][user.first].end());
            }
        }
    }

    for (auto rrset_id : covered_rrsets) {
        rrset_sum += estimate_influence(seedset, c, rrset_id, tmp);
    }
    return rrset_sum;
}

/**
 * Estimate the influence of and seedset according to its rrset coverage.
 *
 * The relation to the rrset coverage is only guaranteed to hold in the if all
 * nodes in the seedsets are of the same color.
 */
double
ris::estimate_influence(std::unordered_map<vertex_descriptor, color> seedset,
                        bool tmp) {
    double influence{0};
    int number_of_colors = get_number_of_colors(g);

    for (color c{0}; c < number_of_colors; c++) {
        influence += estimate_influence(seedset, c, tmp);
    }

    return influence;
}

double
ris::estimate_influence(std::unordered_map<vertex_descriptor, color> seedset,
                        color c, bool tmp) {
    double rrset_sum{0};
    std::unordered_set<int> covered_rrsets;

    for (auto user : seedset) {
        if (user.second == c) {
            if (!tmp) {
                covered_rrsets.insert(color_user_rrsets[c][user.first].begin(),
                                      color_user_rrsets[c][user.first].end());
            } else {
                covered_rrsets.insert(
                    tmp_color_user_rrsets[c][user.first].begin(),
                    tmp_color_user_rrsets[c][user.first].end());
            }
        }
    }

    for (auto rrset_id : covered_rrsets) {
        rrset_sum += estimate_influence(seedset, c, rrset_id, tmp);
    }

    if (!tmp) {
        return color_sumimportances[c] * rrset_sum /
               color_hypergraphs[c].size();
    } else {
        return color_sumimportances[c] * rrset_sum /
               tmp_color_hypergraphs[c].size();
    }
}

/**
 * Returns the probability color c has to influence rrset_id, given the seedset
 * (which also contains the seeds of the other colors)
 *
 * Essentially this computes |S_c ∩ RR_i| / | S ∩ RR_i|.
 * TODO: We should include only the seeds with the same minimal distance to v
 * (where v is the user from which the rrset was sampled)
 */
double
ris::estimate_influence(std::unordered_map<vertex_descriptor, color> seedset,
                        color c, int rrset_id, bool tmp) {
    double c_size{0};
    double all_size{0};
    int closest_seed_distance{std::numeric_limits<int>::max()};

    if (!tmp) {
        for (user_distance &u_d : color_hypergraphs[c][rrset_id]) {
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
    } else {
        for (user_distance &u_d : tmp_color_hypergraphs[c][rrset_id]) {
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
std::vector<std::vector<int>> ris::get_color_user_rrsetdegrees() {
    int number_of_colors = get_number_of_colors(g);
    int number_of_users = num_vertices(g);
    std::vector<std::vector<int>> color_user_rrsetsdegrees(number_of_colors);

    for (color c{0}; c < number_of_colors; c++) {
        std::vector<int> user_rrsetsdegrees(number_of_users);

        for (int i{0}; i < number_of_users; i++) {
            user_rrsetsdegrees[i] = color_user_rrsets[c][i].size();
        }

        color_user_rrsetsdegrees[c] = user_rrsetsdegrees;
    }

    return color_user_rrsetsdegrees;
}

/**
 * Return the user with highest rrsetdegree over all colors.
 *
 * Given that the color still has a budget left.
 */
std::pair<vertex_descriptor, color> ris::get_best_user_color(
    const std::vector<unsigned int> budgets,
    const std::vector<std::vector<int>> color_user_rrsetsdegrees) {
    auto number_of_colors = get_number_of_colors(g);

    color c{-1};
    int current_max = std::numeric_limits<int>::lowest();
    std::vector<int>::const_iterator current_iterator;
    for (color i{0}; i < number_of_colors; i++) {
        auto t = std::max_element(color_user_rrsetsdegrees[i].begin(),
                                  color_user_rrsetsdegrees[i].end());

        // We have to normalize this with Wc and θ
        auto value =
            color_sumimportances[i] * (*t) / color_hypergraphs[i].size();

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
ris::build_seedset(std::vector<unsigned int> budgets) {
    int number_of_colors = get_number_of_colors(g);
    std::unordered_map<vertex_descriptor, color> seedset;

    auto color_user_rrsetsdegrees = get_color_user_rrsetdegrees();

    // We store the rrset_ids which we have processed already
    std::vector<std::vector<bool>> color_seen_rrsets(number_of_colors);
    for (color c{0}; c < number_of_colors; c++) {
        color_seen_rrsets[c] = std::vector<bool>(color_hypergraphs[c].size());
    }

    while (std::accumulate(budgets.begin(), budgets.end(), 0) > 0) {
        vertex_descriptor seed;
        color c;
        std::tie(seed, c) =
            get_best_user_color(budgets, color_user_rrsetsdegrees);

        seedset.insert({seed, c});
        budgets[c]--;

        // If we add a seed we have to adjust the rrsets for all colors
        for (color i{0}; i < number_of_colors; i++) {
            // A user can only be seed once, by setting its degree to 0 we
            // prevent it from being selected again
            color_user_rrsetsdegrees[i][seed] = 0;

            // For all rrsets, that the seed is part of
            for (int rrset_id : color_user_rrsets[i][seed]) {
                // If we have not processed the rrset before
                if (!color_seen_rrsets[i][rrset_id]) {
                    color_seen_rrsets[i][rrset_id] = true;

                    // Decrease the degree of all users, that are in those
                    // rrsets.
                    // Because the rrsets, that our seed is part of are not
                    // interesting anymore,
                    // we therefore "ignore" them by decreasing the degree.
                    for (user_distance &u_d : color_hypergraphs[i][rrset_id]) {
                        color_user_rrsetsdegrees[i][u_d.user]--;
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

        BGL_FORALL_VERTICES(user, g, network) {
            tmp_importances.push_back(g[user].importances[c]);
        }

        color_distributions.insert(
            {c, std::discrete_distribution<>(tmp_importances.begin(),
                                             tmp_importances.end())});
        color_sumimportances.push_back(std::accumulate(
            tmp_importances.begin(), tmp_importances.end(), 0.0));

        // Initialize the hypergraphs
        color_hypergraphs.insert({c, rr_sets{}});
        tmp_color_hypergraphs.insert({c, rr_sets{}});

        // Initialize the mapping from users to the iterations in which they
        // were influencers
        color_user_rrsets.insert({c, vector_vector(number_of_users)});
        tmp_color_user_rrsets.insert({c, vector_vector(number_of_users)});
    }
}
/**
 * Makes sure that the hypergraph for color c contains at least r edges.
 *
 * Returns true, if new edges were added.
 */
bool ris::build_hypergraph_r_color(const int r, const color c, bool tmp) {
    int previous_size = color_hypergraphs[c].size();
    if (tmp) {
        previous_size = tmp_color_hypergraphs[c].size();
    }
    bool edges_added{false};

    for (int i{previous_size}; i < r; i++) {
        edges_added = true;
        add_hypergraph_edge(random_user(c), c, tmp);
    }

    for (int i{previous_size}; i < r; i++) {
        if (!tmp) {
            for (user_distance &u_d : color_hypergraphs[c][i]) {
                color_user_rrsets[c][u_d.user].push_back(i);
            }
        } else {
            for (user_distance &u_d : tmp_color_hypergraphs[c][i]) {
                tmp_color_user_rrsets[c][u_d.user].push_back(i);
            }
        }
    }

    return edges_added;
}

void ris::add_hypergraph_edge(const vertex_descriptor user, const color c,
                              bool tmp) {
    auto influencers_distance = reverse_random_propagation(user);

    if (!tmp) {
        color_hypergraphs[c].push_back(influencers_distance);
    } else {
        tmp_color_hypergraphs[c].push_back(influencers_distance);
    }
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
