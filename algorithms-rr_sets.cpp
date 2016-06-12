#include "algorithms-rr_sets.hpp"

#include "Graph_reader.hpp"

#include <boost/graph/iteration_macros.hpp>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>

namespace gsinfmax {
namespace algorithms {
namespace {
using color = int;
}
std::unordered_map<vertex_descriptor, int>
rr_sets::maximize_influence(std::vector<unsigned int> budgets) {
    double epsilon = 0.1;

    auto opt_lower_bound = estimate_lower_bound(epsilon, budgets);
    theta_based_on_lower_bound(epsilon, opt_lower_bound, budgets);

    auto seedset = build_seedset(budgets);

    return seedset;
}

std::vector<double>
rr_sets::estimate_lower_bound(double epsilon,
                              std::vector<unsigned int> budgets) {
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
bool rr_sets::assure_theta_i(const std::vector<int> &color_i,
                             const std::vector<unsigned int> &budgets,
                             double epsilon_prime) {
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
                       Math::logcnk(color_sumimportances[c], budgets[c]) +
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
void rr_sets::theta_based_on_lower_bound(double epsilon,
                                         std::vector<double> color_OPT_prime,
                                         std::vector<unsigned int> budgets) {
    int number_of_colors = get_number_of_colors(g);
    double e = std::exp(1);

    for (color c{0}; c < number_of_colors; c++) {
        double alpha =
            std::sqrt(std::log(color_sumimportances[c]) + std::log(2));
        double beta = std::sqrt(
            (1 - 1 / e) * (Math::logcnk(color_sumimportances[c], budgets[c]) +
                           std::log(color_sumimportances[c]) + std::log(2)));

        int R = 2.0 * color_sumimportances[c] *
                std::sqrt((1 - 1 / e) * alpha + beta) / color_OPT_prime[c] /
                epsilon / epsilon;

        build_hypergraph_r_color(R, c);
    }
}

/**
 * Estimate the influence of and seedset according to its rrset coverage.
 *
 * The relation to the rrset coverage is only guaranteed to hold in the if all
 * nodes in the seedsets are of the same color.
 */
double rr_sets::estimate_influence(
    std::unordered_map<vertex_descriptor, color> seedset) {
    double influence{0};
    int number_of_colors = get_number_of_colors(g);

    for (color c{0}; c < number_of_colors; c++) {
        influence += estimate_influence(seedset, c);
    }

    return influence;
}

double rr_sets::estimate_influence(
    std::unordered_map<vertex_descriptor, color> seedset, color c) {
    double rrset_sum{0};
    std::unordered_set<int> covered_rrsets;

    for (auto user : seedset) {
        if (user.second == c) {
            covered_rrsets.insert(color_user_rrsets[c][user.second].begin(),
                                  color_user_rrsets[c][user.second].end());
        }
    }

    for (auto rrset_id : covered_rrsets) {
        rrset_sum += estimate_influence(seedset, c, rrset_id);
    }

    return color_sumimportances[c] * rrset_sum / color_hypergraphs[c].size();
}

/**
 * Returns the probability color c has to influence rrset_id, given the seedset
 * (which also contains the seeds of the other colors)
 *
 * Essentially this computes |S_c ∩ RR_i| / | S ∩ RR_i|.
 * TODO: We should include only the seeds with the same minimal distance to v
 * (where v is the user from which the rrset was sampled)
 */
double rr_sets::estimate_influence(
    std::unordered_map<vertex_descriptor, color> seedset, color c,
    int rrset_id) {
    double c_size{0};
    double all_size{0};

    for (auto influencer : color_hypergraphs[c][rrset_id]) {
        auto seed_it = seedset.find(influencer);
        if (seed_it != seedset.end()) {
            if ((*seed_it).second == c) {
                c_size++;
            }
            all_size++;
        }
    }

    // If all_size == 0 also c_size is 0.
    return (all_size == 0) ? 0 : c_size / all_size;
}

/**
 * Extract the number of rr_sets each user is in for each color
 */
std::vector<std::vector<int>> rr_sets::get_color_user_rrsetdegrees() {
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
std::pair<vertex_descriptor, color> rr_sets::get_best_user_color(
    const std::vector<unsigned int> budgets,
    const std::vector<std::vector<int>> color_user_rrsetsdegrees) {
    auto number_of_colors = get_number_of_colors(g);

    color c{-1};
    int current_max = std::numeric_limits<int>::lowest();
    std::vector<int>::const_iterator current_iterator;
    for (color i{0}; i < number_of_colors; i++) {
        auto t = std::max_element(color_user_rrsetsdegrees[i].begin(),
                                  color_user_rrsetsdegrees[i].end());

        if (budgets[i] > 0 && *t > current_max) {
            current_max = *t;
            current_iterator = t;
            c = i;
        }
    }
    assert(c != -1);

    return {current_iterator - color_user_rrsetsdegrees[c].begin(), c};
}

std::unordered_map<vertex_descriptor, int>
rr_sets::build_seedset(std::vector<unsigned int> budgets) {
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
                    for (int node : color_hypergraphs[i][rrset_id]) {
                        color_user_rrsetsdegrees[i][node]--;
                    }
                }
            }
        }
    }

    return seedset;
}

rr_sets::rr_sets(network g) : influence_maximization_algorithm(g) {
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
        color_hypergraphs.insert({c, hypergraph{}});

        // Initialize the mapping from users to the iterations in which they
        // were influencers
        color_user_rrsets.insert({c, hypergraph(number_of_users)});
    }
}
/**
 * Makes sure that the hypergraph for color c contains at least r edges.
 *
 * Returns true, if new edges were added.
 */
bool rr_sets::build_hypergraph_r_color(const int r, const color c) {
    int previous_size = color_hypergraphs[c].size();
    bool edges_added{false};

    for (int i{previous_size}; i < r; i++) {
        edges_added = true;
        add_hypergraph_edge(random_user(c), c);
    }

    int number_of_added_influencers{0};
    for (int i{previous_size}; i < r; i++) {
        for (vertex_descriptor influencer : color_hypergraphs[c][i]) {
            color_user_rrsets[c][influencer].push_back(i);
            number_of_added_influencers++;
        }
    }

    return edges_added;
}

void rr_sets::add_hypergraph_edge(const vertex_descriptor user, const color c) {
    std::unordered_map<vertex_descriptor, color> s;
    s.insert({user, c});

    auto colored_influencers = random_propagation(s);
    std::vector<vertex_descriptor> influencers;
    influencers.reserve(colored_influencers.size());

    std::transform(
        colored_influencers.begin(), colored_influencers.end(),
        std::back_inserter(influencers),
        [](const std::pair<vertex_descriptor, color> &e) { return e.first; });

    color_hypergraphs[c].push_back(influencers);
}

vertex_descriptor rr_sets::random_user(const color c) {
    return color_distributions[c](generator);
}

double rr_sets::Math::logcnk(int n, int k) {
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
