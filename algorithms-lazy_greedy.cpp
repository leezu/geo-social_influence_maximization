#include "algorithms-lazy_greedy.hpp"

#include "Graph_reader.hpp"
#include "analysis.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <boost/graph/iteration_macros.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_map>

namespace gsinfmax {
namespace algorithms {
namespace {
using color = int;
using importance = double;

int special_color{-1};
}

struct node {
    vertex_descriptor user;
    color c;
    double mgs;
    int iteration;

    node(vertex_descriptor v, color c, double mgs, int i)
        : user(v), c(c), mgs(mgs), iteration(i) {}
};

struct compare_node {
    bool operator()(const node &n1, const node &n2) const {
        return n1.mgs < n2.mgs;
    }
};

void lazy_greedy::disable_statusline() { statusline_enabled = false; }
void lazy_greedy::enable_statusline() { statusline_enabled = true; }
void lazy_greedy::disable_generate_statistics() { generate_statistics = false; }
void lazy_greedy::enable_generate_statistics() { generate_statistics = true; }
void lazy_greedy::set_number_of_mc(int mc) { number_of_mc_sim = mc; }

void lazy_greedy::update_statusline_seeds(
    const int current_seedset_size, const std::vector<unsigned int> &budgets,
    const int number_of_mgs_updates_current_iteration,
    const int total_number_of_mgs_updates) {
    int total_number_of_seeds =
        current_seedset_size +
        std::accumulate(budgets.begin(), budgets.end(), 0);

    if (statusline_printed) {
        std::cout << "\r";
    } else {
        statusline_printed = true;
    }

    std::cout << "Found " << current_seedset_size << " of "
              << total_number_of_seeds << " Seeds: Updated "
              << number_of_mgs_updates_current_iteration
              << " marginal gains this iteration (total: "
              << total_number_of_mgs_updates << ")" << std::flush;
}

void lazy_greedy::logger::sigma_set(const Eigen::ArrayXd &sigma_set,
                                    const int &iteration,
                                    const std::string &fname) {
    if (files.find(fname) == files.end()) {
        files.emplace(fname, get_logfile(fname));
    }

    auto f = files.find(fname);

    if (initialized.find(fname) == initialized.end()) {
        // Write header
        (*f).second << "Iteration\tColor\tInfluence\n";
        initialized.insert(fname);
    }

    for (color c{0}; c < sigma_set.size(); c++) {
        (*f).second << iteration << "\t" << c << "\t" << sigma_set[c] << "\n";
    }
}

/**
 * Maximize the influence of the colors according to their budget with the
 * baseline lazy greedy algorithm.
 *
 * budgets should contain the budget of each color.
 *
 * The baseline algorithm optimizes each color individually.
 * The seeds of colors that were optimized previously are known to the color
 * that is optimized in the given iteration.
 */
std::pair<std::vector<double>, std::unordered_map<vertex_descriptor, color>>
lazy_greedy::maximize_influence_baseline(std::vector<unsigned int> budgets) {
    // Sum of budgets should not exceed network size
    assert(std::accumulate(budgets.begin(), budgets.end(), 0) <=
           num_vertices(g));

    int number_of_colors = get_number_of_colors(g);
    std::unordered_map<vertex_descriptor, color> seedset;
    if (budgets.size() != number_of_colors) {
        throw std::runtime_error(
            "The number of budgets does not match the number of colors");
    }

    // Prepare logging
    int total_number_of_mgs_updates{0};
    int number_of_mgs_updates_this_iteration{0};

    int total_iteration{0};
    for (int c{0}; c < number_of_colors; c++) {
        std::unordered_map<vertex_descriptor, color> seedset_c;
        Eigen::ArrayXd sigma_seedset_c = Eigen::ArrayXd::Zero(number_of_colors);

        boost::heap::fibonacci_heap<node, boost::heap::compare<compare_node>>
            heap;
        int iteration{0};

        if (statusline_enabled) {
            update_statusline_seeds(seedset.size(), budgets,
                                    number_of_mgs_updates_this_iteration,
                                    total_number_of_mgs_updates);
        }

        // Initially fill queue
        BGL_FORALL_VERTICES(user, g, network) {
            // Ignore users, that are already seeds of a different color
            if (seedset.find(user) == seedset.end()) {
                auto mgs = marginal_influence_gain(user, seedset_c,
                                                   sigma_seedset_c, seedset);
                heap.push({user, c, mgs[c], iteration});
            }
        }

        while (budgets[c] > 0) {
            node top_node = heap.top();
            heap.pop();

            if (top_node.iteration == iteration) {
                seedset_c.insert({top_node.user, top_node.c});
                sigma_seedset_c = influence(seedset_c, seedset);

                if (generate_statistics) {
                    log.sigma_set(sigma_seedset_c, iteration,
                                  "baseline_greedy_" + std::to_string(c));

                    // Log influence spread including the seeds of other colors
                    std::unordered_map<vertex_descriptor, color> tmp_seedset;
                    tmp_seedset.insert(seedset_c.begin(), seedset_c.end());
                    tmp_seedset.insert(seedset.begin(), seedset.end());

                    auto sigma_tmp_seedset = influence(tmp_seedset);
                    log.sigma_set(sigma_tmp_seedset, total_iteration,
                                  "baseline_greedy");
                }

                budgets[top_node.c]--;
                iteration++;
                total_iteration++;

                number_of_mgs_updates_this_iteration = 0;

                assert(budgets[top_node.c] >= 0);
            } else {
                auto mgs = marginal_influence_gain(top_node.user, seedset_c,
                                                   sigma_seedset_c, seedset);
                heap.push({top_node.user, c, mgs[c], iteration});

                total_number_of_mgs_updates++;
                number_of_mgs_updates_this_iteration++;
            }

            if (statusline_enabled) {
                update_statusline_seeds(seedset.size() + seedset_c.size(),
                                        budgets,
                                        number_of_mgs_updates_this_iteration,
                                        total_number_of_mgs_updates);
            }
        }

        seedset.insert(seedset_c.begin(), seedset_c.end());
    }

    // Statusline will not be updated anymore, therefore print newline
    std::cout << std::endl;

    auto influences = influence(seedset);
    std::vector<double> influences_vector;
    for (color c{0}; c < influences.size(); c++) {
        influences_vector.push_back(influences[c]);
    }

    return {influences_vector, seedset};
}

/**
 * Maximize the influence of the colors according to their budget with the
 * adapted lazy greedy algorithm.
 *
 * budgets should contain the budget of each color.
 *
 * The adapted algorithm optimizes the global target function (sum of all
 * colors).
 */
std::pair<std::vector<double>, std::unordered_map<vertex_descriptor, color>>
lazy_greedy::maximize_influence(std::vector<unsigned int> budgets) {
    // Sum of budgets should not exceed network size
    assert(std::accumulate(budgets.begin(), budgets.end(), 0) <=
           num_vertices(g));

    int number_of_colors = get_number_of_colors(g);
    std::unordered_map<vertex_descriptor, color> seedset;
    Eigen::ArrayXd sigma_seedset = Eigen::ArrayXd::Zero(number_of_colors);

    int iteration{0};
    if (budgets.size() != number_of_colors) {
        throw std::runtime_error(
            "The number of budgets does not match the number of colors");
    }

    boost::heap::fibonacci_heap<node, boost::heap::compare<compare_node>> heap;

    // Prepare logging
    int total_number_of_mgs_updates{0};
    int number_of_mgs_updates_this_iteration{0};

    // Initially fill heap
    auto number_of_users = num_vertices(g);
    int number_of_initialized_users{0};
    BGL_FORALL_VERTICES(user, g, network) {
        std::cout << "\rInitialized " << number_of_initialized_users << " of "
                  << number_of_users << " users." << std::flush;

        auto mgs = marginal_influence_gain(user, seedset, sigma_seedset);
        Eigen::ArrayXd::Index max_index;
        while (true) {
            mgs.maxCoeff(&max_index);

            if (budgets[max_index] > 0) {
                break;
            } else if (mgs[max_index] ==
                       std::numeric_limits<double>::lowest()) {
                // All budgets are 0
                // return std::make_pair(0, seedset);
                return {std::vector<double>(number_of_colors), seedset};
            } else {
                mgs[max_index] = std::numeric_limits<double>::lowest();
            }
        }
        heap.push({user, (color)max_index, mgs[max_index], iteration});
        number_of_initialized_users++;
    }
    std::cout << std::endl;

    if (statusline_enabled) {
        update_statusline_seeds(seedset.size(), budgets,
                                number_of_mgs_updates_this_iteration,
                                total_number_of_mgs_updates);
    }

    while (std::accumulate(budgets.begin(), budgets.end(), 0) > 0) {
        node top_node = heap.top();
        heap.pop();

        if (iteration == top_node.iteration) {
            seedset.insert({top_node.user, top_node.c});
            sigma_seedset = influence(seedset);

            if (generate_statistics) {
                log.sigma_set(sigma_seedset, iteration, "adapted_greedy_");
            }

            budgets[top_node.c]--;
            iteration++;

            number_of_mgs_updates_this_iteration = 0;

            assert(budgets[top_node.c] >= 0);
        } else {
            auto mgs =
                marginal_influence_gain(top_node.user, seedset, sigma_seedset);

            Eigen::ArrayXd::Index max_index;
            while (true) {
                mgs.maxCoeff(&max_index);

                if (budgets[max_index] > 0) {
                    break;
                } else {
                    assert(mgs[max_index] !=
                           std::numeric_limits<double>::lowest());
                    mgs[max_index] = std::numeric_limits<double>::lowest();
                }
            }

            heap.push(
                {top_node.user, (color)max_index, mgs[max_index], iteration});

            total_number_of_mgs_updates++;
            number_of_mgs_updates_this_iteration++;
        }

        if (statusline_enabled) {
            update_statusline_seeds(seedset.size(), budgets,
                                    number_of_mgs_updates_this_iteration,
                                    total_number_of_mgs_updates);
        }
    }

    // Statusline will not be updated anymore, therefore print newline
    std::cout << std::endl;

    auto influences = influence(seedset);
    std::vector<double> influences_vector;
    for (color c{0}; c < influences.size(); c++) {
        influences_vector.push_back(influences[c]);
    }

    return {influences_vector, seedset};
};

/**
 * Estimate the influence of seedset with monte carlo simulations.
 */
Eigen::ArrayXd lazy_greedy::influence(
    const std::unordered_map<vertex_descriptor, color> &seedset,
    const std::unordered_map<vertex_descriptor, color> &ignore) {
    Eigen::ArrayXd influence = Eigen::ArrayXd::Zero(get_number_of_colors(g));

    for (int iteration{0}; iteration < number_of_mc_sim; iteration++) {
        auto importance =
            importance_of_user_set(random_propagation(seedset, ignore));

        influence += importance;
    }

    return influence / number_of_mc_sim;
}

/**
 * Calculate the marginal influence gain toward the global target function of
 * user u compared to set s, ignoring the set ignore.
 *
 * Set s contains colored users, whereas u is uncolored.
 * The returned vector contains the marginal influence gains toward the global
 * target function (sum of all colored influences), that u would provide given
 * that it is added with the respective color.
 */
Eigen::ArrayXd lazy_greedy::marginal_influence_gain(
    const vertex_descriptor u, std::unordered_map<vertex_descriptor, color> s,
    Eigen::ArrayXd sigma_s,
    const std::unordered_map<vertex_descriptor, color> ignore) {

    auto number_of_colors = get_number_of_colors(g);
    Eigen::ArrayXd gain = Eigen::ArrayXd::Zero(number_of_colors);

    for (int iteration{0}; iteration < number_of_mc_sim; iteration++) {
        // The returned vector contains the influence spread,
        // that u would provide given that it is added with the respective
        // color.
        // To distinguish u from already colored seed users we set it to
        // special_color.
        s.insert({u, special_color});
        auto sigma_su =
            global_importance_of_user_set(random_propagation(s, ignore));
        s.erase(u);

        gain += (sigma_su - sigma_s.sum());
    }

    return gain / number_of_mc_sim;
};

/**
 * Computes the global target function.
 *
 * Returns a vector where each element is the global target function given that
 * users of special color are taken as the respective color.
 *
 * The global target function is the sum of the values of the set for each
 * color.
 */
Eigen::ArrayXd
lazy_greedy::global_importance_of_user_set(const std::vector<user_color> &set) {
    std::vector<user_color> special_users;
    std::vector<user_color> colored_users;

    for (const auto &u_c : set) {
        if (u_c.color == special_color) {
            special_users.push_back(u_c);
        } else {
            colored_users.push_back(u_c);
        }
    }

    auto color_importances = importance_of_user_set(colored_users);
    auto special_importances = importance_of_user_set(special_users);

    return special_importances + color_importances.sum();
};

/**
 * Computes the value of a set of users for each color.
 *
 * Returns a vector where each element contains the value of the set for
 * this
 * color.
 *
 * Users of special color are added to the value of each of the normal
 * colors.
 */
Eigen::ArrayXd
lazy_greedy::importance_of_user_set(const std::vector<user_color> &set) {
    int number_of_colors = get_number_of_colors(g);
    Eigen::ArrayXd importances = Eigen::ArrayXd::Zero(number_of_colors);

    for (const auto &u_c : set) {
        if (u_c.color != special_color) {
            importances[u_c.color] += g[u_c.user].importances[u_c.color];
        } else {
            importances += g[u_c.user].importances;
        }
    }

    return importances;
};
}
}
