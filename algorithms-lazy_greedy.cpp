#include "algorithms-lazy_greedy.hpp"

#include "Graph_reader.hpp"
#include "algorithms.hpp"
#include "analysis.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <boost/graph/iteration_macros.hpp>
#include <deque>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace gsinfmax {
namespace algorithms {
namespace {
using color = int;
using importance = double;

int special_color{-1};
}

void lazy_greedy::disable_generate_statistics() { generate_statistics = false; }

void lazy_greedy::enable_generate_statistics() { generate_statistics = true; }

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

/**
 * Maximize the influence of the colors according to their budget with the
 * baseline lazy greedy algorithm.
 *
 * budgets should contain the budget of each color.
 *
 * The baseline algorithm optimizes each color individually.
 * This means, that users that are already selected as seeds of a different
 * color are ignored
 * and especially they expected propagation from those seeds of a different
 * color is not taken
 * into account.
 */
std::unordered_map<vertex_descriptor, color>
lazy_greedy::maximize_influence_baseline(std::vector<unsigned int> budgets) {
    // Sum of budgets should not exceed network size
    assert(std::accumulate(budgets.begin(), budgets.end(), 0) <=
           num_vertices(g));

    int number_of_colors = get_number_of_colors(g);
    std::unordered_map<vertex_descriptor, color> seedset;
    assert(budgets.size() == number_of_colors);

    // Prepare logging
    int total_number_of_mgs_updates{0};
    int number_of_mgs_updates_this_iteration{0};
    logfile = get_logfile("lazy_greedy_baseline");

    for (int c{0}; c < number_of_colors; c++) {
        std::unordered_map<vertex_descriptor, color> seedset_c;
        std::multimap<importance, std::pair<vertex_descriptor, color>> queue_c;
        int iteration{0};
        current_iteration =
            iteration; ///< Store iteration in class scope for logging

        update_statusline_seeds(seedset.size(), budgets,
                                number_of_mgs_updates_this_iteration,
                                total_number_of_mgs_updates);

        // Initially fill queue
        BGL_FORALL_VERTICES(user, g, network) {
            // Ignore users, that are already seeds of a different color
            if (seedset.find(user) == seedset.end()) {
                auto mgs = marginal_influence_gain(user, seedset_c, seedset);
                queue_c.insert({mgs[c], {user, iteration}});
            }
        }

        while (budgets[c] > 0) {
            vertex_descriptor user;
            int user_iteration;
            std::tie(user, user_iteration) = (*queue_c.begin()).second;
            queue_c.erase(queue_c.begin());

            if (user_iteration == iteration) {
                seedset_c.insert({user, c});
                budgets[c]--;
                iteration++;
                current_iteration =
                    iteration; ///< Store iteration in class scope for logging

                number_of_mgs_updates_this_iteration = 0;

                assert(budgets[c] >= 0);
            } else {
                auto mgs = marginal_influence_gain(user, seedset_c, seedset);
                queue_c.insert({mgs[c], {user, iteration}});

                total_number_of_mgs_updates++;
                number_of_mgs_updates_this_iteration++;
            }

            update_statusline_seeds(seedset.size() + seedset_c.size(), budgets,
                                    number_of_mgs_updates_this_iteration,
                                    total_number_of_mgs_updates);
        }

        // Log the influence estimation of the final local seedset
        if (generate_statistics) {
            influence(seedset_c, seedset);
        }

        seedset.insert(seedset_c.begin(), seedset_c.end());
    }

    // Statusline will not be updated anymore, therefore print newline
    std::cout << std::endl;

    // Log the influence estimation of the final seedset
    if (generate_statistics) {
        current_iteration =
            -1; ///< -1 denotes that the final seedset is estimated / logged
        influence(seedset);
    }

    return seedset;
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
std::unordered_map<vertex_descriptor, color>
lazy_greedy::maximize_influence(std::vector<unsigned int> budgets) {
    // Sum of budgets should not exceed network size
    assert(std::accumulate(budgets.begin(), budgets.end(), 0) <=
           num_vertices(g));

    int number_of_colors = get_number_of_colors(g);
    std::unordered_map<vertex_descriptor, color> seedset;
    int iteration{0};

    // One Queue per color
    auto queues = std::vector<
        std::multimap<importance, std::pair<vertex_descriptor, color>>>(
        number_of_colors);
    for (int c{0}; c < number_of_colors; c++) {
        queues[c] =
            std::multimap<importance,
                          std::pair<vertex_descriptor, color>>(); // Queue of
                                                                  // users to
                                                                  // evaluate
                                                                  // for color c
    }

    // Prepare logging
    int total_number_of_mgs_updates{0};
    int number_of_mgs_updates_this_iteration{0};
    current_iteration = iteration;
    logfile = get_logfile("lazy_greedy_adapted");

    // Initially fill Q
    BGL_FORALL_VERTICES(user, g, network) {
        auto mgs = marginal_influence_gain(user);
        for (int c{0}; c < mgs.size(); c++) {
            queues[c].insert({mgs[c], {user, iteration}});
        }
    }

    update_statusline_seeds(seedset.size(), budgets,
                            number_of_mgs_updates_this_iteration,
                            total_number_of_mgs_updates);

    while (std::accumulate(budgets.begin(), budgets.end(), 0) > 0) {
        // Find color c
        // - that still has a budget priority queue
        // - whose priority queue has the maximum first element
        color c{-1};
        importance current_max = std::numeric_limits<importance>::lowest();
        for (int i{0}; i < queues.size(); i++) {
            if (budgets[i] > 0 && (*queues[i].cbegin()).first > current_max) {
                current_max = (*queues[i].cbegin()).first;
                c = i;
            }
        }
        assert(c != -1);

        vertex_descriptor user;
        int user_iteration;
        std::tie(user, user_iteration) = (*queues[c].begin()).second;

        // Delete the user in the queue
        // In other queues the user could be at an arbitrary position
        for (int i{0}; i < queues.size(); i++) {
            for (auto it = queues[i].begin();; ++it) {
                if ((*it).second.first == user) {
                    queues[i].erase(it);
                    break;
                }
            }
        }

        if (iteration == user_iteration) {
            seedset.insert({user, c});
            budgets[c]--;
            iteration++;
            current_iteration =
                iteration; ///< Store iteration in class scope for logging

            number_of_mgs_updates_this_iteration = 0;

            assert(budgets[c] >= 0);
        } else {
            auto mgs = marginal_influence_gain(user, seedset);
            for (int i{0}; i < mgs.size(); i++) {
                // Insert the user with updated values
                queues[i].insert({mgs[i], {user, iteration}});
            }

            total_number_of_mgs_updates++;
            number_of_mgs_updates_this_iteration++;
        }

        update_statusline_seeds(seedset.size(), budgets,
                                number_of_mgs_updates_this_iteration,
                                total_number_of_mgs_updates);
    }

    // Statusline will not be updated anymore, therefore print newline
    std::cout << std::endl;

    // Log the influence estimation of the final seedset
    if (generate_statistics) {
        current_iteration =
            -1; ///< -1 denotes that the final seedset is estimated / logged
        influence(seedset);
    }

    return seedset;
};

/**
 * Perform a single simulation for the each user in the set s.
 *
 * The users in ignore are ignored (as if they wouldn't exist in the graph).
 */
std::unordered_map<vertex_descriptor, color> lazy_greedy::random_propagation(
    const std::unordered_map<vertex_descriptor, color> &s,
    const std::unordered_map<vertex_descriptor, color> &ignore) {
    std::unordered_map<vertex_descriptor, color> propagation_set(s);
    std::deque<std::pair<vertex_descriptor, color>> queue;
    for (auto &user_color_pair : s) {
        queue.push_back(user_color_pair);
    }

    // Shuffle queue
    std::shuffle(queue.begin(), queue.end(), generator);

    while (!queue.empty()) {
        auto user_color_pair = queue.front();

        auto neighbors = random_neighbors(user_color_pair.first);
        for (auto &neighbor : neighbors) {
            if (propagation_set.find(neighbor) == propagation_set.end() &&
                ignore.find(neighbor) == ignore.end()) {
                queue.push_back({neighbor, user_color_pair.second});
                propagation_set.insert({neighbor, user_color_pair.second});
            }
        }
        queue.pop_front();
    }
    return propagation_set;
}

/**
 * Return a random set of neighbors of user according to the edge propabilities.
 */
std::unordered_set<vertex_descriptor>
lazy_greedy::random_neighbors(const vertex_descriptor user) {
    std::unordered_set<vertex_descriptor> neighbors;

    BGL_FORALL_OUTEDGES(user, friendship, g, network) {
        assert(0 <= g[friendship].weight <= 1);
        std::bernoulli_distribution bernoulli(g[friendship].weight);

        // Does user influence his friend?
        if (bernoulli(generator)) {
            neighbors.insert(target(friendship, g));
        }
    }

    return neighbors;
}

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

        // Log all MC simulated influence estimations
        if (generate_statistics) {
            logfile << current_iteration << "\t" << iteration << "\t";
            for (int i{0}; i < importance.size(); i++) {
                logfile << importance[i] << "\t";
            }
            logfile << "\n";
        }

        influence += importance;
    }

    return influence / number_of_mc_sim;
}

/**
 * Calculate the marginal influence gain of user u compared to set s, ignoring
 * the set ignore.
 *
 * Set s contains colored users, whereas u is uncolored.
 * The returned vector contains the marginal influences,
 * that u would provide given that it is added with the respective color.
 */
Eigen::ArrayXd lazy_greedy::marginal_influence_gain(
    const vertex_descriptor u, std::unordered_map<vertex_descriptor, color> s,
    const std::unordered_map<vertex_descriptor, color> ignore) {
    Eigen::ArrayXd gain = Eigen::ArrayXd::Zero(get_number_of_colors(g));

    for (int iteration{0}; iteration < number_of_mc_sim; iteration++) {
        auto s_importance =
            importance_of_user_set(random_propagation(s, ignore));

        // Log all MC simulated influence estimations
        // Skip logging if s is empty
        if (generate_statistics && !s.empty()) {
            logfile << current_iteration << "\t" << iteration << "\t";
            for (int i{0}; i < s_importance.size(); i++) {
                logfile << s_importance[i] << "\t";
            }
            logfile << "\n";
        }

        // special_color is a special value for importance_of_user_set
        // The returned vector contains the marginal influences,
        // that u would provide given that it is added with the respective
        // color.
        s.insert({u, special_color});
        auto su_importance =
            importance_of_user_set(random_propagation(s, ignore));
        s.erase(u);

        gain += (su_importance - s_importance);
    }

    return gain / number_of_mc_sim;
};

/**
 * Computes the importance of a set of users.
 *
 * Users of special color are added to the importance of each of the normal
 * colors.
 */
Eigen::ArrayXd lazy_greedy::importance_of_user_set(const auto &set) {
    int number_of_colors = get_number_of_colors(g);
    Eigen::ArrayXd importances = Eigen::ArrayXd::Zero(number_of_colors);

    for (const auto &user_color : set) {
        vertex_descriptor user = user_color.first;
        color c = user_color.second;

        if (c != special_color) {
            importances[c] += g[user].importances[c];
        } else {
            importances += g[user].importances;
        }
    }

    return importances;
};
}
}
