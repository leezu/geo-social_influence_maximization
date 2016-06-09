#include "algorithms-lazy_greedy.hpp"

#include "Graph_reader.hpp"
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

        std::multimap<importance, std::pair<vertex_descriptor, color>,
                      std::greater<importance>>
            queue_c;
        int iteration{0};

        update_statusline_seeds(seedset.size(), budgets,
                                number_of_mgs_updates_this_iteration,
                                total_number_of_mgs_updates);

        // Initially fill queue
        BGL_FORALL_VERTICES(user, g, network) {
            // Ignore users, that are already seeds of a different color
            if (seedset.find(user) == seedset.end()) {
                auto mgs = marginal_influence_gain(user, seedset_c,
                                                   sigma_seedset_c, seedset);
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

                budgets[c]--;
                iteration++;
                total_iteration++;

                number_of_mgs_updates_this_iteration = 0;

                assert(budgets[c] >= 0);
            } else {
                auto mgs = marginal_influence_gain(user, seedset_c,
                                                   sigma_seedset_c, seedset);
                queue_c.insert({mgs[c], {user, iteration}});

                total_number_of_mgs_updates++;
                number_of_mgs_updates_this_iteration++;
            }

            update_statusline_seeds(seedset.size() + seedset_c.size(), budgets,
                                    number_of_mgs_updates_this_iteration,
                                    total_number_of_mgs_updates);
        }

        seedset.insert(seedset_c.begin(), seedset_c.end());
    }

    // Statusline will not be updated anymore, therefore print newline
    std::cout << std::endl;

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
    Eigen::ArrayXd sigma_seedset = Eigen::ArrayXd::Zero(number_of_colors);

    int iteration{0};
    if (budgets.size() != number_of_colors) {
        throw std::runtime_error(
            "The number of budgets does not match the number of colors");
    }

    // One Queue per color
    auto queues = std::vector<
        std::multimap<importance, std::pair<vertex_descriptor, color>,
                      std::greater<importance>>>(number_of_colors);
    for (int c{0}; c < number_of_colors; c++) {
        queues[c] =
            std::multimap<importance, std::pair<vertex_descriptor, color>,
                          std::greater<importance>>(); // Queue of users to
                                                       // evaluate for color c
    }

    // Prepare logging
    int total_number_of_mgs_updates{0};
    int number_of_mgs_updates_this_iteration{0};

    // Initially fill Q
    auto number_of_users = num_vertices(g);
    int number_of_initialized_users{0};
    BGL_FORALL_VERTICES(user, g, network) {
        std::cout << "\rInitialized " << number_of_initialized_users << " of "
                  << number_of_users << " users." << std::flush;

        auto mgs = marginal_influence_gain(user, seedset, sigma_seedset);
        for (int c{0}; c < mgs.size(); c++) {
            queues[c].insert({mgs[c], {user, iteration}});
        }

        number_of_initialized_users++;
    }
    std::cout << std::endl;

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
            sigma_seedset = influence(seedset);

            if (generate_statistics) {
                log.sigma_set(sigma_seedset, iteration, "adapted_greedy_");
            }

            budgets[c]--;
            iteration++;

            number_of_mgs_updates_this_iteration = 0;

            assert(budgets[c] >= 0);
        } else {
            auto mgs = marginal_influence_gain(user, seedset, sigma_seedset);
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
Eigen::ArrayXd lazy_greedy::global_importance_of_user_set(const auto &set) {
    std::unordered_map<vertex_descriptor, color> special_users;
    std::unordered_map<vertex_descriptor, color> colored_users;

    for (const auto &user_color : set) {
        if (user_color.second == special_color) {
            special_users.insert(user_color);
        } else {
            colored_users.insert(user_color);
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
