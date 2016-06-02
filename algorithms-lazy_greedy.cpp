#include "algorithms-lazy_greedy.hpp"

#include "algorithms.hpp"
#include "Graph_reader.hpp"

#include <boost/graph/iteration_macros.hpp>
#include <Eigen/Core>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <random>
#include <algorithm>
#include <numeric>
#include <functional>

namespace gsinfmax { namespace algorithms {
	namespace {
		using color = int;
		using importance = double;

		int special_color {-1};
	}

	/**
	 * Maximize the influence of the colors according to their budget with the baseline lazy greedy algorithm.
	 *
	 * budgets should contain the budget of each color.
	 *
	 * The baseline algorithm optimizes each color individually.
	 * This means, that users that are already selected as seeds of a different color are ignored
	 * and especially they expected propagation from those seeds of a different color is not taken
	 * into account.
	 */
	std::unordered_map<vertex_descriptor, color> lazy_greedy::maximize_influence_baseline(std::vector<unsigned int> budgets) {
		// Sum of budgets should not exceed network size
		assert(std::accumulate(budgets.begin(), budgets.end(), 0) <= num_vertices(g));

		int number_of_colors = get_number_of_colors(g);
		std::unordered_map<vertex_descriptor, color> seedset;
		assert(budgets.size() == number_of_colors);

		for (int c {0}; c < number_of_colors; c++) {
			std::unordered_map<vertex_descriptor, color> seedset_c;
			std::multimap<importance, std::pair<vertex_descriptor, color>> queue_c;
			int iteration {0};

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

					assert(budgets[c] >= 0);
				} else {
					auto mgs = marginal_influence_gain(user, seedset_c, seedset);
					queue_c.insert({mgs[c], {user, iteration}});
				}
			}

			seedset.insert(seedset_c.begin(), seedset_c.end());
		}

		return seedset;
	}

	/**
	 * Maximize the influence of the colors according to their budget with the adapted lazy greedy algorithm.
	 *
	 * budgets should contain the budget of each color.
	 *
	 * The adapted algorithm optimizes the global target function (sum of all colors).
	 */
	std::unordered_map<vertex_descriptor, color> lazy_greedy::maximize_influence(std::vector<unsigned int> budgets) {
		// Sum of budgets should not exceed network size
		assert(std::accumulate(budgets.begin(), budgets.end(), 0) <= num_vertices(g));

		int number_of_colors = get_number_of_colors(g);
		std::unordered_map<vertex_descriptor, color> seedset;
		int iteration {0};

		// One Queue per color
		auto queues = std::vector<std::multimap<importance, std::pair<vertex_descriptor, color>>>(number_of_colors);
		for (int c {0}; c < number_of_colors; c++) {
			queues[c] = std::multimap<importance, std::pair<vertex_descriptor, color>>(); // Queue of users to evaluate for color c
		}

		// Initially fill Q
		BGL_FORALL_VERTICES(user, g, network) {
			auto mgs = marginal_influence_gain(user);
			for (int c {0}; c < mgs.size(); c++) {
				queues[c].insert({mgs[c], {user, iteration}});
			}
		}

		while (std::accumulate(budgets.begin(), budgets.end(), 0) > 0) {
			// Find color c
			// - that still has a budget priority queue
			// - whose priority queue has the maximum first element
			color c {-1};
			importance current_max = std::numeric_limits<importance>::lowest();
			for (int i {0}; i < queues.size(); i++) {
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
			for (int i {0}; i < queues.size(); i++) {
				for (auto it = queues[i].begin(); ; ++it) {
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

				assert(budgets[c] >= 0);
			} else {
				auto mgs = marginal_influence_gain(user, seedset);
				for (int i {0}; i < mgs.size(); i++) {
					// Insert the user with updated values
					queues[i].insert({mgs[i], {user, iteration}});
				}
			}
		}

		return seedset;
	};

	/**
	 * Perform a single simulation for the each user in the set s.
	 *
	 * The users in ignore are ignored (as if they wouldn't exist in the graph).
	 */
	std::unordered_map<vertex_descriptor, color> lazy_greedy::random_propagation(const std::unordered_map<vertex_descriptor, color>& s,
			const std::unordered_map<vertex_descriptor, color>& ignore) {
		std::unordered_map<vertex_descriptor, color> propagation_set(s);
		std::queue<std::pair<vertex_descriptor, color>> queue;
		for (auto& user_color_pair : s) {
			queue.push(user_color_pair);
		}

		while(! queue.empty()) {
			auto user_color_pair = queue.front();

			auto neighbors = random_neighbors(user_color_pair.first);
			for (auto& neighbor : neighbors) {
				if (propagation_set.find(neighbor) == propagation_set.end() &&
						ignore.find(neighbor) == ignore.end()) {
					queue.push({neighbor, user_color_pair.second});
					propagation_set.insert({neighbor, user_color_pair.second});
				}
			}
			queue.pop();
		}
		return propagation_set;
	}

	/**
	 * Return a random set of neighbors of user according to the edge propabilities.
	 */
	std::unordered_set<vertex_descriptor> lazy_greedy::random_neighbors(const vertex_descriptor user) {
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
	 * Calculate the marginal influence gain of user u compared to set s, ignoring the set ignore.
	 *
	 * Set s contains colored users, whereas u is uncolored.
	 * The returned vector contains the marginal influences,
	 * that u would provide given that it is added with the respective color.
	 */
	Eigen::ArrayXd lazy_greedy::marginal_influence_gain(const vertex_descriptor u,
			std::unordered_map<vertex_descriptor, color> s,
			const std::unordered_map<vertex_descriptor, color> ignore) {
		Eigen::ArrayXd gain = Eigen::ArrayXd::Zero(get_number_of_colors(g));

		for (int iteration {0}; iteration < number_of_mc_sim; iteration++) {
			auto s_importance = importance_of_user_set(random_propagation(s, ignore));

			// special_color is a special value for importance_of_user_set
			// The returned vector contains the marginal influences,
			// that u would provide given that it is added with the respective color.
			s.insert({u, special_color});
			auto su_importance = importance_of_user_set(random_propagation(s, ignore));
			s.erase(u);

			gain += (su_importance - s_importance);
		}

		return gain / number_of_mc_sim;
	};

	/**
	 * Computes the importance of a set of users.
         *
         * Users of special color are added to the importance of each of the normal colors.
	 */
	Eigen::ArrayXd lazy_greedy::importance_of_user_set(const auto& set) {
		int number_of_colors = get_number_of_colors(g);
		Eigen::ArrayXd importances = Eigen::ArrayXd::Zero(number_of_colors);
		
		for (const auto& user_color: set) {
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
}}
