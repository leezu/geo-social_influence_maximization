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

	std::unordered_map<vertex_descriptor, color> lazy_greedy::maximize_influence_baseline(std::vector<int> budgets) {
		// Sum of budgets should not exceed network size
		assert(std::accumulate(budgets.begin(), budgets.end(), 0) <= num_vertices(g));

		int number_of_parties = get_number_of_parties(g);
		std::unordered_map<vertex_descriptor, color> seedset;
		assert(budgets.size() == number_of_parties);

		for (int c {0}; c < number_of_parties; c++) {
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



	std::unordered_map<vertex_descriptor, color> lazy_greedy::maximize_influence(std::vector<int> budgets) {
		// Sum of budgets should not exceed network size
		assert(std::accumulate(budgets.begin(), budgets.end(), 0) <= num_vertices(g));

		int number_of_parties = get_number_of_parties(g);
		std::unordered_map<vertex_descriptor, color> Seedset;
		int iteration {0};

		// One Queue per color
		auto Qs = std::vector<std::multimap<importance, std::pair<vertex_descriptor, int>>>(number_of_parties);
		for (int i {0}; i < number_of_parties; i++) {
			Qs[i] = std::multimap<importance, std::pair<vertex_descriptor, int>>(); // Queue of users to evaluate for color i
		}

		// Initially fill Q
		BGL_FORALL_VERTICES(user, g, network) {
			auto mgs = marginal_influence_gain(user);
			for (int i {0}; i < mgs.size(); i++) {
				Qs[i].insert({mgs[i], {user, iteration}});
			}
		}

		while (std::accumulate(budgets.begin(), budgets.end(), 0) > 0) {
			// Find color c
			// - that still has a budget priority queue
			// - whose priority queue has the maximum first element
			color c {-1};
			importance current_max = std::numeric_limits<importance>::lowest();
			for (int i {0}; i < Qs.size(); i++) {
				if (budgets[i] > 0 && (*Qs[i].cbegin()).first > current_max) {
					current_max = (*Qs[i].cbegin()).first;
					c = i;
				}
			}
			assert(c != -1);

			auto mg_user_iteration_it = Qs[c].begin();
			auto mg = (*mg_user_iteration_it).first;
			vertex_descriptor user = (*mg_user_iteration_it).second.first;
			int users_iteration = (*mg_user_iteration_it).second.second;

			// Delete the user in the queue
			// In other queues the user could be at an arbitrary position
			for (int i {0}; i < Qs.size(); i++) {
				for (auto it = Qs[i].begin(); ; ++it) {
					if ((*it).second.first == user) {
						Qs[i].erase(it);
						break;
					}
				}
			}

			if (iteration == users_iteration) {
				Seedset.insert({user, c});
				budgets[c]--;
				iteration++;

				assert(budgets[c] >= 0);
			} else {
				auto mgs = marginal_influence_gain(user, Seedset);
				for (int i {0}; i < mgs.size(); i++) {
					// Insert the user with updated values
					Qs[i].insert({mgs[i], {user, iteration}});
				}
			}
		}

		return Seedset;
	};

	/**
	 * Perform a single simulation for the each user in the set S.
	 *
	 * The users in ignored (as if they wouldn't exist in the graph).
	 */
	std::unordered_map<vertex_descriptor, color> lazy_greedy::random_propagation(const std::unordered_map<vertex_descriptor, color>& S,
			const std::unordered_map<vertex_descriptor, color>& ignore) {
		std::unordered_map<vertex_descriptor, color> propagation_set(S);
		std::queue<std::pair<vertex_descriptor, color>> Q;
		for (auto& s : S) {
			Q.push(s);
		}

		while(! Q.empty()) {
			auto user_color_pair = Q.front();

			auto neighbors = random_neighbors(user_color_pair.first);
			for (auto& neighbor : neighbors) {
				if (propagation_set.find(neighbor) == propagation_set.end() &&
						ignore.find(neighbor) == ignore.end()) {
					Q.push({neighbor, user_color_pair.second});
					propagation_set.insert({neighbor, user_color_pair.second});
				}
			}
			Q.pop();
		}
		return propagation_set;
	}

	/**
	 * Return a random set of neighbors of user according to the edge propabilities.
	 *
	 * The returned set does not include any users that are also present in ignore.
	 */
	std::unordered_set<vertex_descriptor> lazy_greedy::random_neighbors(auto user) {
		auto neighbors = std::unordered_set<vertex_descriptor>();

		BGL_FORALL_OUTEDGES(user, friendship, g, network) {
			assert(0 <= g[friendship].weight <= 1);
			std::bernoulli_distribution bernoulli(g[friendship].weight);

			// Is the friend already influenced (or ignored)? If no, does user influence his friend?
			if (bernoulli(generator)) {
				neighbors.insert(target(friendship, g));
			}
		}

		return neighbors;
	}

	/**
	 * Calculate the marginal influence gain of set u compared to set S.
	 *
	 * Set s contains colored users, whereas u is uncolored.
	 * The returned vector contains the marginal influences,
	 * that u would provide given that it is added with the respective color.
	 */
	Eigen::ArrayXd lazy_greedy::marginal_influence_gain(const vertex_descriptor u,
			std::unordered_map<vertex_descriptor, color> s,
			const std::unordered_map<vertex_descriptor, color> ignore) {
		Eigen::ArrayXd gain = Eigen::ArrayXd::Zero(get_number_of_parties(g));

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
	 */
	Eigen::ArrayXd lazy_greedy::importance_of_user_set(const auto& set) {
		int number_of_parties = get_number_of_parties(g);
		Eigen::ArrayXd importances = Eigen::ArrayXd::Zero(number_of_parties);
		
		if (set.empty()) {
			return importances;
		}

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
