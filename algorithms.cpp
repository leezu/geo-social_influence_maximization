#include "algorithms.hpp"
#include "Graph_reader.hpp"

#include <algorithm>
#include <boost/graph/iteration_macros.hpp>
#include <unordered_map>
#include <vector>

namespace gsinfmax {
namespace algorithms {

/**
 * Performs a reverse random propagation from user.
 *
 * Returns a map from users to the iteration in which they were added.
 */
std::vector<user_distance>
influence_maximization_algorithm::reverse_random_propagation(
    vertex_descriptor user) {
    std::vector<user_distance> rr_set;
    rr_set.push_back({user, 0});

    std::deque<user_distance> queue;
    queue.push_back({user, 0});

    while (!queue.empty()) {
        auto u_d = queue.front();

        auto neighbors = reverse_random_neighbors(u_d.user);
        for (auto &neighbor : neighbors) {
            auto same_user = [&neighbor](auto element) {
                return element.user == neighbor;
            };

            if (std::find_if(rr_set.begin(), rr_set.end(), same_user) ==
                rr_set.end()) {
                queue.push_back({neighbor, u_d.distance + 1});
                rr_set.push_back({neighbor, u_d.distance + 1});
            }
        }
        queue.pop_front();
    }

    return rr_set;
}
/**
 * Perform a single simulation for the each user in the set s.
 *
 * The users in ignore are ignored (as if they wouldn't exist in the graph).
 */
std::vector<user_color> influence_maximization_algorithm::random_propagation(
    const std::unordered_map<vertex_descriptor, color> &s,
    const std::unordered_map<vertex_descriptor, color> &ignore) {
    std::vector<user_color> propagation_set(s.begin(), s.end());
    std::deque<user_color> queue;
    for (auto &u_c : propagation_set) {
        queue.push_back(u_c);
    }

    // Shuffle queue
    std::shuffle(queue.begin(), queue.end(), this->generator);

    while (!queue.empty()) {
        auto user_color_pair = queue.front();

        auto neighbors = random_neighbors(user_color_pair.user);
        for (auto &neighbor : neighbors) {
            auto same_user = [&neighbor](auto element) {
                return element.user == neighbor;
            };

            if (std::find_if(propagation_set.begin(), propagation_set.end(),
                             same_user) == propagation_set.end() &&
                ignore.find(neighbor) == ignore.end()) {
                queue.push_back({neighbor, user_color_pair.color});
                propagation_set.push_back({neighbor, user_color_pair.color});
            }
        }
        queue.pop_front();
    }
    return propagation_set;
}

/**
 * Return a random set of neighbors of user according to the edge propabilities.
 */
std::vector<vertex_descriptor>
influence_maximization_algorithm::random_neighbors(
    const vertex_descriptor user) {
    std::vector<vertex_descriptor> neighbors;

    BGL_FORALL_OUTEDGES(user, friendship, g, network) {
        assert(0 <= this->g[friendship].weight <= 1);
        std::bernoulli_distribution bernoulli(this->g[friendship].weight);

        // Does user influence his friend?
        if (bernoulli(this->generator)) {
            neighbors.push_back(target(friendship, this->g));
        }
    }

    return neighbors;
}

/**
 * Return a random set of in-neighbors of user according to the edge
 * propabilities.
 */
std::vector<vertex_descriptor>
influence_maximization_algorithm::reverse_random_neighbors(
    const vertex_descriptor user) {
    std::vector<vertex_descriptor> neighbors;

    BGL_FORALL_INEDGES(user, friendship, g, network) {
        assert(0 <= this->g[friendship].weight <= 1);
        std::bernoulli_distribution bernoulli(this->g[friendship].weight);

        // Does user influence his friend?
        if (bernoulli(this->generator)) {
            neighbors.push_back(source(friendship, this->g));
        }
    }

    return neighbors;
}
}
}
