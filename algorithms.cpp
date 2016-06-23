#include "algorithms.hpp"
#include "Graph_reader.hpp"

#include <boost/graph/iteration_macros.hpp>
#include <unordered_map>

namespace gsinfmax {
namespace algorithms {

/**
 * Performs a reverse random propagation from user.
 *
 * Returns a map from users to the iteration in which they were added.
 */
std::unordered_map<vertex_descriptor, int>
influence_maximization_algorithm::reverse_random_propagation(
    vertex_descriptor user) {
    std::unordered_map<vertex_descriptor, int> propagation_set;
    propagation_set.insert({user, 0});

    std::deque<std::pair<vertex_descriptor, int>> queue;
    queue.push_back({user, 0});

    while (!queue.empty()) {
        auto user_iteration_pair = queue.front();

        auto neighbors = reverse_random_neighbors(user_iteration_pair.first);
        for (auto &neighbor : neighbors) {
            if (propagation_set.find(neighbor) == propagation_set.end()) {
                queue.push_back({neighbor, user_iteration_pair.second + 1});
                propagation_set.insert(
                    {neighbor, user_iteration_pair.second + 1});
            }
        }
        queue.pop_front();
    }

    return propagation_set;
}
/**
 * Perform a single simulation for the each user in the set s.
 *
 * The users in ignore are ignored (as if they wouldn't exist in the graph).
 */
std::unordered_map<vertex_descriptor, influence_maximization_algorithm::color>
influence_maximization_algorithm::random_propagation(
    const std::unordered_map<vertex_descriptor, color> &s,
    const std::unordered_map<vertex_descriptor, color> &ignore) {
    std::unordered_map<vertex_descriptor, color> propagation_set(s);
    std::deque<std::pair<vertex_descriptor, color>> queue;
    for (auto &user_color_pair : s) {
        queue.push_back(user_color_pair);
    }

    // Shuffle queue
    std::shuffle(queue.begin(), queue.end(), this->generator);

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
