#include <boost/graph/adjacency_list.hpp>
#include <fstream>
#include <iostream>
#include <string>

#include "docopt.h"

#include "Graph_reader.hpp"
#include "algorithms-lazy_greedy.hpp"
#include "analysis.hpp"

static const char USAGE[] =
    R"(Geo-social influence maximization

    Usage:
      gsinfmax gowalla <edges> <locations> (--random-weights | --edge-weights=<w>) [--random-weights | --edge-weights=<w>] [--statistics]
      gsinfmax gowalla_austin_dallas <edges> <locations> <events> (--baseline | --adapted) [--random-weights | --edge-weights=<w>] [--statistics]
      gsinfmax (-h | --help)
      gsinfmax --version

    Options:
      -h --help                 Show this screen.
      --version                 Show version.
      --edge-weights=<w>        Set edge weights to W (if not specified in dataset)
                                [default: 0.1]
      --random-weights          Don't use egde-weights, but generate edge weights uniformly randomly.
      --statistics              Print out statistics for the graph being processed.
)";

using namespace gsinfmax;

void apply_reader_settings(auto &args, auto &reader) {
    if (args.at("--random-weights").asLong()) {
        reader.set_random_weights();
    }
    reader.set_weights(std::stod(args.at("--edge-weights").asStringList()[0]));
}

network get_network(auto args) {
    if (args.at("gowalla").asBool()) {
        reader::gowalla reader;
        apply_reader_settings(args, reader);
        return reader.read_network(args["<edges>"].asString(),
                                   args["<locations>"].asString());
    } else if (args.at("gowalla_austin_dallas").asBool()) {
        reader::gowalla_austin_dallas reader;
        apply_reader_settings(args, reader);

        // reader.set_no_drop_users_without_friends();

        return reader.read_network(args["<edges>"].asString(),
                                   args["<locations>"].asString(),
                                   args["<events>"].asString());
    } else {
        throw std::runtime_error("Couldn't get network for given argument");
    }
}

void print_seedset(auto seedset) {
    std::cout << "The seedset contains: ";
    for (auto &s : seedset) {
        std::cout << s.first << "<" << s.second << "> ";
    }
    std::cout << "\n";
}

int main(int argc, char *argv[]) {
    std::map<std::string, docopt::value> args = docopt::docopt(
        USAGE, {argv + 1, argv + argc},
        true,                                    // show help if requested
        "geo-social influence maximzation 0.1"); // version string

    network g = get_network(args);

    std::vector<unsigned int> budgets{1, 2, 3};
    if (args.at("--statistics").asBool()) {
        write_node_degrees(g);
        std::cout << "Average node degree: " << get_average_node_degree(g)
                  << std::endl;

        std::cout << "Graph has " << boost::num_vertices(g) << " vertices"
                  << " and " << boost::num_edges(g) << " edges" << std::endl;
    }

    auto lazy_greedy = algorithms::lazy_greedy(g);
    lazy_greedy.enable_generate_statistics();

    if (args.at("--adapted").asBool()) {
        std::cout << "Running adapted lazy greedy\n";
        auto seedset = lazy_greedy.maximize_influence(budgets);
        print_seedset(seedset);
    } else if (args.at("--baseline").asBool()) {
        std::cout << "Running baseline lazy greedy\n";
        auto seedset = lazy_greedy.maximize_influence_baseline(budgets);
        print_seedset(seedset);
    }

    return 0;
}
