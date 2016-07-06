#include <boost/chrono.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include "docopt.h"

#include "Graph_reader.hpp"
#include "algorithms-lazy_greedy.hpp"
#include "algorithms-rr_sets.hpp"
#include "algorithms.hpp"
#include "analysis.hpp"

static const char USAGE[] =
    R"(Geo-social influence maximization

    Usage:
      gsinfmax gowalla [options] <edges> <locations>
      gsinfmax gowalla_austin_dallas [options] <edges> <locations> <events>
      gsinfmax synthetic [options] <edges> <locations> <events>
      gsinfmax (-h | --help)
      gsinfmax --version

    Options:
      -h, --help                Show this screen.
      --version                 Show version.
      -A, --algorithm=<a>       Use algorithm a. Must be one of ris, classic
                                [default: ris]
      -V, --algo-version=<v>    Use the specified version of the algorithm. I.e. exact or heuristic version of ris.
                                Adapted or naive version of classic algorithm.
                                [default: 0]
      --edge-weights=<w>        Set edge weights to w (if not specified in dataset)
                                [default: 0.1]
      --random-weights          Don't use egde-weights, but generate edge weights uniformly randomly.
      --budget=<k>              Set the budget for each color to k
                                [default: 3]
      --epsilon=<e>             Set epsilon for the algorithms that support it.
                                [default: 0.1]
      --delta=<d>               Set delta for the algorithms that support it.
                                [default: 0.1]
      --no-statusline
      --statistics              Print out statistics for the graph being processed.
      --export-network          Save network as graphviz file to ./graph. Do not run algorithm.
      -N, --dataset-name=<n>    Use this name to identify the dataset used in logfiles.
                                [default: default]
      -C, --num-colors=<n>      Limit the number of colors to <n>. 0 is unlimited.
                                [default: 0]
      -M, --num-mc=<n>         Number of monte carlo iterations. Only applicable to lazy greedy algorithm.
                                [default: 100]
)";
using namespace gsinfmax;

void apply_reader_settings(auto &args, auto &reader) {
    if (args.at("--random-weights").asBool()) {
        reader.set_random_weights();
    }
    reader.limit_colors(args.at("--num-colors").asLong());
    reader.set_weights(std::stod(args.at("--edge-weights").asString()));
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
    } else if (args.at("synthetic").asBool()) {
        reader::synthetic reader;
        apply_reader_settings(args, reader);

        return reader.read_network(args["<edges>"].asString(),
                                   args["<locations>"].asString(),
                                   args["<events>"].asString());
    } else {
        throw std::runtime_error("Couldn't get network for given argument");
    }
}

void print_seedset(auto seedset) {
    std::cout << "The seedset contains: ";
    for (auto &s : seedset.second) {
        std::cout << s.first << "<" << s.second << "> ";
    }

    std::cout << "\nImportances: ";
    double sum{0};
    for (auto &s : seedset.first) {
        std::cout << s << " ";
        sum += s;
    }
    std::cout << "(sum: " << sum << ")" << std::endl;
}

int main(int argc, char *argv[]) {
    std::map<std::string, docopt::value> args = docopt::docopt(
        USAGE, {argv + 1, argv + argc},
        true,                                    // show help if requested
        "geo-social influence maximzation 0.1"); // version string

    network g = get_network(args);

    if (args.at("--export-network").asBool()) {
        std::cout << "Exporting network" << std::endl;

        auto graphfile = get_logfile("graph");
        write_graphviz(graphfile, g);
    }

    if (args.at("--statistics").asBool()) {
        write_node_degrees(g);
        std::cout << "Average node degree: " << get_average_node_degree(g)
                  << std::endl;

        std::cout << "Graph has " << boost::num_vertices(g) << " vertices"
                  << " and " << boost::num_edges(g) << " edges" << std::endl;
    }

    // Set budgets
    auto number_of_colors = get_number_of_colors(g);
    int budget_per_color = args.at("--budget").asLong();
    std::vector<unsigned int> budgets(number_of_colors, budget_per_color);

    // Get algorithm
    auto algorithm_name = args.at("--algorithm").asString();
    auto dataset_name = args.at("--dataset-name").asString();
    auto edge_weights = args.at("--edge-weights").asString();
    if (args.at("--random-weights").asBool()) {
        edge_weights = "random";
    }

    auto start_time = boost::chrono::process_user_cpu_clock::now();

    if (algorithm_name == "classic") {
        auto algorithm = algorithms::lazy_greedy(g);

        if (args.at("--statistics").asBool()) {
            algorithm.enable_generate_statistics();
        }
        if (args.at("--no-statusline").asBool()) {
            algorithm.disable_statusline();
        }

        algorithm.set_number_of_mc(args.at("--num-mc").asLong());

        std::pair<std::vector<double>,
                  std::unordered_map<vertex_descriptor, int>>
            seedset;
        if (args.at("--algo-version").asLong() == 0) {
            seedset = algorithm.maximize_influence(budgets);
        } else {
            seedset = algorithm.maximize_influence_baseline(budgets);
        }
        print_seedset(seedset);

        // Create logfile with influences and parameters of the program
        auto influence_file = get_logfile("influences-classic");
        for (int i{0}; i < seedset.first.size(); i++) {
            influence_file << budget_per_color << "\t" << number_of_colors
                           << "\t" << args.at("--num-mc").asLong() << "\t"
                           << edge_weights << "\t" << dataset_name << "\t"
                           << args.at("--algo-version").asLong() << "\t" << i
                           << "\t" << seedset.first[i] << "\n";
        }

        // Log time spent
        auto time_file = get_logfile("time-classic");
        auto runtime =
            boost::chrono::round<boost::chrono::microseconds>(
                boost::chrono::process_user_cpu_clock::now() - start_time)
                .count();
        time_file << budget_per_color << "\t" << number_of_colors << "\t"
                  << args.at("--num-mc").asLong() << "\t" << edge_weights
                  << "\t" << dataset_name << "\t"
                  << args.at("--algo-version").asLong() << "\t" << runtime
                  << "\n";
    } else if (algorithm_name == "ris") {
        auto algorithm = algorithms::ris(g);
        algorithm.use_version(args.at("--algo-version").asLong());

        double epsilon = std::stod(args.at("--epsilon").asString());
        double delta = std::stod(args.at("--delta").asString());

        auto seedset = algorithm.maximize_influence(budgets, epsilon, delta);
        print_seedset(seedset);

        // Create logfile with influences and parameters of the program
        auto influence_file = get_logfile("influences-ris");
        for (int i{0}; i < seedset.first.size(); i++) {
            influence_file << budget_per_color << "\t" << number_of_colors
                           << "\t" << edge_weights << "\t" << epsilon << "\t"
                           << delta << "\t" << dataset_name << "\t"
                           << args.at("--algo-version").asLong() << "\t" << i
                           << "\t" << seedset.first[i] << "\n";
        }

        // Log time spent
        auto time_file = get_logfile("time-ris");
        auto runtime =
            boost::chrono::round<boost::chrono::microseconds>(
                boost::chrono::process_user_cpu_clock::now() - start_time)
                .count();
        time_file << budget_per_color << "\t" << number_of_colors << "\t"
                  << edge_weights << "\t" << epsilon << "\t" << delta << "\t"
                  << dataset_name << "\t" << args.at("--algo-version").asLong()
                  << "\t" << runtime << "\n";
    }

    return 0;
}
