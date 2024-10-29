#include <iostream>
#include <cstdlib>
#include "dot_graph.hpp"
#include "cmdline/cmdline.h"
#include <filesystem>

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("dot", 'd', "graph.dot file under LJA output", true);
    argParser.add<std::string>("fasta", 'f', "the graph.fasta file under LJA output", true);
    argParser.add<std::string>("nodes", 'n', "the mdbg_vertex_seqs.fasta file under LJA output", true);
    // argParser.add<std::string>("restart", 'r', "restart from", true);
    argParser.add<std::string>("output", 'o', "the output directory (should be new)", true);

    argParser.parse_check(argc, argv);
    std::string graph_dot = argParser.get<std::string>("dot");
    std::string graph_fasta = argParser.get<std::string>("fasta");
    std::string nodes_fasta = argParser.get<std::string>("nodes");
    std::string restart_from = ""; // TODO: enable program restart at middle
    std::string output = argParser.get<std::string>("output");

    // define count variables
    unsigned removed_paths = 1;
    unsigned removed_whirls = 1;
    unsigned removed_bulges = 1;
    int decoupled = 1;
    unsigned cnt_rounds = 0;
    int total_removed = 0;
    unsigned total_whirls = 0;
    unsigned removed_edges = 1;
    unsigned removed_tips = 1;

    // Step 1 Read graph
    std::cout << "----------Read graph----------" << std::endl;
    Graph graph;
    graph.read_graph(output, restart_from, graph_dot, graph_fasta, nodes_fasta);

    removed_bulges = 1;
    cnt_rounds = 0;
    while (removed_bulges) {
        std::cout << "----------Stage 1: Simple bulge collapsing round " << ++cnt_rounds << "----------" << std::endl;
        graph.multi_bulge_removal(removed_bulges);
        std::cout << "Removed " << removed_bulges << " bulges" << std::endl;
    }
    graph.write_graph(output + "/graph.bulge_removel");
    return 0;
}
