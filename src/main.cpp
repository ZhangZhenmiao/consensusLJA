#include <iostream>
#include <cstdlib>
#include "dot_graph.hpp"
#include "cmdline.h"
#include <filesystem>
#include <unistd.h>

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("dot", 'd', "graph.dot file under LJA output", true);
    argParser.add<std::string>("fasta", 'f', "the graph.fasta file under LJA output", true);
    argParser.add<std::string>("graph_dbg", 'g', "the graph.dot file under 01_TopologyBasedCorrection for getting multiplicity", true);
    argParser.add<std::string>("paths_dbg", 'p', "the paths file for getting multiplicity", true);
    argParser.add<std::string>("nodes", 'n', "the mdbg_vertex_seqs.fasta file under LJA output", true);
    argParser.add<std::string>("output", 'o', "the output directory (should be new)", true);

    argParser.parse_check(argc, argv);
    std::string graph_dot = argParser.get<std::string>("dot");
    std::string graph_fasta = argParser.get<std::string>("fasta");
    std::string nodes_fasta = argParser.get<std::string>("nodes");
    std::string output = argParser.get<std::string>("output");
    std::string graph_dbg = argParser.get<std::string>("graph_dbg");
    std::string paths_dbg = argParser.get<std::string>("paths_dbg");

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
    graph.read_graph(output, graph_dot, graph_fasta, nodes_fasta, graph_dbg, paths_dbg);

    std::cout << "----------Stage 0: clean graph----------" << std::endl;
    graph.remove_low_coverage_edges(removed_edges);
    std::cout << "Removed " << removed_edges << " low-coverage edges" << std::endl;
    // graph.write_graph(output + "/graph.cleaned");
    // graph.get_annotation(output + "/graph.cleaned");
    // graph.write_graph_contracted(output + "/graph.cleaned.contracted.10k");
    // graph.write_graph_contracted(output + "/graph.cleaned.contracted.20k", 20000);

    std::cout << "----------Stage 1: simple bulge collapsing----------" << std::endl;
    removed_bulges = 1;
    total_removed = 0;
    while (removed_bulges) {
        graph.multi_bulge_removal(removed_bulges);
        total_removed += removed_bulges;
    }
    std::cout << "Removed " << total_removed << " simple bulges" << std::endl;
    graph.write_graph(output + "/graph.bulge_removel");
    // graph.get_annotation(output + "/graph.bulge_removel");
    // graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.10k");
    // graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.20k", 20000);

    std::cout << "----------Stage 2: whirl removal----------" << std::endl;
    removed_whirls = 1;
    total_removed = 0;
    while (removed_whirls) {
        graph.general_whirl_removal(removed_whirls);
        total_removed += removed_whirls;
    }
    std::cout << "Removed " << total_removed << " general whirls" << std::endl;

    std::cout << "----------Stage 3: N-M bulge collapsing----------" << std::endl;
    // secure N-M bulges, len=3, sim=0.8
    removed_paths = 1;
    total_removed = 0;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 3, 0.8, true, 3);
        total_removed += removed_paths;
    }
    // semi-secure N-M bulges, len=4, sim=0.8
    removed_paths = 1;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 4, 0.8, true, 2);
        total_removed += removed_paths;
    }
    // semi-secure N-M bulges, len=5, sim=0.8
    removed_paths = 1;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.8, true, 2);
        total_removed += removed_paths;
    }

    // semi-secure N-M bulges, len=5, sim=0
    removed_paths = 1;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.6, true, 2);
        total_removed += removed_paths;
    }
    std::cout << "Removed complex bulges: " << total_removed << std::endl;
    graph.write_graph(output + "/graph.complex_bulge");
    graph.get_annotation(output + "/graph.complex_bulge");
    graph.write_graph_contracted(output + "/graph.complex_bulge.contracted.10k");
    graph.write_graph_contracted(output + "/graph.complex_bulge.contracted.20k", 20000);

    std::cout << "----------Stage 4: merge tips into edges----------" << std::endl;
    removed_tips = 1;
    total_removed = 0;
    while (removed_tips) {
        graph.merge_tips_into_edges(removed_tips);
        total_removed += removed_tips;
    }

    std::cout << "Removed tips: " << total_removed << std::endl;
    graph.write_graph(output + "/graph.remove_tips");
    // graph.get_annotation(output + "/graph.remove_tips");
    // graph.write_graph_contracted(output + "/graph.remove_tips.contracted.10k");
    // graph.write_graph_contracted(output + "/graph.remove_tips.contracted.20k", 20000);

    std::cout << "----------Stage 5: decoupling strands----------" << std::endl;
    decoupled = 1;
    total_removed = 0;
    while (decoupled) {
        graph.resolve_edges_in_reverse_complement(decoupled);
        total_removed += decoupled;
    }

    std::cout << "Removed 2-in-2-out: " << total_removed << std::endl;
    graph.write_graph(output + "/graph.decoupling");
    graph.get_annotation(output + "/graph.decoupling");
    graph.write_graph_contracted(output + "/graph.decoupling.contracted.10k");
    graph.write_graph_contracted(output + "/graph.decoupling.contracted.20k", 20000);

    std::cout << "----------Stage 6: for complex components----------" << std::endl;

    removed_tips = 1;
    total_removed = 0;
    while (removed_tips) {
        graph.merge_tips_into_edges(removed_tips, 0.2);
        total_removed += removed_tips;
        removed_paths = 1;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0.6, true, 2);
            std::cout << "Removed complex bulges: " << removed_paths << std::endl;
        }
    }
    std::cout << "Removed tips: " << total_removed << std::endl;

    removed_whirls = 1;
    while (removed_whirls) {
        graph.general_whirl_removal(removed_whirls, false, true);
        graph.merge_non_branching_paths(true);
        std::cout << "Removed whirls: " << removed_whirls << std::endl;
    }

    decoupled = 1;
    while (decoupled) {
        graph.resolve_edges_in_reverse_complement(decoupled);
        std::cout << "Decoupled strands: " << decoupled << std::endl;
    }

    graph.write_graph(output + "/graph.complex_comp");
    graph.get_annotation(output + "/graph.complex_comp");
    graph.write_graph_contracted(output + "/graph.complex_comp.contracted.10k");
    graph.write_graph_contracted(output + "/graph.complex_comp.contracted.20k", 20000);

    return 0;
}
