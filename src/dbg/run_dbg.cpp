#include <iostream>
#include <cstdlib>
#include "dot_graph.hpp"
#include "cmdline.h"
#include <filesystem>
#include "run_dbg.hpp"

using namespace dbg;

void DBGRunner::cleanDBG() {
    std::string restart_from = "";

    unsigned removed_edges = 1;

    std::cout << "----------Clean DBG----------" << std::endl;
    Graph graph;
    graph.read_graph(output, restart_from, graph_dot, graph_fasta, graph_aln);

    // std::cout << "----------Stage 1: Remove low-coverage edges ----------" << std::endl;
    // // removed_edges = 1;
    // // while (removed_edges) {
    // //     graph.remove_low_coverage_edges(removed_edges, graph.error_peak);
    // //     std::cout << "Removed " << removed_edges << " low-coverage edges" << std::endl;
    // // }
    // removed_edges = 1;
    // // 10 for hist figure, auto for results
    // while (removed_edges) {
    //     graph.remove_low_coverage_edges(removed_edges, 10, true);
    //     std::cout << "Removed " << removed_edges << " low-coverage tips" << std::endl;
    // }
    // this->first_peak = graph.error_peak;
    // this->first_minima = graph.first_minima;
    // graph.write_graph(output + "/graph.remove_low");
    // graph.write_graph_gfa(output + "/graph.remove_low");

    // std::unordered_set<std::string> nodes_linked_to_low;
    // for (auto&& node : graph.graph) {
    //     for (auto&& n_out : node.second.outgoing_edges) {
    //         for (auto&& e : n_out.second) {
    //             if (e.multiplicity <= 10) {
    //                 nodes_linked_to_low.insert(node.first);
    //                 nodes_linked_to_low.insert(n_out.first);
    //             }
    //         }
    //     }
    // }

    // graph.get_annotation(output + "/graph.remove_low");
    // graph.write_graph(output + "/graph.remove_low_only", 1000000, false, true, nodes_linked_to_low);

    // // std::cout << "----------Stage 2: Remove chimeric edges ----------" << std::endl;
    // // graph.detect_chimeric_reads();
    // // removed_edges = 1;
    // // while (removed_edges) {
    // //     graph.remove_low_coverage_edges(removed_edges, 0);
    // //     std::cout << "Removed " << removed_edges << " chimeric edges" << std::endl;
    // // }
    // // graph.write_graph(output + "/graph.cleaned");
    // // graph.write_graph_gfa(output + "/graph.cleaned");

    removed_edges = 1;
    while (removed_edges) {
        graph.multi_bulge_removal(removed_edges);
        std::cout << "Removed " << removed_edges << " simple bulges" << std::endl;
    }
    graph.write_graph(output + "/graph.simple_bulge");
    graph.get_annotation(output + "/graph.simple_bulge");
    graph.write_graph_contracted(output + "/graph.simple_bulge.20k", 20000);
    std::cout << "----------Clean DBG finished----------" << std::endl;
}
