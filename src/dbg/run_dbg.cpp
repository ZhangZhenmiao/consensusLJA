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
    this->first_peak = graph.error_peak;
    this->first_minima = graph.first_minima;
    this->mean_cov = graph.mean_cov;

    std::cout << "----------Stage 1: Remove low-coverage edges ----------" << std::endl;
    removed_edges = 1;
    while (removed_edges) {
        graph.remove_low_coverage_edges(removed_edges, graph.first_minima, true);
        std::cout << "Removed " << removed_edges << " low-coverage tips" << std::endl;
    }
    // removed_edges = 1;
    // while (removed_edges) {
    //     graph.remove_low_coverage_edges(removed_edges, graph.first_minima, false);
    //     std::cout << "Removed " << removed_edges << " low-coverage edges" << std::endl;
    // }
    this->first_peak = graph.error_peak;
    this->first_minima = graph.first_minima;
    graph.write_graph(output + "/graph.remove_low");
    graph.write_graph_gfa(output + "/graph.remove_low");

    std::cout << "----------Stage 2: Remove chimeric edges ----------" << std::endl;
    graph.detect_chimeric_reads();
    removed_edges = 1;
    while (removed_edges) {
        graph.remove_low_coverage_edges(removed_edges, 0, false, true);
        std::cout << "Removed " << removed_edges << " 0-coverage edges" << std::endl;
    }
    graph.write_graph(output + "/graph.cleaned");
    graph.write_graph_gfa(output + "/graph.cleaned");
    std::cout << "----------Clean DBG finished----------" << std::endl;
}
