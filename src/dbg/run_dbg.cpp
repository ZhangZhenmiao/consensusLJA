#include <iostream>
#include <cstdlib>
#include "dot_graph.hpp"
#include "cmdline.h"
#include <filesystem>
#include "run_dbg.hpp"

using namespace dbg;

void DBGRunner::cleanDBG() {
    std::string graph_dot = dot;
    std::string graph_fasta = fasta;
    std::string graph_aln = aln;
    double coverage = low;
    std::string restart_from = "";
    std::string output = output;

    unsigned removed_edges = 1;

    std::cout << "----------Clean DBG----------" << std::endl;
    Graph graph;
    graph.read_graph(output, restart_from, graph_dot, graph_fasta, graph_aln);

    std::cout << "----------Stage 1: Remove low-coverage edges ----------" << std::endl;
    removed_edges = 1;
    while (removed_edges) {
        graph.remove_low_coverage_edges(removed_edges, coverage);
        std::cout << "Removed " << removed_edges << " low-coverage edges" << std::endl;
    }
    graph.write_graph(output + "/graph.remove_low");
    graph.write_graph_gfa(output + "/graph.remove_low");

    std::cout << "----------Stage 2: Remove chimeric edges ----------" << std::endl;
    graph.detect_chimeric_reads();
    removed_edges = 1;
    while (removed_edges) {
        graph.remove_low_coverage_edges(removed_edges, coverage);
        std::cout << "Removed " << removed_edges << " low-coverage edges" << std::endl;
    }
    graph.write_graph(output + "/graph.remove_chimeric");
    graph.write_graph_gfa(output + "/graph.remove_chimeric");
    std::cout << "----------Clean DBG finished----------" << std::endl;
}
