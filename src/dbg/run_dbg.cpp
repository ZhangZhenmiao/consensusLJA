#include <iostream>
#include <cstdlib>
#include "dot_graph.hpp"
#include "cmdline.h"
#include <filesystem>
#include "run_dbg.hpp"
#include "utils.hpp"

using namespace dbg;

bool DBGRunner::correctHigh() {
    std::cout << "----------Correct reads in high-coverage regions----------" << std::endl;
    Graph graph;
    std::string restart_from = "";
    unsigned removed_edges = 1;

    graph.read_graph(output, restart_from, graph_dot, graph_fasta, graph_aln);
    int max_read_length = graph.write_reads(output + "/corrected_reads");
    std::cout << "Max compressed read: " << max_read_length << std::endl;
    this->mean_cov = graph.mean_cov;
    graph.pause_rerouting_reads = true;

    graph.write_graph(output + "/graph.ori");
    graph.get_annotation(output + "/graph.ori");
    graph.write_graph_gfa(output + "/graph.ori");

    graph.remove_low_coverage_edges(removed_edges, mean_cov * 10, false, true);

    unsigned removed_bulges = 1;
    unsigned removed_whirls = 1;
    unsigned removed_tips = 1;

    while (true) {
        bool flag = true;
        removed_bulges = 1;
        while (removed_bulges) {
            graph.multi_bulge_removal(removed_bulges);
            if (removed_bulges)
                flag = false;
        }
        removed_tips = 1;
        while (removed_tips) {
            graph.merge_tips_into_edges(removed_tips);
            if (removed_tips)
                flag = false;
        }
        removed_whirls = 1;
        while (removed_whirls) {
            graph.merge_tips_into_edges(removed_tips);
            graph.general_whirl_removal(removed_whirls);
            if (removed_whirls)
                flag = false;
        }
        if (flag)
            break;
    }

    graph.write_graph(output + "/graph.ori.only_high", 1000000, false, true);
    graph.append_linear_to_circular_genome(output + "/graph.ori.only_high", max_read_length);
    graph.get_annotation(output + "/graph.ori.only_high");

    if (graph.get_num_nodes() == 0)
        return false;
    else
        return true;
}

void DBGRunner::cleanDBG() {
    std::string restart_from = "";

    unsigned removed_edges = 1;
    unsigned removed_bulges = 1;
    unsigned removed_tips = 1;

    std::cout << "----------Clean DBG----------" << std::endl;
    Graph graph;
    graph.read_graph(output, restart_from, graph_dot, graph_fasta, graph_aln);
    this->first_peak = graph.error_peak;
    this->first_minima = graph.first_minima;
    this->mean_cov = graph.mean_cov;

    graph.write_graph(output + "/graph.ori");
    graph.get_annotation(output + "/graph.ori");

    std::cout << "----------Stage 1: Remove low-coverage edges ----------" << std::endl;
    removed_edges = 1;
    while (removed_edges) {
        graph.remove_low_coverage_edges(removed_edges, graph.first_minima, true);
        std::cout << "Removed " << removed_edges << " low-coverage tips" << std::endl;
    }
    removed_edges = 1;
    while (removed_edges) {
        graph.remove_low_coverage_edges(removed_edges, graph.first_minima, false);
        std::cout << "Removed " << removed_edges << " low-coverage edges" << std::endl;
    }
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

    graph.write_graph(output + "/graph.remove_chimeric");
    graph.get_annotation(output + "/graph.remove_chimeric");
    graph.write_graph_gfa(output + "/graph.remove_chimeric");

    graph.gluing_broken_bulges(removed_bulges);

    graph.write_graph(output + "/graph.cleaned");
    graph.get_annotation(output + "/graph.cleaned");
    graph.write_graph_gfa(output + "/graph.cleaned");

    std::cout << "----------Clean DBG finished----------" << std::endl;
}
