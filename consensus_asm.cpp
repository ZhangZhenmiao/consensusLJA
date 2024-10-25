#include <iostream>
#include <cstdlib>
#include "dot_graph.hpp"
#include "cmdline/cmdline.h"
#include <filesystem>

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("dot", 'd', "graph.dot file under LJA output", true);
    argParser.add<std::string>("fasta", 'f', "the graph.fasta file under LJA output", true);
    argParser.add<std::string>("restart", 'r', "restart from", true);
    argParser.add<std::string>("multidbg", 'm', "the multidbg edges under LJA output", true);
    argParser.add<std::string>("output", 'o', "the output directory (should be new)", true);

    argParser.parse_check(argc, argv);
    std::string graph_dot = argParser.get<std::string>("dot");
    std::string graph_fasta = argParser.get<std::string>("fasta");
    std::string restart_from = argParser.get<std::string>("restart");; // TODO: enable program restart at middle
    std::string output = argParser.get<std::string>("output");
    std::string multidbg = argParser.get<std::string>("multidbg");

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
    graph.read_graph(output, restart_from, graph_dot, graph_fasta);

    if (restart_from != "graph.2_in_2_out_round1") {
        // Step 2 Process bulges and whirls
        removed_bulges = 1;
        cnt_rounds = 0;
        while (removed_bulges) {
            std::cout << "----------Stage 1: Simple bulge collapsing round " << ++cnt_rounds << "----------" << std::endl;
            graph.multi_bulge_removal(removed_bulges);
            std::cout << "Removed " << removed_bulges << " bulges" << std::endl;
        }
        graph.write_graph(output + "/graph.bulge_removel");

        removed_whirls = 1;
        cnt_rounds = 0;
        while (removed_whirls) {
            std::cout << "----------Stage 1: General whirl removal round " << ++cnt_rounds << "----------" << std::endl;
            graph.general_whirl_removal(removed_whirls);
        }
        graph.write_graph(output + "/graph.genral_whirls");

        // Step 3 Complex bulge collapsing
        removed_paths = 1;
        total_removed = 0;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 3, 0.8, true, 3);
            total_removed += removed_paths;
        }
        graph.write_graph(output + "/graph.complex_bulge_stage3.1");
        std::cout << "Removed complex bulges: " << total_removed << std::endl;

        // Step 3.2 Collapse paths < 4 edges, semi-secure, do not allow reverse complementary
        removed_paths = 1;
        total_removed = 0;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 4, 0.8, true, 2);
            total_removed += removed_paths;
        }
        graph.write_graph(output + "/graph.complex_bulge_stage3.2");
        std::cout << "Removed complex bulges: " << total_removed << std::endl;

        // Step 3.3 Collapse paths < 5 edges, semi-secure, do not allow reverse complementary
        removed_paths = 1;
        total_removed = 0;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.8, true, 2);
        }
        graph.write_graph(output + "/graph.complex_bulge_stage3.3");
        std::cout << "Removed complex bulges: " << total_removed << std::endl;

        // Step 4 Decouple with reverse complementary in nodes and out nodes
        decoupled = 1;
        while (decoupled) {
            graph.decoupling(decoupled, true);
        }
        graph.write_graph(output + "/graph.decoupling_rc");

        // Step 5 Broken bulges and tips
        removed_bulges = 1;
        while (removed_bulges) {
            graph.gluing_broken_bulges(removed_bulges);
        }
        graph.write_graph(output + "/graph.tips_processed");

        removed_edges = 1;
        while (removed_edges)
            graph.remove_low_coverage_edges(removed_edges);

        removed_paths = 1;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 2);
            graph.merge_non_branching_paths(true);
        }

        graph.write_graph(output + "/graph.tips_processed");

        // Step 6 Decouple with MultiDBG
        std::cout << "----------Stage 6: MultiDBG----------" << std::endl;
        // TODO: will this module to not depend on minimap2
        graph.resolve_2_in_2_out_mdbg(multidbg, output + "/graph.multidbg.round1");
        removed_paths = 1;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 2);
            graph.merge_non_branching_paths(true);
        }
        graph.write_graph(output + "/graph.2_in_2_out_round1");

        // TODO: the second round may not be necessary
        graph.resolve_2_in_2_out_mdbg(multidbg, output + "/graph.multidbg.round2");
        removed_paths = 1;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 2);
            graph.merge_non_branching_paths(true);
        }
        graph.write_graph(output + "/graph.2_in_2_out_round2");
    }

    // Step 7 Further simplifications without considering reverse complementary nodes
    std::cout << "----------Stage 7: Further simplifications (no RC)----------" << std::endl;
    // remove low coverage edges: cov <= 10
    removed_edges = 1;
    while (removed_edges)
        graph.remove_low_coverage_edges(removed_edges);

    removed_tips = 1;
    while (removed_tips) {
        graph.merge_tips_into_edges(removed_tips);
    }

    // remove low coverage edges: cov <= 15
    removed_edges = 1;
    while (removed_edges)
        graph.remove_low_coverage_edges(removed_edges, 15);
    removed_tips = 1;
    while (removed_tips) {
        graph.merge_tips_into_edges(removed_tips);
    }

    // complex bulges with at most 10 edges for each path
    removed_paths = 1;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 10, 0, true, 2);
        graph.merge_non_branching_paths(true);
    }

    decoupled = 1;
    while (decoupled) {
        graph.decoupling(decoupled, true);
    }

    graph.general_whirl_removal(cnt_rounds, false, true);

    graph.merge_tips(removed_tips);
    graph.write_graph(output + "/before_complex_rc");

    // Step 8 Further simplifications considering reverse complementary nodes
    std::cout << "----------Stage 8: Further simplifications (with RC)----------" << std::endl;
    decoupled = 1;
    while (decoupled) {
        graph.decoupling(decoupled);
    }

    removed_paths = 1;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 2, true);
        graph.merge_non_branching_paths(true);
    }

    removed_paths = 1;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 10, 0, true, 2, true);
        graph.merge_non_branching_paths(true);
    }

    decoupled = 1;
    while (decoupled) {
        graph.decoupling(decoupled);
    }

    removed_tips = 1;
    while (removed_tips) {
        graph.merge_tips_into_edges(removed_tips);
    }
    graph.merge_tips(removed_tips, false);

    decoupled = 1;
    while (decoupled) {
        graph.decoupling(decoupled);
        graph.merge_tips(removed_tips, false);
    }

    removed_paths = 1;
    while (removed_paths) {
        decoupled = 1;
        while (decoupled) {
            graph.decoupling(decoupled);
            graph.merge_tips(removed_tips, false);
            graph.merge_non_branching_paths(true);
        }
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 10, 0, true, 2, true);
        graph.merge_non_branching_paths(true);
    }

    removed_tips = 1;
    while (removed_tips) {
        decoupled = 1;
        while (decoupled) {
            graph.decoupling(decoupled);
        }
        graph.merge_tips(removed_tips, false);
        graph.merge_non_branching_paths(true);
        total_whirls = 1;
        while (total_whirls)
            graph.general_whirl_removal(total_whirls, false, true);
        graph.merge_non_branching_paths(true);

    }
    removed_tips = 1;
    while (removed_tips) {
        decoupled = 1;
        while (decoupled) {
            graph.decoupling(decoupled);
        }
        graph.merge_tips(removed_tips, true);
        graph.merge_non_branching_paths(true);
        total_whirls = 1;
        while (total_whirls)
            graph.general_whirl_removal(total_whirls, false, true);
        graph.merge_non_branching_paths(true);
    }

    graph.write_graph(output + "/graph.final");
    std::cout << "----------cLJA finished----------" << std::endl;
    return 0;
}
