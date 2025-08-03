#include <iostream>
#include <cstdlib>
#include "dot_graph.hpp"
#include "cmdline.h"
#include <filesystem>
#include <unistd.h>
#include "run_multidbg.hpp"
#include "dbg/dot_graph.hpp"

using namespace multidbg;

void MDBGRunner::simplifyMDBG() {
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

    // // Step 1 Read graph
    // std::cout << "----------Read graph----------" << std::endl;
    // dbg::Graph g;
    // std::string initial_fasta = initial_dbg.substr(0, initial_dbg.rfind(".")) + ".fasta";
    // g.read_from_dot(initial_dbg, initial_fasta);

    // this->mean_cov = g.mean_cov;
    // this->first_minima = g.first_minima;
    // this->first_peak = g.error_peak;

    // Graph graph;
    // graph.read_graph(output, graph_dot, graph_fasta, nodes_fasta, graph_dbg, paths_dbg);
    // graph.write_graph(output + "/graph.ori");
    // graph.write_graph_gfa(output + "/graph.ori");

    // // std::cout << "First minima: " << first_minima << "; average coverage: " << mean_cov << std::endl;

    // graph.mean_cov = mean_cov;
    // graph.first_minima = first_minima;
    // graph.error_peak = first_peak;

    // std::cout << "----------Stage 0: clean graph----------" << std::endl;
    // removed_edges = 1;
    // while (removed_edges) {
    //     graph.remove_low_coverage_edges(removed_edges, this->first_minima, true);
    //     std::cout << "Removed " << removed_edges << " low-coverage tips" << std::endl;
    // }
    // // removed_edges = 1;
    // // while (removed_edges) {
    // //     graph.remove_low_coverage_edges(removed_edges, this->first_minima, false);
    // //     std::cout << "Removed " << removed_edges << " low-coverage edges" << std::endl;
    // // }

    // graph.write_graph(output + "/graph.cleaned");
    // graph.write_graph_gfa(output + "/graph.cleaned");
    // // graph.get_annotation(output + "/graph.cleaned");
    // // graph.write_graph_contracted(output + "/graph.cleaned.contracted.10k");
    // graph.write_graph_contracted(output + "/graph.cleaned.contracted.20k", 20000);

    // std::string prefix = output + "/graph.cleaned";
    // execute_command(remove_chimeric + " " + reads + " " + prefix + " " + prefix + " " + compress + " " + analyze_chimeric);
    // graph.remove_chimeric_edge(prefix + ".chimeric.txt");

    // removed_tips = 1;
    // total_removed = 0;
    // while (removed_tips) {
    //     graph.merge_tips_into_edges(removed_tips, 0.8, true);
    //     graph.merge_non_branching_paths();
    //     total_removed += removed_tips;
    // }
    // std::cout << "Removed tips: " << total_removed << std::endl;
    // graph.write_graph(output + "/graph.chimeric_removed");

    // std::cout << "----------Stage 1: simple bulge collapsing----------" << std::endl;
    // removed_bulges = 1;
    // total_removed = 0;
    // while (removed_bulges) {
    //     graph.multi_bulge_removal(removed_bulges);
    //     total_removed += removed_bulges;
    // }
    // std::cout << "Removed " << total_removed << " simple bulges" << std::endl;
    // graph.write_graph(output + "/graph.bulge_removel");
    // // graph.get_annotation(output + "/graph.bulge_removel");
    // // graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.10k");
    // graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.20k", 20000);

    // std::cout << "----------Stage 2: whirl removal----------" << std::endl;
    // removed_whirls = 1;
    // total_removed = 0;
    // while (removed_whirls) {
    //     graph.general_whirl_removal(removed_whirls);
    //     total_removed += removed_whirls;
    // }
    // std::cout << "Removed " << total_removed << " general whirls" << std::endl;

    // std::cout << "----------Stage 3: N-M bulge collapsing----------" << std::endl;
    // // secure N-M bulges, len=3, sim=0.8
    // removed_paths = 1;
    // total_removed = 0;
    // while (removed_paths) {
    //     graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 3, 0.8, true, 3);
    //     total_removed += removed_paths;
    // }
    // // semi-secure N-M bulges, len=4, sim=0.8
    // removed_paths = 1;
    // while (removed_paths) {
    //     graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 4, 0.8, true, 2);
    //     total_removed += removed_paths;
    // }
    // // semi-secure N-M bulges, len=5, sim=0.8
    // removed_paths = 1;
    // while (removed_paths) {
    //     graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.8, true, 2);
    //     total_removed += removed_paths;
    // }

    // // semi-secure N-M bulges, len=5, sim=0
    // removed_paths = 1;
    // while (removed_paths) {
    //     graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.6, true, 2);
    //     total_removed += removed_paths;
    // }
    // std::cout << "Removed complex bulges: " << total_removed << std::endl;
    // graph.write_graph(output + "/graph.complex_bulge");
    // // graph.get_annotation(output + "/graph.complex_bulge");
    // // graph.write_graph_contracted(output + "/graph.complex_bulge.contracted.10k");
    // graph.write_graph_contracted(output + "/graph.complex_bulge.contracted.20k", 20000);

    // // ensure all below outputting graph have no simple bulges, or the program will fail
    // std::cout << "----------Stage 4: decoupling strands----------" << std::endl;
    // decoupled = 1;
    // total_removed = 0;
    // while (decoupled) {
    //     graph.resolve_edges_in_reverse_complement(decoupled);
    //     total_removed += decoupled;
    // }

    // std::cout << "Removed 2-in-2-out: " << total_removed << std::endl;
    // // removed_bulges = 1;
    // // while (removed_bulges) {
    // //     graph.merge_non_branching_paths(true);
    // //     graph.multi_bulge_removal(removed_bulges, false);
    // // }
    // graph.write_graph(output + "/graph.decoupling");
    // // graph.get_annotation(output + "/graph.decoupling");
    // graph.write_graph_contracted(output + "/graph.decoupling.contracted.20k", 20000);


    // std::cout << "----------Stage 5: merge tips into edges----------" << std::endl;

    // total_removed = 0;
    // removed_tips = 1;
    // while (removed_tips) {
    //     graph.merge_tips_into_edges(removed_tips);
    //     total_removed += removed_tips;
    //     if (removed_tips > 0)
    //         std::cout << "Merged " << removed_tips << " tips to edges" << std::endl;
    // }

    // std::cout << "Removed tips: " << total_removed << std::endl;
    // // removed_bulges = 1;
    // // while (removed_bulges) {
    // //     graph.merge_non_branching_paths(true);
    // //     graph.multi_bulge_removal(removed_bulges, false);
    // // }
    // graph.write_graph(output + "/graph.remove_tips");
    // // graph.get_annotation(output + "/graph.remove_tips");
    // // graph.write_graph_contracted(output + "/graph.remove_tips.contracted.10k");
    // graph.write_graph_contracted(output + "/graph.remove_tips.contracted.20k", 20000);

    // std::cout << "----------Stage 6: for complex components----------" << std::endl;

    // while (true) {
    //     bool flag = true;
    //     removed_paths = 1;
    //     while (removed_paths) {
    //         graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0.6, true, 2);
    //         if (removed_paths)
    //             flag = false;
    //         if (removed_paths > 0)
    //             std::cout << "Detoured " << removed_paths << " paths" << std::endl;
    //     }

    //     decoupled = 1;
    //     while (decoupled) {
    //         graph.resolve_edges_in_reverse_complement(decoupled);
    //         if (decoupled)
    //             flag = false;
    //         if (decoupled > 0)
    //             std::cout << "Decoupled " << decoupled << " strands" << std::endl;
    //     }

    //     removed_tips = 1;
    //     while (removed_tips) {
    //         graph.merge_tips_into_edges(removed_tips);
    //         if (removed_tips)
    //             flag = false;
    //         if (removed_tips > 0)
    //             std::cout << "Merged " << removed_tips << " tips to edges" << std::endl;
    //     }

    //     removed_whirls = 1;
    //     while (removed_whirls) {
    //         graph.general_whirl_removal(removed_whirls);
    //         graph.merge_non_branching_paths(true);
    //         if (removed_whirls)
    //             flag = false;
    //         if (removed_whirls > 0)
    //             std::cout << "Removed " << removed_whirls << " whirls" << std::endl;
    //     }

    //     removed_bulges = 1;
    //     while (removed_bulges) {
    //         graph.multi_bulge_removal(removed_bulges);
    //         graph.merge_non_branching_paths(true);
    //         if (removed_bulges)
    //             flag = false;
    //         if (removed_bulges > 0)
    //             std::cout << "Removed " << removed_bulges << " bulges" << std::endl;
    //     }

    //     if (flag)
    //         break;
    // }

    // graph.write_graph(output + "/graph.complex_comp");
    // // graph.get_annotation(output + "/graph.complex_comp");

    // std::cout << "----------Stage 7: decoupling further----------" << std::endl;
    // while (true)
    // {
    //     bool flag = true;
    //     removed_paths = 1;
    //     while (removed_paths) {
    //         graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0.6, true, 2, true);
    //         if (removed_paths)
    //             flag = false;
    //         if (removed_paths > 0)
    //             std::cout << "Detoured " << removed_paths << " paths" << std::endl;
    //     }

    //     decoupled = 1;
    //     while (decoupled) {
    //         graph.resolve_edges_in_reverse_complement(decoupled);
    //         if (decoupled)
    //             flag = false;
    //         if (decoupled > 0)
    //             std::cout << "Decoupled " << decoupled << " strands" << std::endl;
    //     }

    //     if (flag)
    //         break;
    // }

    // graph.write_graph(output + "/graph.decoupling_further");
    // // graph.get_annotation(output + "/graph.decoupling_further");

    // // for debugging
    // Graph graph;
    // graph.restart_from_dot(output + "/graph.decoupling_further.dot", output + "/graph.decoupling_further.fasta");
    // graph.write_graph_colored_from_bam(output + "/graph.decoupling_further" + ".color", output + "/graph.decoupling_further" + ".ref.bam.stats");
    // // for debugging

    // std::cout << "----------Stage 8: contract graph----------" << std::endl;
    // graph.write_graph_contracted(output + "/graph.complex_comp.contracted.20k", 20000);
    // graph.write_graph_contracted(output + "/graph.complex_comp_simplify.contracted.20k", 20000, true);

    // while (true) {
    //     bool flag = true;
    //     removed_paths = 1;
    //     while (removed_paths) {
    //         graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0, true, 2, true);
    //         if (removed_paths)
    //             flag = false;
    //         if (removed_paths > 0)
    //             std::cout << "Removed " << removed_paths << " paths" << std::endl;
    //     }
    //     decoupled = 1;
    //     while (decoupled) {
    //         graph.resolve_edges_in_reverse_complement(decoupled);
    //         if (decoupled)
    //             flag = false;
    //         if (decoupled > 0)
    //             std::cout << "Decoupled " << decoupled << " strands" << std::endl;
    //     }

    //     removed_tips = 1;
    //     while (removed_tips) {
    //         graph.merge_tips_into_edges(removed_tips, 0.8);
    //         if (removed_tips)
    //             flag = false;
    //         if (removed_tips > 0)
    //             std::cout << "Merged " << removed_tips << " tips to edges" << std::endl;
    //     }

    //     removed_whirls = 1;
    //     while (removed_whirls) {
    //         graph.general_whirl_removal(removed_whirls);
    //         graph.merge_non_branching_paths(true);
    //         if (removed_whirls)
    //             flag = false;
    //         if (removed_whirls > 0)
    //             std::cout << "Removed " << removed_whirls << " whirls" << std::endl;
    //     }

    //     removed_bulges = 1;
    //     while (removed_bulges) {
    //         graph.multi_bulge_removal(removed_bulges);
    //         graph.merge_non_branching_paths(true);
    //         if (removed_bulges)
    //             flag = false;
    //         if (removed_bulges > 0)
    //             std::cout << "Removed " << removed_bulges << " bulges" << std::endl;
    //     }

    //     if (flag)
    //         break;
    // }

    // graph.write_graph(output + "/graph.before_final", false, true);
    // graph.write_graph_colored_from_bam(output + "/graph.before_final" + ".color", output + "/graph.before_final" + ".ref.bam.stats");
    // // graph.get_annotation(output + "/graph.before_final");

    Graph graph;
    graph.restart_from_dot(output + "/graph.before_final.dot", output + "/graph.before_final.fasta");
    // graph.write_graph_colored_from_bam(output + "/graph.before_final" + ".color", output + "/graph.before_final" + ".ref.bam.stats");

    while (true) {
        bool flag = true;

        removed_paths = 1;
        while (removed_paths) {
            graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.9, true, 2, true, true);
            if (removed_paths)
                flag = false;
            if (removed_paths > 0)
                std::cout << "Removed " << removed_paths << " paths" << std::endl;
        }

        decoupled = 1;
        while (decoupled) {
            graph.resolve_edges_in_reverse_complement(decoupled);
            if (decoupled)
                flag = false;
            if (decoupled > 0)
                std::cout << "Decoupled " << decoupled << " strands" << std::endl;
        }

        decoupled = 1;
        while (decoupled) {
            graph.resolving_complex_palindromic_bulges(decoupled, 8);
            if (decoupled)
                flag = false;
            if (decoupled > 0)
                std::cout << "Removed " << decoupled << " complex palindromic bulges" << std::endl;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph.merge_tips_into_edges(removed_tips, 0.8, false, true);
            if (removed_tips)
                flag = false;
            if (removed_tips > 0)
                std::cout << "Merged " << removed_tips << " tips to edges" << std::endl;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph.merge_tips(removed_tips);
            if (removed_tips)
                flag = false;
            if (removed_tips > 0)
                std::cout << "Merged " << removed_tips << " tips to tips" << std::endl;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph.merge_tips_into_edges_further(removed_tips);
            if (removed_tips)
                flag = false;
            if (removed_tips > 0)
                std::cout << "Merged " << removed_tips << " tips to paths" << std::endl;
        }

        removed_whirls = 1;
        while (removed_whirls) {
            graph.general_whirl_removal(removed_whirls, false, true);
            graph.merge_non_branching_paths(true);
            if (removed_whirls)
                flag = false;
            if (removed_whirls > 0)
                std::cout << "Removed " << removed_whirls << " whirls" << std::endl;
        }

        removed_bulges = 1;
        while (removed_bulges) {
            graph.multi_bulge_removal(removed_bulges);
            graph.merge_non_branching_paths(true);
            if (removed_bulges)
                flag = false;
            if (removed_bulges > 0)
                std::cout << "Removed " << removed_bulges << " bulges" << std::endl;
        }

        removed_edges = 1;
        while (removed_edges) {
            graph.merge_secondary_edges(removed_edges);
            if (removed_edges)
                flag = false;
            if (removed_edges > 0)
                std::cout << "Merged " << removed_edges << " deadend edges" << std::endl;
        }

        if (flag)
            break;
    }

    graph.write_graph(output + "/graph.before_removing_contained", false, true);
    // graph.get_annotation(output + "/graph.before_removing_contained");

    // Graph graph;
    // graph.restart_from_dot(output + "/graph.before_removing_contained.dot", output + "/graph.before_removing_contained.fasta");
    // graph.write_graph_colored_from_bam(output + "/graph.before_removing_contained" + ".color", output + "/graph.before_removing_contained" + ".ref.bam.stats");

    removed_edges = 1;
    int cnt_round = 1;
    while (removed_edges) {
        graph.write_prefix_siffux_linear_edges(output + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), jumbodbg, threads, 501, removed_edges);
        std::cout << "Glued " << removed_edges << " edges" << std::endl;
        graph.write_graph(output + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), 1000000, false, true);

        removed_paths = 1;
        graph.remove_contained_contigs_minimap(output + "/graph.remove_contained_r" + std::to_string(cnt_round), threads, removed_paths);
        std::cout << "Removed " << removed_paths << " edges" << std::endl;
        graph.write_graph(output + "/graph.remove_contained_r" + std::to_string(cnt_round), 1000000, false, true);

        cnt_round += 1;
    }

    removed_edges = 1;
    while (removed_edges) {
        graph.write_prefix_siffux_linear_edges(output + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), jumbodbg, threads, 301, removed_edges);
        std::cout << "Glued " << removed_edges << " edges" << std::endl;
        graph.write_graph(output + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), 1000000, false, true);

        removed_paths = 1;
        graph.remove_contained_contigs_minimap(output + "/graph.remove_contained_r" + std::to_string(cnt_round), threads, removed_paths);
        std::cout << "Removed " << removed_paths << " edges" << std::endl;
        graph.write_graph(output + "/graph.remove_contained_r" + std::to_string(cnt_round), 1000000, false, true);

        cnt_round += 1;
    }

    graph.write_graph(output + "/graph.final", 1000000, false, true);
    graph.get_annotation(output + "/graph.final");

    graph.write_graph_gfa(output + "/graph.final");
}
