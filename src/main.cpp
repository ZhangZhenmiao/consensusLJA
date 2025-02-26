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

    graph.write_graph(output + "/original");
    removed_bulges = 1;
    cnt_rounds = 0;
    while (removed_bulges) {
        std::cout << "----------Stage 1: Simple bulge collapsing round " << ++cnt_rounds << "----------" << std::endl;
        graph.multi_bulge_removal(removed_bulges);
        std::cout << "Removed " << removed_bulges << " bulges" << std::endl;
    }
    graph.write_graph(output + "/graph.bulge_removel");
    std::string prefix = output + "/graph.bulge_removel";
    std::string ref_seq = "/Poppy/zmzhang/Rust_fungi/genome/reference.compressed.only_chrs.fasta";
    if (system(("minimap2 -ax asm20 " + ref_seq + " " + prefix + ".fasta -t 100 | grep -v '^@' > " + prefix + ".ref.sam").c_str()) != 0) {
        exit(1);
    }
    if (!std::filesystem::exists(ref_seq + ".fai")) {
        if (system(("samtools faidx " + ref_seq).c_str()) != 0)
            exit(1);
    }
    if (system(("cut -f1,2 " + ref_seq + ".fai | awk " + R"('{print "@SQ\tSN:"$1"\tLN:"$2}')" + " > " + prefix + ".ref.header.sam").c_str()) != 0)
        exit(1);
    if (system(("cat " + prefix + ".ref.header.sam " + prefix + ".ref.sam | samtools sort -@ 50 -o " + prefix + ".ref.bam").c_str()) != 0) {
        exit(1);
    }
    std::string exeDir = graph.getExecutablePath();
    if (system((exeDir + "/../src/scripts/get_reference.py -o " + prefix + ".ref.bam.stats " + prefix + ".ref.bam " + prefix + ".fasta").c_str()) != 0)
        exit(1);
    graph.write_graph_colored_from_bam(output + "/graph.bulge_removel.color", prefix + ".ref.bam.stats");
    graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.10k");
    graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.20k", 20000);

    removed_whirls = 1;
    total_removed = 0;
    while (removed_whirls) {
        graph.general_whirl_removal(removed_whirls);
        std::cout << "Removed " << removed_whirls << " general whirls" << std::endl;
        total_removed += removed_whirls;
    }
    std::cout << "Removed " << total_removed << " general whirls in total" << std::endl;

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
        total_removed += removed_paths;
    }
    graph.write_graph(output + "/graph.complex_bulge_stage3.3");
    std::cout << "Removed complex bulges: " << total_removed << std::endl;

    total_removed = 1;
    while (total_removed) {
        graph.resolve_edges_in_reverse_complement(total_removed, true);
    }
    graph.write_graph(output + "/graph.resolve_edges_in_reverse_complement_rc");

    // Step 5 Broken bulges and tips
    removed_bulges = 1;
    while (removed_bulges) {
        graph.gluing_broken_bulges(removed_bulges);
    }
    graph.write_graph(output + "/graph.tips_processed");

    removed_paths = 1;
    total_removed = 0;
    while (removed_paths) {
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 2);
        total_removed += removed_paths;
        graph.merge_non_branching_paths(true);
    }
    std::cout << "Removed complex bulges: " << total_removed << std::endl;

    graph.write_graph(output + "/graph.final");

    prefix = output + "/graph.final";
    ref_seq = "/Poppy/zmzhang/Rust_fungi/genome/reference.compressed.only_chrs.fasta";
    if (system(("minimap2 -ax asm20 " + ref_seq + " " + prefix + ".fasta -t 100 | grep -v '^@' > " + prefix + ".ref.sam").c_str()) != 0) {
        exit(1);
    }
    if (!std::filesystem::exists(ref_seq + ".fai")) {
        if (system(("samtools faidx " + ref_seq).c_str()) != 0)
            exit(1);
    }
    if (system(("cut -f1,2 " + ref_seq + ".fai | awk " + R"('{print "@SQ\tSN:"$1"\tLN:"$2}')" + " > " + prefix + ".ref.header.sam").c_str()) != 0)
        exit(1);
    if (system(("cat " + prefix + ".ref.header.sam " + prefix + ".ref.sam | samtools sort -@ 50 -o " + prefix + ".ref.bam").c_str()) != 0) {
        exit(1);
    }
    exeDir = graph.getExecutablePath();
    if (system((exeDir + "/../src/scripts/get_reference.py -o " + prefix + ".ref.bam.stats " + prefix + ".ref.bam " + prefix + ".fasta").c_str()) != 0)
        exit(1);
    graph.write_graph_colored_from_bam(output + "/graph.final.color", prefix + ".ref.bam.stats");
    graph.write_graph_contracted(output + "/graph.final.contracted.color.10k");
    graph.write_graph_contracted(output + "/graph.final.contracted.color.15k", 15000);
    graph.write_graph_contracted(output + "/graph.final.contracted.color.20k", 20000);
    return 0;
}
