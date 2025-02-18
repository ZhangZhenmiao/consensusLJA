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
    if (system((exeDir + "/get_reference.py -o " + prefix + ".ref.bam.stats " + prefix + ".ref.bam " + prefix + ".fasta").c_str()) != 0)
        exit(1);
    graph.write_graph_colored_from_bam(output + "/graph.bulge_removel.color", prefix + ".ref.bam.stats");
    graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.10k");
    graph.write_graph_contracted(output + "/graph.bulge_removel.contracted.20k", 20000);

    removed_paths = 1;
    total_removed = 0;
    while (removed_paths) {
        removed_bulges = 1;
        while (removed_bulges) {
            graph.multi_bulge_removal(removed_bulges);
        }
        graph.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 4);
        total_removed += removed_paths;
    }
    std::cout << "Removed " << total_removed << " complex bulges" << std::endl;

    graph.write_graph(output + "/graph.detouring_1_5");

    prefix = output + "/graph.detouring_1_5";
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
    if (system((exeDir + "/get_reference.py -o " + prefix + ".ref.bam.stats " + prefix + ".ref.bam " + prefix + ".fasta").c_str()) != 0)
        exit(1);
    graph.write_graph_colored_from_bam(output + "/graph.detouring_1_5.color", prefix + ".ref.bam.stats");
    graph.write_graph_contracted(output + "/graph.detouring_1_5.contracted.color.10k");
    graph.write_graph_contracted(output + "/graph.detouring_1_5.contracted.color.20k", 20000);
    return 0;
}
