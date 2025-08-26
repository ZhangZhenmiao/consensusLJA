#include <iostream>
#include <dbg/run_dbg.hpp>
#include <multidbg/run_multidbg.hpp>
#include <filesystem>
#include "cmdline.h"
#include <iostream>
#include <cstdlib>
#include "multidbg/dot_graph.hpp"
#include "cmdline.h"
#include <filesystem>
#include <unistd.h>
#include "dbg/dot_graph.hpp"

namespace fs = std::filesystem;

using namespace multidbg;

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("reads", 'r', "path to reads", true);
    argParser.add<std::string>("LJA", 'l', "path to the LJA directory", true);
    argParser.add<std::string>("output", 'o', "the output directory", true);
    argParser.add<int>("threads", 't', "number of threads", false, 50);
    argParser.parse_check(argc, argv);

    //output
    std::string output = argParser.get<std::string>("output");

    fs::path output_all = fs::path(output);
    if (!fs::is_directory(output_all)) {
        fs::create_directory(output_all);
    }

    //input
    fs::path lja_path = fs::path(argParser.get<std::string>("LJA"));
    fs::path mdbg_dir = lja_path / "02_MDBG";
    fs::path dbg_dir = lja_path / "01_TopologyBasedCorrection";

    //third-party
    std::string mdbg = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "multiDBG";
    std::string align_and_print = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "align_and_print";
    std::string compress = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "compress";
    std::string remove_chimeric = fs::path(getExecutablePath()).parent_path() / "src" / "scripts" / "remove_chimeric.sh";
    std::string analyze_chimeric = fs::path(getExecutablePath()).parent_path() / "src" / "scripts" / "remove_chimeric.py";
    std::string jumbodbg = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "jumboDBG";
    std::string polisher = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "run_polishing";

    int threads = argParser.get<int>("threads");
    //input for simplification
    std::string graph_dot, graph_fasta, nodes_fasta, graph_dbg, paths_dbg;
    std::string reads = argParser.get<std::string>("reads");

    fs::path align_dbg = output_all / "mdbg.align";
    if (!fs::is_directory(align_dbg)) {
        bool flag = true;
        int cnt = 0;
        while (flag) {
            if (execute_command(align_and_print + " --dbg " + (dbg_dir / "final_dbg.gfa").string() + " --paths " + (mdbg_dir / "mdbg_edge_seqs.fasta").string() + " --k-mer-size 5001 --output-dir " + align_dbg.string()) == 0)
                flag = false;
            if (++cnt == 10)
                throw std::runtime_error("Failed to execute align_and_print");
        }
    }

    graph_dot = mdbg_dir / "mdbg.hpc.dot";
    graph_fasta = mdbg_dir / "mdbg_edge_seqs.fasta";;
    nodes_fasta = mdbg_dir / "mdbg_vertex_seqs.fasta";
    std::string output_connected = output_all / "connecting_using_dbg_and_spanning_reads";
    graph_dbg = dbg_dir / "final_dbg.dot";
    paths_dbg = align_dbg / "alignments.txt";

    fs::path gfa_path = fs::path(output_connected) / "graph.final.gfa";
    if (!fs::is_regular_file(gfa_path)) {
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
        dbg::Graph g;
        std::string dbg_fasta = graph_dbg.substr(0, graph_dbg.rfind(".")) + ".fasta";
        g.read_from_dot(graph_dbg, dbg_fasta);

        Graph graph;
        graph.read_graph(output_connected, graph_dot, graph_fasta, nodes_fasta, graph_dbg, paths_dbg);
        graph.write_graph(output_connected + "/graph.ori");
        graph.write_graph_gfa(output_connected + "/graph.ori");

        graph.mean_cov = g.mean_cov;
        graph.first_minima = g.first_minima;
        graph.error_peak = g.error_peak;

        removed_edges = 1;
        int cnt_round = 1;
        while (removed_edges) {
            graph.write_prefix_siffux_linear_edges(output_connected + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), jumbodbg, threads, 501, removed_edges);
            std::cout << "Glued " << removed_edges << " edges" << std::endl;
            graph.write_graph(output_connected + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), 1000000, false, true);

            cnt_round += 1;
        }

        removed_edges = 1;
        while (removed_edges) {
            graph.write_prefix_siffux_linear_edges(output_connected + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), jumbodbg, threads, 301, removed_edges);
            std::cout << "Glued " << removed_edges << " edges" << std::endl;
            graph.write_graph(output_connected + "/graph.glue_linear_edges_r" + std::to_string(cnt_round), 1000000, false, true);

            cnt_round += 1;
        }

        graph.connect_linear_and_tips_using_spanning_reads(output_connected + "/graph.spanning_reads", threads, reads, 0.95);

        graph.write_graph(output_connected + "/graph.final", 1000000, false, true);
        graph.write_graph_gfa(output_connected + "/graph.final");
    }

    // polishing
    fs::path corrected_reads_path = fs::path(output_connected) / "corrected_reads.paths";
    std::ofstream corrected_reads(corrected_reads_path);
    corrected_reads << (dbg_dir / "final_dbg.gfa").string() << std::endl;
    corrected_reads << (dbg_dir / "corrected_reads.aln").string() << std::endl;
    corrected_reads << 5001 << std::endl;
    corrected_reads.close();

    fs::path polisher_out = output_all / "polishing";

    if (execute_command(polisher + " --output-dir " + polisher_out.string() + " --graph " + gfa_path.string() + " --corrected_reads " + corrected_reads_path.string() + " --reads " + reads + " -t " + std::to_string(threads)) != 0)
        throw std::runtime_error("Failed to execute polisher: " + polisher + " --output-dir " + polisher_out.string() + " --graph " + gfa_path.string() + " --corrected_reads " + corrected_reads_path.string() + " --reads " + reads + " -t " + std::to_string(threads));
    return 0;
}
