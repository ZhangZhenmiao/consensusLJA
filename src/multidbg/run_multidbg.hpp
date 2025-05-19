#include <string>
#include <filesystem>
#include <utils.hpp>
namespace fs = std::filesystem;

class MDBGRunner {
public:
    //input
    int threads = 0;
    //input for multidbg
    std::string graph_gfa, graph_aln;
    //input for simplification
    std::string graph_dot, graph_fasta, nodes_fasta, graph_dbg, paths_dbg;
    std::string reads;

    //output
    std::string mdbg_dir, output;
    //third-party
    std::string mdbg = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "multiDBG";
    std::string align_and_print = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "align_and_print";
    std::string compress = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "compress";
    std::string remove_chimeric = fs::path(getExecutablePath()).parent_path() / "src" / "scripts" / "remove_chimeric.sh";
    std::string analyze_chimeric = fs::path(getExecutablePath()).parent_path() / "src" / "scripts" / "analyze_chimeric.py";

    MDBGRunner(fs::path output_all, std::string reads, int threads) {
        this->threads = threads;
        this->reads = reads;
        fs::path mdbg_dir = output_all / "2_multiplex_DBG";
        if (!fs::is_directory(mdbg_dir)) {
            bool flag = true;
            int cnt = 0;
            while (flag) {
                if (execute_command(mdbg + " -g " + (output_all / "1_clean_DBG" / "graph.cleaned.gfa").string() + " -a " + (output_all / "1_clean_DBG" / "graph.cleaned.aln").string() + " -t " + std::to_string(threads) + " -k 5001 -o " + mdbg_dir.string() + " --diploid") == 0)
                    flag = false;
                if (++cnt == 10)
                    throw std::runtime_error("Failed to execute multidbg");
            }
        }
        fs::path align_dbg = mdbg_dir / "mdbg.align";
        if (!fs::is_directory(align_dbg)) {
            bool flag = true;
            int cnt = 0;
            while (flag) {
                if (execute_command(align_and_print + " --dbg " + (output_all / "1_clean_DBG" / "graph.cleaned.gfa").string() + " --paths " + (mdbg_dir / "mdbg_edge_seqs.fasta").string() + " --k-mer-size 5001 --output-dir " + align_dbg.string()) == 0)
                    flag = false;
                if (++cnt == 10)
                    throw std::runtime_error("Failed to execute align_and_print");
            }
        }

        graph_dot = mdbg_dir / "mdbg.hpc.dot";
        graph_fasta = mdbg_dir / "mdbg_edge_seqs.fasta";;
        nodes_fasta = mdbg_dir / "mdbg_vertex_seqs.fasta";
        output = output_all / "3_simplify_multiDBG";
        graph_dbg = output_all / "1_clean_DBG" / "graph.cleaned.dot";
        paths_dbg = align_dbg / "alignments.txt";

        simplifyMDBG();
    }
    void simplifyMDBG();
};