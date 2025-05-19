#include <string>
#include <filesystem>
#include "utils.hpp"

namespace fs = std::filesystem;

class DBGRunner {
public:
    //input
    int threads = 0;
    std::string graph_dot, graph_fasta, graph_aln;
    //output
    std::string dbg_dir, output;
    int first_peak = 0, first_minima = 0;
    //third-party
    std::string lja = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "lja";

    DBGRunner(std::string reads, fs::path output_all, int threads) {
        this->threads = threads;

        fs::create_directory(output_all);
        fs::path dbg_dir = output_all / "0_condensed_dbg";
        this->dbg_dir = dbg_dir;

        if (!fs::is_directory(dbg_dir)) {
            execute_command(lja + " -t " + std::to_string(threads) + " --reads " + reads + " --output-dir " + dbg_dir.string() + " --diploid");
        }

        graph_dot = dbg_dir / "01_TopologyBasedCorrection" / "final_dbg.dot";
        graph_fasta = dbg_dir / "01_TopologyBasedCorrection" / "final_dbg.fasta";
        graph_aln = dbg_dir / "01_TopologyBasedCorrection" / "final_dbg.aln";

        this->output = output_all / "1_clean_DBG";

        if (!fs::is_directory(output))
            cleanDBG();
    }
    void cleanDBG();
};