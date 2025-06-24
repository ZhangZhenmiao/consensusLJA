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
    int first_peak = 0, first_minima = 0, mean_cov = 0;
    //third-party
    std::string lja = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "lja";
    std::string lja_cdb = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "lja_cdb";
    std::string compress = fs::path(getExecutablePath()).parent_path() / "lib" / "LJA" / "bin" / "compress";
    std::string correct_reads_script = fs::path(getExecutablePath()).parent_path() / "src" / "scripts" / "correct_reads_high.sh";
    std::string correct_reads_py = fs::path(getExecutablePath()).parent_path() / "src" / "scripts" / "correct_reads_high.py";

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

        this->output = output_all / "1_correct_high_reads";
        std::string high_contigs = output + "/graph.ori.only_high.fasta";
        std::string correct_reads_lja = output + "/corrected_reads.fasta";
        std::string correct_reads_high = output + "/reads_all.corrected.fasta";

        bool correct_high = false;
        if (!fs::is_regular_file(high_contigs))
            correct_high = correctHigh();

        dbg_dir = output_all / "1_correct_high_reads" / "condensed_dbg";
        this->dbg_dir = dbg_dir;
        if (correct_high) {
            if (!fs::is_regular_file(correct_reads_high)) {
                if (execute_command(correct_reads_script + " " + reads + " " + correct_reads_lja + " " + high_contigs + " " + output + "/reads_all" + " " + compress + " " + correct_reads_py) != 0) {
                    throw std::runtime_error("Failed to execute " + correct_reads_script);
                }
            }

            if (!fs::is_directory(dbg_dir)) {
                execute_command(lja_cdb + " -t " + std::to_string(threads) + " --reads " + correct_reads_high + " --output-dir " + dbg_dir.string() + " --diploid");
            }
        }
        else {
            if (!fs::is_regular_file(correct_reads_high)) {
                if (!fs::is_symlink(dbg_dir)) {
                    execute_command("ln -s `realpath " + (output_all / "0_condensed_dbg").string() + "` " + dbg_dir.string());
                }
            }
            else if (!fs::is_directory(dbg_dir / "01_TopologyBasedCorrection")) {
                execute_command(lja_cdb + " -t " + std::to_string(threads) + " --reads " + correct_reads_high + " --output-dir " + dbg_dir.string() + " --diploid");
            }
        }

        graph_dot = dbg_dir / "01_TopologyBasedCorrection" / "final_dbg.dot";
        graph_fasta = dbg_dir / "01_TopologyBasedCorrection" / "final_dbg.fasta";
        graph_aln = dbg_dir / "01_TopologyBasedCorrection" / "final_dbg.aln";

        this->output = output_all / "2_clean_DBG";

        if (!fs::is_directory(output))
            cleanDBG();
    }
    bool correctHigh();
    void cleanDBG();
};