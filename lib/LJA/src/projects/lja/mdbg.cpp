//
// Created by anton on 17.07.2020.
//
#define _GLIBCXX_PARALLEL
#include "error_correction/diploidy_analysis.hpp"
#include "error_correction/mult_correction.hpp"
#include "error_correction/parameter_estimator.hpp"
#include "dbg/visualization.hpp"
#include "dbg/graph_algorithms.hpp"
#include "dbg/dbg_construction.hpp"
#include "dbg/dbg_disjointigs.hpp"
#include "dbg/minimizer_selection.hpp"
#include "dbg/sparse_dbg.hpp"
#include "common/rolling_hash.hpp"
#include "common/hash_utils.hpp"
#include "common/hash_utils.hpp"
#include "error_correction/tournament_correction.hpp"
#include "sequences/seqio.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include "../dbg/graph_printing.hpp"
#include "dbg/dbg_graph_aligner.hpp"
#include <iostream>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <wait.h>
#include <common/id_index.hpp>
#include "repeat_resolution/repeat_resolution.hpp"

using namespace dbg;

// copied directly from repeat_resolution/mdbg_stage.hpp
std::unordered_map<std::string, std::experimental::filesystem::path> MDBGConstruction(
    logging::Logger& logger, size_t threads, size_t k, size_t kmdbg, size_t unique_threshold, bool diploid,
    const std::experimental::filesystem::path& dir,
    const std::experimental::filesystem::path& graph_gfa,
    const std::experimental::filesystem::path& read_paths, bool debug) {
    logger.info() << "Performing repeat resolution by transforming de Bruijn graph into Multiplex de Bruijn graph" << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1)
            << " to make it odd" << std::endl;
        k += 1;
    }
    hashing::RollingHash hasher(k);
    logger.info() << "Starting LoadDBGFromEdgeSequences" << std::endl;
    SparseDBG dbg = dbg::LoadDBGFromEdgeSequences(logger, threads, { graph_gfa }, hasher);
    logger.info() << "Done" << std::endl;

    logger.info() << "Initializing index" << std::endl;
    IdIndex<Vertex> index(dbg.vertices().begin(), dbg.vertices().end());
    logger.info() << "Done" << std::endl;

    size_t extension_size = 10000000;
    logger.info() << "Initializing readStorage" << std::endl;
    dbg::ReadAlignmentStorage readStorage(dbg, 0, extension_size, true, debug);
    logger.info() << "Done" << std::endl;
    logger.info() << "Initializing extra_reads" << std::endl;
    dbg::ReadAlignmentStorage extra_reads(dbg, 0, extension_size, false, debug);
    logger.info() << "Done" << std::endl;
    logger.info() << "loading all reads" << std::endl;
    ag::LoadAllReads<DBGTraits>(read_paths, { &readStorage, &extra_reads }, index);
    logger.info() << "Done" << std::endl;
    for (Vertex& v : dbg.verticesUnique()) {
        if (v.inDeg() == 1 && v.outDeg() == 1) {
            VERIFY(v.back() == v.rc().back().rc() || (v.back() == v.back().rc() && v.rc().back() == v.rc().back().rc()));
        }
    }
    repeat_resolution::RepeatResolver rr(dbg, &readStorage, { &extra_reads },
        k, kmdbg, dir, unique_threshold,
        diploid, debug, logger);
    return rr.ResolveRepeats(logger, threads);
}


std::string constructMessage() {
    std::stringstream ss;
    ss << "multiplexDBG\n";
    ss << "Usage: multiplexDBG [options] -o <output-dir> -i <mowerdbg-dir> -k <int>\n\n";
    ss << "Options:\n";
    ss << "  -o <file_name> (or --output-dir <file_name>)  Name of output folder. multiplexDBG outputs will be stored here.\n";
    ss << "  -g <file_name> (or --input-gfa <file_name>)   Path to input folder containing mowerDBG outputs (i.e. 01_TopologyBasedCorrection).\n";
    ss << "  -a <file_name> (or --input-aln <file_name>)   Path to input folder containing mowerDBG outputs (i.e. 01_TopologyBasedCorrection).\n";
    ss << "  -k <int>                                      Big k-mer size that was used in the final mowerDBG step (probably 5001).\n";
    ss << "  --max-k <int>                                 Max k. Default 40000.\n";
    ss << "  --unique-threshold <int>                      (From Marcus: I don't actually know what this is.) Default 40000.\n";
    ss << "  --diploid                                     Indicates whether the genome is diploid. We just pass this directly to multiplexDBG.\n";
    ss << "  -h (or --help)                                Print this help message.\n";
    ss << "  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.\n";
    return ss.str();
}

int main(int argc, char** argv) {
    AlgorithmParameters parameters({ "vertices=none", "unique=none", "coverages=none", "dbg=none", "output-dir=", "input-gfa=", "input-aln=",
                                    "max-k=40000", "unique-threshold=40000",
                                   "threads=16", "k-mer-size=", "base=239", "debug", "disjointigs=none", "reference=none",
                                   "simplify", "coverage", "cov-threshold=2", "rel-threshold=10", "tip-correct",
                                   "initial-correct", "mult-correct", "mult-analyse", "compress", "dimer-compress=1000000000,1000000000,1", "help", "genome-path",
                                   "dump", "extension-size=none", "print-all", "extract-subdatasets", "print-alignments", "subdataset-radius=10000",
    "split", "diploid" }, { "reads", "pseudo-reads", "align", "paths", "print-segment" }, constructMessage());
    CLParser parser(parameters, { "h=help", "o=output-dir", "t=threads", "k=k-mer-size","g=input-gfa","a=input-aln" });
    AlgorithmParameterValues params = parser.parseCL(argc, argv);
    if (params.getCheck("help")) {
        std::cout << params.helpMessage() << std::endl;
        return 0;
    }
    if (!params.checkMissingValues().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << params.checkMissingValues() << "\n" << std::endl;
        std::cout << params.helpMessage() << std::endl;
        return 1;
    }

    bool debug = params.getCheck("debug");
    bool diploid = params.getCheck("diploid");
    StringContig::homopolymer_compressing = params.getCheck("compress");
    StringContig::SetDimerParameters(params.getValue("dimer-compress"));
    const std::experimental::filesystem::path dir(params.getValue("output-dir"));
    const std::experimental::filesystem::path gfa(params.getValue("input-gfa"));
    const std::experimental::filesystem::path aln(params.getValue("input-aln"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for (size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    size_t k = std::stoi(params.getValue("k-mer-size"));
    size_t max_k = std::stoi(params.getValue("max-k"));
    size_t unique_threshold = std::stoi(params.getValue("unique-threshold"));
    logger << std::endl;
    logger.info() << "Hello! You are running multiplexDBG.\n";
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    size_t threads = std::stoi(params.getValue("threads"));
    omp_set_num_threads(threads);

    MDBGConstruction(
        logger, threads, k, max_k, unique_threshold, diploid,
        dir, gfa, aln, debug
    );

    logger.info() << "multiplexDBG finished." << std::endl;
    return 0;
}