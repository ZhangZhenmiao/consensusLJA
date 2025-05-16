#include <iostream>
#include <dbg/run_dbg.hpp>
#include <multidbg/run_multidbg.hpp>
#include "cmdline.h"

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("dot", 'd', "graph.dot file under LJA output", true);
    argParser.add<std::string>("fasta", 'f', "the graph.fasta file under LJA output", true);
    argParser.add<std::string>("aln", 'a', "the graph.aln file under LJA output", true);
    argParser.add<std::string>("output", 'o', "the output directory (should be new)", true);

    argParser.parse_check(argc, argv);

    DBGRunner(argParser.get<std::string>("dot"), argParser.get<std::string>("fasta"), argParser.get<std::string>("aln"), argParser.get<std::string>("output"));
    return 0;
}
