#include <iostream>
#include <dbg/run_dbg.hpp>
#include <multidbg/run_multidbg.hpp>
#include "cmdline.h"

int main(int argc, char* argv[]) {
    cmdline::parser argParser;
    argParser.add<std::string>("reads", 'r', "path to reads", true);
    argParser.add<std::string>("output", 'o', "the output directory", true);
    argParser.add<int>("threads", 't', "number of threads", false, 50);

    argParser.parse_check(argc, argv);

    DBGRunner dbgrunner = DBGRunner(argParser.get<std::string>("reads"), argParser.get<std::string>("output"), argParser.get<int>("threads"));
    // MDBGRunner(argParser.get<std::string>("output"), argParser.get<std::string>("reads"), argParser.get<int>("threads"));
    return 0;
}
