#include <regex>
#include <random>
#include "utils.hpp"
#include <cstdlib>
#include <fstream>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <iostream>
#include <limits.h>

int execute_command(const std::string& command) {
    struct CommandResult {
        bool success;
        int exit_code;
    };
    CommandResult result;
    result.success = false;
    result.exit_code = -1;

    result.exit_code = system(command.c_str());

    // Check if command executed successfully
    if (WIFEXITED(result.exit_code)) {
        result.exit_code = WEXITSTATUS(result.exit_code);
        result.success = (result.exit_code == 0);
    }

    std::string status;

    status += "Command " + command + " " + std::string(result.success ? "succeeded" : "failed");
    status += " with exit code " + std::to_string(result.exit_code) + "\n";
    // std::cout << status << std::endl;

    if (!result.success)
        std::cout << "Failed to execute " + command << std::endl;

    return result.exit_code;
}

std::string getExecutablePath() {
    char buffer[1024];
    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len != -1) {
        buffer[len] = '\0';
        std::string execPath = std::string(buffer);
        return execPath.substr(0, execPath.find_last_of("/"));
    }
    return "";
}


std::string replace_N(const std::string& seq) {
    static const char nucleotides[] = { 'A', 'C', 'G', 'T' };
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, 3);
    std::string out = seq;
    for (char& c : out) {
        if (c == 'N' || c == 'n') c = nucleotides[dis(gen)];
    }
    return out;
}

double parse_cigar_identity(const std::string& cigar, size_t len1, size_t len2) {
    std::regex re("(\\d+)([MIDNSHP=X])");
    auto begin = std::sregex_iterator(cigar.begin(), cigar.end(), re);
    auto end = std::sregex_iterator();
    size_t matches = 0;
    for (auto i = begin; i != end; ++i) {
        int length = std::stoi((*i)[1]);
        char op = (*i)[2].str()[0];
        if (op == 'M' || op == '=') matches += length;
    }
    size_t shorter_len = std::min(len1, len2);
    return shorter_len > 0 ? (double)matches / shorter_len : 0.0;
}

double unialigner_identity(const std::string& seq1, const std::string& seq2, int prefix) {
    // Preprocess: replace N, keep first 1Mbp, reverse complement second
    std::string s1 = replace_N(seq1.substr(0, prefix));
    std::string s2 = replace_N(seq2.substr(0, prefix));

    // Write temp FASTA files
    std::ofstream f1("seq1_tmp.fasta");
    f1 << ">seq1\n" << s1 << "\n";
    f1.close();
    std::ofstream f2("seq2_tmp.fasta");
    f2 << ">seq2\n" << s2 << "\n";
    f2.close();

    // Run unialigner
    execute_command("mkdir -p unialigner_out");
    std::string cmd = "/Poppy/zmzhang/software/unialigner_new/tandem_aligner/build/bin/tandem_aligner --first seq1_tmp.fasta --second seq2_tmp.fasta -o unialigner_out > /dev/null 2>&1";
    int ret = execute_command(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Unialigner failed\n";
        return -1.0;
    }

    // Read CIGAR
    std::ifstream cigar_file("unialigner_out/cigar.txt");
    if (!cigar_file) {
        std::cerr << "CIGAR file not found\n";
        return -1.0;
    }
    std::string cigar;
    std::getline(cigar_file, cigar);
    cigar_file.close();

    std::remove("seq1_tmp.fasta");
    std::remove("seq2_tmp.fasta");
    execute_command("rm -rf unialigner_out");

    return parse_cigar_identity(cigar, s1.size(), s2.size());
}