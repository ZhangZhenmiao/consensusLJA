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
#include <sstream>

namespace {

std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (char c : value) {
        if (c == '\'')
            quoted += "'\\''";
        else
            quoted += c;
    }
    quoted += "'";
    return quoted;
}

} // namespace

int execute_command(const std::string& command, bool exit_when_fail,
                    bool mute_output, const std::string& log_file) {
    std::string actual_command = command;
    if (!log_file.empty()) {
        std::ofstream log(log_file, std::ios::app);
        if (!log)
            throw std::runtime_error("Failed to open command log: " + log_file);
        log << "\n[CMD] " << command << '\n';
        log.close();
        actual_command += " >> " + shell_quote(log_file) + " 2>&1";
        std::cout << "[CMD] Log: " << log_file << std::endl;
    }
    else if (mute_output) {
        if (actual_command.find(">") == std::string::npos)
            actual_command += " > /dev/null 2>&1";
        else
            actual_command += " 2>/dev/null";
    }
    // std::cout << "[CMD] Executing command: " << command << (mute_output ? " (output muted)" : ", logs are below:") << std::endl;
    struct CommandResult {
        bool success;
        int exit_code;
    };
    CommandResult result;
    result.success = false;
    result.exit_code = -1;

    const int wait_status = system(actual_command.c_str());

    // Check if command executed successfully
    if (wait_status == -1) {
        result.exit_code = -1;
    }
    else if (WIFEXITED(wait_status)) {
        result.exit_code = WEXITSTATUS(wait_status);
        result.success = (result.exit_code == 0);
    }
    else if (WIFSIGNALED(wait_status)) {
        result.exit_code = 128 + WTERMSIG(wait_status);
    }

    std::string status;

    // status += "Command " + command + " " + std::string(result.success ? "succeeded" : "failed");
    status += " with exit code " + std::to_string(result.exit_code) + "\n";
    // std::cout << status << std::endl;

    if (!result.success && exit_when_fail) {
        std::ostringstream message;
        message << "Failed to execute command (exit code " << result.exit_code << "): "
                << command;
        if (!log_file.empty())
            message << "\nSee log: " << log_file;
        throw std::runtime_error(message.str());
    }

    // std::cout << "[CMD] " << command << " " << (result.success ? "succeeded" : "failed") << std::endl;

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

std::string replaceNsWithRandomBases(const std::string& seq) {
    static const char bases[] = { 'A', 'C', 'G', 'T' };
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, 3);

    std::string result = seq;
    for (char& c : result) {
        if (c == 'N' || c == 'n') {
            c = bases[dis(gen)];
        }
    }
    return result;
}

std::string format_with_commas(int value) {
    std::string num = std::to_string(value);
    int insertPosition = num.length() - 3;

    while (insertPosition > 0) {
        num.insert(insertPosition, ",");
        insertPosition -= 3;
    }

    return num;
}

std::pair<double, double> calculate_identities_from_cigar(const std::string& cigar, int gap_threshold) {
    int matches = 0;
    int mismatches = 0;
    int insertions = 0;
    int deletions = 0;
    int long_gaps = 0;

    int query_len = 0;
    int ref_len = 0;
    int aligned_bases = 0;

    std::regex cigar_regex(R"((\d+)([MIDNSHP=X]))");
    auto words_begin = std::sregex_iterator(cigar.begin(), cigar.end(), cigar_regex);
    auto words_end = std::sregex_iterator();

    for (auto it = words_begin; it != words_end; ++it) {
        int length = std::stoi((*it)[1]);
        char op = (*it)[2].str()[0];

        if (op == 'M' || op == '=' || op == 'X') {
            query_len += length;
            ref_len += length;
            aligned_bases += length;

            if (op == '=') {
                matches += length;
            }
            else if (op == 'X') {
                mismatches += length;
            }
            else if (op == 'M') {
                matches += length;  // Assumes M = match
            }
        }
        else if (op == 'I') {
            query_len += length;
            insertions += length;
            aligned_bases += length;
            if (length >= gap_threshold) {
                long_gaps += length;
            }
        }
        else if (op == 'D') {
            ref_len += length;
            deletions += length;
            aligned_bases += length;
            if (length >= gap_threshold) {
                long_gaps += length;
            }
        }
        // S, H, N, P are ignored
    }

    double denom = matches + mismatches + insertions + deletions;
    double denom_no_gap = denom - long_gaps;

    double identity = (denom > 0) ? matches / denom : 0.0;
    double identity_nogap = (denom_no_gap > 0) ? matches / denom_no_gap : 0.0;

    return std::make_pair(identity, identity_nogap);
}
