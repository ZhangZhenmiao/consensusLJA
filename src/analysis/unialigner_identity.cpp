#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <regex>
#include <string>
#include <unordered_map>

#include "edlib.h"

struct Identities {
    double pi1;
    double pi2;
};

Identities calculate_identities_from_cigar(const std::string& cigar,
                                            int gap_threshold = 10) {
    int matches = 0;
    int mismatches = 0;
    int insertions = 0;
    int deletions = 0;
    int long_gaps = 0;

    const std::regex cigar_regex(R"((\d+)([MIDNSHP=X]))");
    for (auto it = std::sregex_iterator(cigar.begin(), cigar.end(), cigar_regex);
         it != std::sregex_iterator(); ++it) {
        const int length = std::stoi((*it)[1]);
        const char op = (*it)[2].str()[0];

        if (op == '=' || op == 'M') {
            // Some aligners emit M without distinguishing matches from mismatches.
            matches += length;
        } else if (op == 'X') {
            mismatches += length;
        } else if (op == 'I') {
            insertions += length;
            if (length >= gap_threshold) long_gaps += length;
        } else if (op == 'D') {
            deletions += length;
            if (length >= gap_threshold) long_gaps += length;
        }
        // S, H, N and P do not contribute to these identity definitions.
    }

    const double denominator = matches + mismatches + insertions + deletions;
    const double denominator_without_long_gaps = denominator - long_gaps;
    return {
        denominator > 0 ? matches / denominator : 0.0,
        denominator_without_long_gaps > 0
            ? matches / denominator_without_long_gaps
            : 0.0
    };
}

int count_matches(const std::string& cigar) {
    int matches = 0;
    int length = 0;
    for (const char c : cigar) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            length = length * 10 + (c - '0');
        } else {
            if (c == '=' || c == 'M') matches += length;
            length = 0;
        }
    }
    return matches;
}

void report_identities(const std::string& cigar, size_t seq1_len, size_t seq2_len) {
    const Identities identities = calculate_identities_from_cigar(cigar);
    const size_t shorter_len = std::min(seq1_len, seq2_len);
    const double pi3 = shorter_len > 0
        ? static_cast<double>(count_matches(cigar)) / shorter_len
        : 0.0;

    std::cout << "PI1: " << identities.pi1 << '\n'
              << "PI2: " << identities.pi2 << '\n'
              << "PI3: " << pi3 << '\n';
}

bool run_unialigner(const std::string& seq1, const std::string& seq2) {
    if (system("mkdir -p unialigner_out") != 0) {
        std::cerr << "Failed to create unialigner_out\n";
        return false;
    }

    std::ofstream f1("unialigner_out/seq1_tmp.fasta");
    std::ofstream f2("unialigner_out/seq2_tmp.fasta");
    if (!f1 || !f2) {
        std::cerr << "Failed to create temporary FASTA files\n";
        return false;
    }
    f1 << ">seq1\n" << seq1 << '\n';
    f2 << ">seq2\n" << seq2 << '\n';
    f1.close();
    f2.close();

    const std::string command =
        "/Poppy/zmzhang/software/unialigner_new/tandem_aligner/build/bin/"
        "tandem_aligner --first unialigner_out/seq1_tmp.fasta --second "
        "unialigner_out/seq2_tmp.fasta -o unialigner_out > /dev/null 2>&1";
    if (system(command.c_str()) != 0) {
        std::cerr << "Unialigner failed\n";
        return false;
    }

    std::ifstream cigar_file("unialigner_out/cigar.txt");
    std::string cigar;
    if (!cigar_file || !std::getline(cigar_file, cigar) || cigar.empty()) {
        std::cerr << "Cannot read a CIGAR from unialigner_out/cigar.txt\n";
        return false;
    }
    report_identities(cigar, seq1.size(), seq2.size());
    return true;
}

bool run_edlib(const std::string& seq1, const std::string& seq2) {
    if (seq1.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        seq2.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        std::cerr << "Sequences are too long for edlib\n";
        return false;
    }

    EdlibAlignResult result = edlibAlign(
        seq1.c_str(), static_cast<int>(seq1.size()),
        seq2.c_str(), static_cast<int>(seq2.size()),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));
    if (result.status != EDLIB_STATUS_OK || result.editDistance < 0 ||
        result.alignment == nullptr) {
        std::cerr << "Edlib alignment failed\n";
        edlibFreeAlignResult(result);
        return false;
    }

    char* cigar_ptr = edlibAlignmentToCigar(
        result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
    if (cigar_ptr == nullptr) {
        std::cerr << "Edlib failed to generate a CIGAR\n";
        edlibFreeAlignResult(result);
        return false;
    }
    const std::string cigar(cigar_ptr);
    std::free(cigar_ptr);
    edlibFreeAlignResult(result);

    report_identities(cigar, seq1.size(), seq2.size());
    return true;
}

std::unordered_map<std::string, std::string> read_fasta_map(
        const std::string& filename) {
    std::unordered_map<std::string, std::string> contigs;
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "Cannot open file: " << filename << '\n';
        return contigs;
    }

    std::string line;
    std::string name;
    std::string sequence;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) contigs[name] = sequence;
            name = line.substr(1);
            sequence.clear();
        } else {
            sequence += line;
        }
    }
    if (!name.empty()) contigs[name] = sequence;
    return contigs;
}

int main(int argc, char* argv[]) {
    if (argc != 3 && argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <fasta1> <fasta2> [--aligner unialigner|edlib]\n"
                  << "  Each FASTA should contain a single sequence.\n";
        return 1;
    }

    std::string aligner = "unialigner";
    if (argc == 5) {
        if (std::string(argv[3]) != "--aligner") {
            std::cerr << "Expected --aligner before the aligner name\n";
            return 1;
        }
        aligner = argv[4];
        if (aligner != "unialigner" && aligner != "edlib") {
            std::cerr << "Unknown aligner: " << aligner << '\n';
            return 1;
        }
    }

    const auto contigs1 = read_fasta_map(argv[1]);
    const auto contigs2 = read_fasta_map(argv[2]);
    if (contigs1.empty() || contigs2.empty()) {
        std::cerr << "Could not find a sequence in one or both FASTA files.\n";
        return 1;
    }
    if (contigs1.size() != 1 || contigs2.size() != 1) {
        std::cerr << "Each FASTA must contain exactly one sequence.\n";
        return 1;
    }

    const std::string& seq1 = contigs1.begin()->second;
    const std::string& seq2 = contigs2.begin()->second;
    const bool success = aligner == "edlib"
        ? run_edlib(seq1, seq2)
        : run_unialigner(seq1, seq2);
    return success ? 0 : 1;
}
