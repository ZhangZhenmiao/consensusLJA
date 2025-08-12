#include "edlib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <unordered_map>
#include <random>

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

std::pair<double, double> calculate_identities_from_cigar(const std::string& cigar, int gap_threshold = 10) {
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

std::string reverse_complement(std::string seq) {
    std::string out;
    for (int i = seq.size() - 1; i >= 0; --i) {
        if (seq.at(i) == 'A')
            out += 'T';
        else if (seq.at(i) == 'T')
            out += 'A';
        else if (seq.at(i) == 'C')
            out += 'G';
        else if (seq.at(i) == 'G')
            out += 'C';
        else
            out += seq.at(i);
    }
    return std::move(out);
}

int count_matches(std::string cigar) {
    int matches = 0;
    int num = 0;
    for (char c : cigar) {
        if (std::isdigit(c))
            num = num * 10 + (c - '0');
        else {
            if (c == '=')
                matches += num;
            num = 0;
        }
    }
    return matches;
}

void edlib_identity(std::string sequence1, std::string sequence2) {
    EdlibAlignResult result = edlibAlign(sequence1.c_str(), sequence1.size(), sequence2.c_str(), sequence2.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    std::string cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
    edlibFreeAlignResult(result);
    auto idts = calculate_identities_from_cigar(cigar);
    std::cout << "Edlib identity: " << idts.first << " " << idts.second << std::endl;
    std::cout << 1.0 * count_matches(cigar) / std::min(sequence1.size(), sequence2.size()) << std::endl;
}


void unialigner_identity(const std::string& seq1, const std::string& seq2) {
    std::string s1 = replace_N(seq1);
    std::string s2 = replace_N(seq2);

    system("mkdir -p unialigner_out");

    // Write temp FASTA files
    std::ofstream f1("unialigner_out/seq1_tmp.fasta");
    f1 << ">seq1\n" << s1 << "\n";
    f1.close();
    std::ofstream f2("unialigner_out/seq2_tmp.fasta");
    f2 << ">seq2\n" << s2 << "\n";
    f2.close();

    // Run unialigner
    std::string cmd = "/Poppy/zmzhang/software/unialigner_new/tandem_aligner/build/bin/tandem_aligner --first unialigner_out/seq1_tmp.fasta --second unialigner_out/seq2_tmp.fasta -o unialigner_out > /dev/null 2>&1";
    system(cmd.c_str());

    // Read CIGAR
    std::ifstream cigar_file("unialigner_out/cigar.txt");
    std::string cigar;
    std::getline(cigar_file, cigar);
    cigar_file.close();

    auto idts = calculate_identities_from_cigar(cigar);
    std::cout << "Unialigner identity: " << idts.first << " " << idts.second << std::endl;
}

// Read all contigs from a fasta into a map
std::unordered_map<std::string, std::string> read_fasta_map(const std::string& filename) {
    std::unordered_map<std::string, std::string> contigs;
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        exit(1);
    }
    std::string line, name, seq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) contigs[name] = seq;
            name = line.substr(1);
            seq.clear();
        }
        else {
            seq += line;
        }
    }
    if (!name.empty()) contigs[name] = seq;
    return contigs;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fasta1> <fasta2>" << std::endl;
        std::cerr << "  Each fasta should contain a single sequence." << std::endl;
        return 1;
    }
    std::string fasta1 = argv[1];
    std::string fasta2 = argv[2];
    auto contigs1 = read_fasta_map(fasta1);
    auto contigs2 = read_fasta_map(fasta2);
    if (contigs1.empty() || contigs2.empty()) {
        std::cerr << "Could not find sequence in one or both fasta files." << std::endl;
        return 1;
    }
    // Use the first sequence in each fasta
    std::string seq1 = contigs1.begin()->second;
    std::string seq2 = contigs2.begin()->second;

    // // Take 1Mbp prefix if desired
    // seq1 = seq1.substr(0, 100000);
    // seq2 = seq2.substr(0, 100000);

    unialigner_identity(seq1, seq2);
    edlib_identity(seq1, seq2);
    return 0;
}
