#include "edlib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
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

double matches_by_edlib(std::string sequence1, std::string sequence2, bool lcs = false) {
    // std::string seq_short, seq_long;
    // if (sequence1.size() <= sequence2.size())
    //     sequence2 = sequence2.substr(0, sequence1.size());
    // else
    //     sequence1 = sequence1.substr(0, sequence2.size());
    EdlibAlignResult result = edlibAlign(sequence1.c_str(), sequence1.size(), sequence2.c_str(), sequence2.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        std::string cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
        edlibFreeAlignResult(result);
        int lcs_len = count_matches(cigar);
        if (lcs)
            return lcs_len;
        else
            return 1.0 * lcs_len / std::min(sequence1.size(), sequence1.size());
    }
    else {
        std::cout << "edlib failed" << std::endl;
        return 0;
    }
}


double unialigner_identity(const std::string& seq1, const std::string& seq2) {
    // Preprocess: replace N, keep first 1Mbp, reverse complement second
    std::string s1 = replace_N(seq1);
    std::string s2 = replace_N(seq2);
    // if (seq1.size() <= seq2.size())
    //     s2 = s2.substr(0, seq1.size());
    // else
    //     s1 = s1.substr(0, seq2.size());

    // Write temp FASTA files
    std::ofstream f1("seq1_tmp.fasta");
    f1 << ">seq1\n" << s1 << "\n";
    f1.close();
    std::ofstream f2("seq2_tmp.fasta");
    f2 << ">seq2\n" << s2 << "\n";
    f2.close();

    // Run unialigner
    system("mkdir -p unialigner_out");
    std::string cmd = "/Poppy/zmzhang/software/unialigner_new/tandem_aligner/build/bin/tandem_aligner --first seq1_tmp.fasta --second seq2_tmp.fasta -o unialigner_out > /dev/null 2>&1";
    int ret = system(cmd.c_str());
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

    // std::remove("seq1_tmp.fasta");
    // std::remove("seq2_tmp.fasta");
    // system("rm -rf unialigner_out");

    return parse_cigar_identity(cigar, s1.size(), s2.size());
}

// Read all contigs from a fasta into a map
#include <unordered_map>
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

    // Take 1Mbp prefix if desired
    // seq1 = seq1.substr(0, 1000000);
    // seq2 = seq2.substr(0, 1000000);

    double idt = unialigner_identity(seq1, seq2);
    std::cout << "Unialigner identity: " << idt << std::endl;
    idt = matches_by_edlib(seq1, seq2);
    std::cout << "Edlib identity: " << idt << std::endl;
    return 0;
}
