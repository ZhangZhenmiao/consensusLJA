#include "dot_graph.hpp"
#include <sstream>
#include <iomanip>
#include <filesystem>

std::string Graph::doubleToString(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << value;
    return stream.str();
}

// Get reverse complementary edge_sequence
std::string Graph::reverse_complementary(std::string& seq) {
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

std::string Graph::reverse_complementary_node(std::string node) {
    return nodeid2Rev.at(node);
}

int Graph::count_matches(std::string cigar) {
    int matches = 0;
    int num = 0;
    for (char c : cigar) {
        if (std::isdigit(c))
            num = num * 10 + (c - '0');
        else {
            if (c == 'M')
                matches += num;
            num = 0;
        }
    }
    return matches;
}

void Graph::get_annotation(std::string prefix) {
    std::string ref_seq = "/Poppy/zmzhang/Rust_fungi/genome/reference.compressed.only_chrs.fasta";
    if (system(("minimap2 -ax asm20 " + ref_seq + " " + prefix + ".fasta -t 100 | grep -v '^@' > " + prefix + ".ref.sam").c_str()) != 0) {
        exit(1);
    }
    if (!std::filesystem::exists(ref_seq + ".fai")) {
        if (system(("samtools faidx " + ref_seq).c_str()) != 0)
            exit(1);
    }
    if (system(("cut -f1,2 " + ref_seq + ".fai | awk " + R"('{print "@SQ\tSN:"$1"\tLN:"$2}')" + " > " + prefix + ".ref.header.sam").c_str()) != 0)
        exit(1);
    if (system(("cat " + prefix + ".ref.header.sam " + prefix + ".ref.sam | samtools sort -@ 50 -o " + prefix + ".ref.bam").c_str()) != 0) {
        exit(1);
    }
    std::string exeDir = getExecutablePath();
    if (system((exeDir + "/../src/scripts/get_reference.py -o " + prefix + ".ref.bam.stats " + prefix + ".ref.bam " + prefix + ".fasta").c_str()) != 0)
        exit(1);
    write_graph_colored_from_bam(prefix + ".color", prefix + ".ref.bam.stats");
}