#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>
#include <htslib/sam.h>
#include <cassert>
#include <algorithm>
#include <numeric>

// Calculate percent identity from a BAM alignment
double calculate_identity(const bam1_t* aln) {
    uint32_t* cigar = bam_get_cigar(aln);
    int n_cigar = aln->core.n_cigar;
    int matches = 0, mismatches = 0, insertions = 0, deletions = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH) matches += len;
        else if (op == BAM_CEQUAL) matches += len;
        else if (op == BAM_CDIFF) mismatches += len;
        else if (op == BAM_CINS) insertions += len;
        else if (op == BAM_CDEL) deletions += len;
    }
    int aligned_len = matches + mismatches;
    int total_bases = aligned_len + insertions + deletions;
    return total_bases > 0 ? (double)matches / total_bases : 0.0;
}

// Parse node lengths from DOT file
std::unordered_map<std::string, int> parse_dot_file(const std::string& dot_file_path) {
    std::unordered_map<std::string, int> node2len;
    std::ifstream dot_file(dot_file_path);
    std::string line;
    while (std::getline(dot_file, line)) {
        if (line.find("label") == std::string::npos) continue;
        if (line.find("->") != std::string::npos) continue;
        std::string node = line.substr(0, line.find(' '));
        if (node.find('+') == std::string::npos) {
            size_t lpos = line.find('L');
            size_t qpos = line.find('"', lpos + 1);
            std::string node_len = line.substr(lpos + 1, qpos - lpos - 1);
            node2len[node] = std::stoi(node_len);
        }
    }
    return node2len;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <bam_file> <dot_file> <output_file>\n";
        return 1;
    }
    std::string bam_file_path = argv[1];
    std::string dot_file_path = argv[2];
    std::string output_file = argv[3];

    auto node2len = parse_dot_file(dot_file_path);

    samFile* fp = sam_open(bam_file_path.c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    bam1_t* aln = bam_init1();

    std::unordered_map<std::string, std::vector<std::tuple<std::string, int, int, int, int, double>>> alignments;
    std::set<std::string> chimeric_edges;

    while (sam_read1(fp, hdr, aln) >= 0) {
        if (aln->core.flag & BAM_FUNMAP) continue;
        std::string ref_name = hdr->target_name[aln->core.tid];
        double idt = calculate_identity(aln);
        if (idt > 0.99) {
            alignments[ref_name].emplace_back(
                bam_get_qname(aln),
                aln->core.pos,
                bam_endpos(aln),
                hdr->target_len[aln->core.tid],
                aln->core.l_qseq,
                idt
            );
        }
    }

    // Chimeric detection logic (simplified, adapt as needed)
    int tolerant_size = 1000;
    for (const auto& kv : alignments) {
        const std::string& r = kv.first;
        const auto& alns = kv.second;
        std::vector<std::string> nodes;
        size_t uscore = r.find('_');
        if (uscore == std::string::npos) continue;
        nodes.push_back(r.substr(0, uscore));
        nodes.push_back(r.substr(uscore + 1));
        std::string node1 = nodes[0].substr(0, nodes[0].find('.'));
        std::string node2 = nodes[1].substr(0, nodes[1].find('.'));
        int contig_len = hdr->target_len[sam_hdr_name2tid(hdr, r.c_str())];
        int split_coordinate1 = std::min(node2len[node1], contig_len - node2len[node2]);
        int split_coordinate2 = std::max(node2len[node1], contig_len - node2len[node2]);
        if (split_coordinate1 <= 2 * tolerant_size || split_coordinate2 >= contig_len - 2 * tolerant_size) continue;
        if (node2len[node1] + node2len[node2] - contig_len >= 10000 - tolerant_size) continue;

        int cnt_f = 0, cnt_r = 0, cnt_all = 0;
        for (const auto& aln_tuple : alns) {
            int start = std::get<1>(aln_tuple);
            int end = std::get<2>(aln_tuple);
            if (start <= std::max(split_coordinate1 - tolerant_size, 0) && end > node2len[node1]) cnt_f++;
            if (start < contig_len - node2len[node2] && end >= std::min(split_coordinate2 + tolerant_size, contig_len)) cnt_r++;
            cnt_all++;
        }
        if (cnt_f <= 0 || cnt_r <= 0) {
            std::cout << "contig name " << r << " contig len " << contig_len << " len 1 " << node2len[node1] << " len 2 " << node2len[node2]
                << " supporting reads " << cnt_f << " " << cnt_r << " " << cnt_all << std::endl;
            std::cout << r << " is chimeric" << std::endl;
            chimeric_edges.insert(node1);
            chimeric_edges.insert(node2);
        }
    }

    // Internal chimeric detection
    tolerant_size = 100;
    for (const auto& kv : alignments) {
        const std::string& r = kv.first;
        int ref_len = hdr->target_len[sam_hdr_name2tid(hdr, r.c_str())];
        std::vector<int> coverages(ref_len, 0);

        // Build coverage array
        for (const auto& aln_tuple : kv.second) {
            int start = std::get<1>(aln_tuple);
            int end = std::get<2>(aln_tuple);
            if (end - tolerant_size > start) {
                int cov_start = std::max(start, 0);
                int cov_end = std::max(end - tolerant_size, cov_start);
                for (int i = cov_start; i < cov_end && i < ref_len; ++i) {
                    coverages[i]++;
                }
            }
        }

        // Find potential positions with low coverage
        std::vector<int> potential_pos;
        double avg_cov = std::accumulate(coverages.begin(), coverages.end(), 0.0) / ref_len;
        if (avg_cov >= 5) {
            for (int i = 0; i < ref_len; ++i) {
                if (coverages[i] <= 1 &&
                    i >= 5001 - tolerant_size &&
                    i <= ref_len - 5001 + tolerant_size) {
                    potential_pos.push_back(i);
                }
            }
        }

        // Find reads supporting chimeric positions
        std::set<std::string> reads;
        for (const auto& aln_tuple : kv.second) {
            int aln_start = std::get<1>(aln_tuple);
            int aln_end = std::get<2>(aln_tuple);
            std::string aln_name = std::get<0>(aln_tuple);
            for (int j : potential_pos) {
                if (aln_start <= j - 5001 + tolerant_size &&
                    aln_end >= j + 5001 - tolerant_size) {
                    reads.insert(aln_name);
                }
            }
        }

        // Print IDs of chimeric reads
        if (!reads.empty()) {
            std::cout << r << " is chimeric (internal)" << std::endl;
            size_t uscore = r.find('_');
            if (uscore != std::string::npos) {
                chimeric_edges.insert(r.substr(0, uscore));
                chimeric_edges.insert(r.substr(uscore + 1));
            }
        }
    }

    // Output chimeric edges
    std::ofstream fout(output_file);
    for (const auto& e : chimeric_edges) {
        fout << e << std::endl;
    }

    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}