#include <bits/stdc++.h>
#include <regex>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

static const size_t SEG_LEN = 2000000; // 2 Mb

string rc(const string& s) {
    string out; out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        char c = *it;
        switch (c) {
            case 'A': out.push_back('T'); break;
            case 'C': out.push_back('G'); break;
            case 'G': out.push_back('C'); break;
            case 'T': out.push_back('A'); break;
            case 'a': out.push_back('t'); break;
            case 'c': out.push_back('g'); break;
            case 'g': out.push_back('c'); break;
            case 't': out.push_back('a'); break;
            default:  out.push_back('N'); break;
        }
    }
    return out;
}

// Very simple FASTA parser that yields pairs (header, sequence)
vector<pair<string,string>> read_fasta(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open fasta: " + path);
    vector<pair<string,string>> out;
    string line;
    string header;
    string seq;
    while (std::getline(in, line)) {
        if (line.size() && line[0] == '>') {
            if (!header.empty()) {
                out.emplace_back(header, seq);
            }
            header = line.substr(1);
            // trim header at first whitespace to get contig name only (common choice)
            auto pos = header.find_first_of(" \t");
            if (pos != string::npos) header = header.substr(0, pos);
            seq.clear();
        } else {
            // append sequence; remove whitespace
            for (char c : line) if (!isspace((unsigned char)c)) seq.push_back(c);
        }
    }
    if (!header.empty()) out.emplace_back(header, seq);
    return out;
}

// replace_N as provided
string replace_N(const string& seq) {
    static const char nucleotides[] = { 'A', 'C', 'G', 'T' };
    static std::mt19937 gen(2025);
    static std::uniform_int_distribution<> dis(0, 3);
    string out = seq;
    for (char& c : out) {
        if (c == 'N' || c == 'n') c = nucleotides[dis(gen)];
    }
    return out;
}

pair<double,double> calculate_identities_from_cigar(const string& cigar, int gap_threshold = 10) {
    int matches=0, mismatches=0, insertions=0, deletions=0, long_gaps=0;
    int query_len=0, ref_len=0, aligned_bases=0;
    regex cigar_regex(R"((\d+)([MIDNSHP=X]))");
    auto words_begin = sregex_iterator(cigar.begin(), cigar.end(), cigar_regex);
    auto words_end = sregex_iterator();

    for (auto it = words_begin; it != words_end; ++it) {
        int length = stoi((*it)[1]);
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

    double identity = (denom > 0) ? (double)matches / denom : 0.0;
    double identity_nogap = (denom_no_gap > 0) ? (double)matches / denom_no_gap : 0.0;

    return make_pair(identity, identity_nogap);
}

string run_unialigner_and_get_cigar(const string& s1, const string& s2, string out_dir) {

    std::string f1 = out_dir + "/" + s1 + ".fasta";
    std::string f2 = out_dir + "/" + s2 + ".fasta";
    std::string UNIALIGNER_CMD = "/Poppy/zmzhang/software/unialigner_new/tandem_aligner/build/bin/tandem_aligner --first " + f1 + " --second " + f2 + " -o " + out_dir + "/" + s1 + "_" + s2 + ".unialigner > /dev/null 2>&1";
    // Run command
    int ret = system(UNIALIGNER_CMD.c_str());
    (void)ret; // we don't treat nonzero specially here; we will check for cigar file

    // Read cigar
    string cigar;
    ifstream cigar_file(out_dir + "/" + s1 + "_" + s2 + ".unialigner/cigar.txt");
    if (!cigar_file) {
        // no cigar file found -> return empty string
        return string();
    }
    getline(cigar_file, cigar);
    // Trim whitespace
    while (!cigar.empty() && isspace((unsigned char)cigar.back())) cigar.pop_back();
    return cigar;
}

std::string contig_name(std::string header) {
    auto pos = header.find("_start");
    if (pos != string::npos) header = header.substr(0, pos);
    else {
        pos = header.find("_end");
        if (pos != string::npos) header = header.substr(0, pos);
    }
    return header;
}

// --- helper to read already completed results
unordered_set<string> load_done_pairs(const string& out_tsv) {
    unordered_set<string> done;
    ifstream fin(out_tsv);
    if (!fin) return done; // file not existing
    string line;
    getline(fin, line); // skip header if any
    while (getline(fin, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string s1, s2;
        ss >> s1 >> s2;
        if (!s1.empty() && !s2.empty()) {
            done.insert(s1 + "\t" + s2);
        }
    }
    return done;
}

struct SeqEntry { string name; string seq; };

vector<SeqEntry> load_seqs_from_dir(const string& out_dir) {
    vector<SeqEntry> seqs;
    for (auto& entry : fs::directory_iterator(out_dir)) {
        if (entry.path().extension() == ".fasta") {
            auto parsed = read_fasta(entry.path().string());
            for (auto &p : parsed) {
                seqs.push_back({p.first, p.second});
            }
        }
    }
    return seqs;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " input.fasta output\n";
        return 1;
    }
    string fasta_in = argv[1];
    string out_dir = argv[2];
    fs::create_directories(out_dir);
    string out_tsv = out_dir + "/end_to_end_identities.tsv";

    vector<SeqEntry> seqs;

    if (fasta_in == "RESTART") {
        seqs = load_seqs_from_dir(out_dir);
        if (seqs.empty()) {
            cerr << "Restart mode: no fasta files found in " << out_dir << "\n";
            return 2;
        }
        cerr << "Restart mode: loaded " << seqs.size() << " sequences from " << out_dir << "\n";
    } else {
        // normal path: read original fasta and generate start/end FASTAs
        // (keep your original 4N-sequence building + writing code here)

        // Read fasta
        vector<pair<string,string>> contigs;
        try {
            contigs = read_fasta(fasta_in);
        } catch (exception& e) {
            std::cout << "Error reading fasta: " << e.what() << "\n";
            return 2;
        }

        if (contigs.empty()) {
            std::cout << "No contigs read from fasta.\n";
            return 3;
        }
        seqs.reserve(contigs.size() * 4);

        for (auto &p : contigs) {
            const string &name = p.first;
            const string &seq = p.second;
            size_t L = seq.size();

            // start
            string sstart = replace_N(seq.substr(0, min<size_t>(SEG_LEN, L)));
            string send = (L <= SEG_LEN) ? replace_N(seq) : replace_N(seq.substr(L - SEG_LEN, SEG_LEN));

            string sstart_rc = rc(sstart);
            string send_rc = rc(send);

            // names follow: contig_name_start_direct, contig_name_start_reverse, contig_name_end_direct, contig_name_end_reverse
            seqs.push_back({name + "_start_direct", sstart});
            seqs.push_back({name + "_start_reverse", sstart_rc});
            seqs.push_back({name + "_end_direct",   send});
            seqs.push_back({name + "_end_reverse",   send_rc});

            // Write temp fasta files
            if (fs::exists(out_dir + "/" + name + "_start_direct" + ".fasta")) {
                ofstream f1(out_dir + "/" + name + "_start_direct" + ".fasta");
                f1 << ">" + name + "_start_direct" + "\n" << sstart << "\n";
                f1.close();
            }
            if (fs::exists(out_dir + "/" + name + "_start_reverse" + ".fasta")) {
                ofstream f2(out_dir + "/" + name + "_start_reverse" + ".fasta");
                f2 << ">" + name + "_start_reverse" + "\n" << sstart_rc << "\n";
                f2.close();
            }
            if (fs::exists(out_dir + "/" + name + "_end_direct" + ".fasta")) {
                ofstream f3(out_dir + "/" + name + "_end_direct" + ".fasta");
                f3 << ">" + name + "_end_direct" + "\n" << send << "\n";
                f3.close();
            }
            if (fs::exists(out_dir + "/" + name + "_end_reverse" + ".fasta")) {
                ofstream f4(out_dir + "/" + name + "_end_reverse" + ".fasta");
                f4 << ">" + name + "_end_reverse" + "\n" << send_rc << "\n";
                f4.close();
            }
        }
    }

    // Prepare output file
    // Prepare output file
    unordered_set<string> done_pairs;
    bool append_mode = fs::exists(out_tsv);

    if (append_mode) {
        done_pairs = load_done_pairs(out_tsv);
    }

    ofstream out(out_tsv, append_mode ? ios::app : ios::out);
    if (!out) {
        std::cout << "Cannot open output file: " << out_tsv << "\n";
        return 4;
    }
    if (!append_mode) {
        out << "seq1\tseq2\tidentity1\tidentity2\n";
    }
    out << fixed << setprecision(6);

    size_t M = seqs.size();
    std::cout << "Total sequences (4N): " << M << " -> pairs: " << (M * (M - 1) / 2) << "\n";
    // Loop over pairs
    #pragma omp parallel for num_threads(10)
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = i + 1; j < M; ++j) {
            if (contig_name(seqs[i].name) == contig_name(seqs[j].name)) {
                // skip pairs from the same contig
                continue;
            }
            string key = seqs[i].name + "\t" + seqs[j].name;
            if (done_pairs.count(key)) {
#pragma omp critical
                {
                    std::cout << "Skipping already done pair: " << key << "\n" << std::flush;
                }
                continue;
            }
            string cigar = run_unialigner_and_get_cigar(seqs[i].name, seqs[j].name, out_dir);
            if (cigar.empty()) {
#pragma omp critical
                {
                    // fallback: output zeros
                    out << seqs[i].name << "\t" << seqs[j].name << "\t0\t0\n" << std::flush;
                    std::cout << "Warning: no cigar for pair " << seqs[i].name << " vs " << seqs[j].name << "\n";
                }
            }
            else {
                auto [identity, identity_nogap] = calculate_identities_from_cigar(cigar);
#pragma omp critical
                {
                    out << seqs[i].name << "\t" << seqs[j].name << "\t"
                        << identity << "\t"
                        << identity_nogap << "\n";
                    out << seqs[j].name << "\t" << seqs[i].name << "\t"
                        << identity << "\t"
                        << identity_nogap << "\n" << std::flush;
                }
            }
            system(("rm -rf " + out_dir + "/" + seqs[i].name + "_" + seqs[j].name + ".unialigner").c_str());
        }
    }

    out.close();
    std::cout << "Done. Results written to " << out_tsv << "\n";
    return 0;
}
