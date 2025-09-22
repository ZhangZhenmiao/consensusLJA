#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <zlib.h>
#include <stdexcept>

// --- Dynamic line reading from gzFile ---
std::string gz_readline(gzFile file) {
    if (!file) throw std::runtime_error("Invalid gzFile");
    std::string line;
    char c;
    while (true) {
        int ret = gzread(file, &c, 1);
        if (ret < 0) throw std::runtime_error("gzread error");
        if (ret == 0 || c == '\n') break; // EOF or newline
        line += c;
    }
    return line;
}

void gz_writeline(gzFile file, const std::string &line) {
    if (!file) throw std::runtime_error("Invalid gzFile");
    std::string out = line + "\n";
    if (gzwrite(file, out.data(), out.size()) != (int)out.size()) {
        throw std::runtime_error("gzwrite failed");
    }
}

// --- Replace non-ATCG bases ---
bool replace_non_atcg(std::string &seq, std::mt19937 &rng) {
    bool changed = false;
    std::uniform_int_distribution<int> dist(0, 3);
    const char bases[] = {'A', 'T', 'C', 'G'};
    for (char &c : seq) {
        if (c != 'A' && c != 'T' && c != 'C' && c != 'G' &&
            c != 'a' && c != 't' && c != 'c' && c != 'g') {
            c = bases[dist(rng)];
            changed = true;
        }
    }
    return changed;
}

// --- Open gzipped or plain input ---
std::istream* open_input_file(const std::string &filename, std::ifstream &file, gzFile &gzfile) {
    if (filename.size() >= 3 && filename.substr(filename.size()-3) == ".gz") {
        gzfile = gzopen(filename.c_str(), "rb");
        if (!gzfile) throw std::runtime_error("Failed to open gzip file: " + filename);
        return nullptr; // use gzfile
    } else {
        file.open(filename);
        if (!file.is_open()) throw std::runtime_error("Failed to open file: " + filename);
        return &file;
    }
}

// --- Open gzipped or plain output ---
std::ostream* open_output_file(const std::string &filename, std::ofstream &file, gzFile &gzfile) {
    if (filename.size() >= 3 && filename.substr(filename.size()-3) == ".gz") {
        gzfile = gzopen(filename.c_str(), "wb");
        if (!gzfile) throw std::runtime_error("Failed to open gzip file: " + filename);
        return nullptr; // use gzfile
    } else {
        file.open(filename, std::ios::out);
        if (!file.is_open()) throw std::runtime_error("Failed to open file: " + filename);
        return &file;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.fastq[.gz] output.fastq[.gz] [seed]\n";
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];
    int seed = (argc >= 4) ? std::stoi(argv[3]) : 2025;

    std::mt19937 rng(seed);

    std::ifstream fin_file;
    std::ofstream fout_file;
    gzFile fin_gz = nullptr, fout_gz = nullptr;

    std::istream *fin = open_input_file(input_file, fin_file, fin_gz);
    std::ostream *fout = open_output_file(output_file, fout_file, fout_gz);

    size_t total_reads = 0;
    size_t changed_reads = 0;

    while (true) {
        std::string id, seq, plus, qual;
        try {
            if (fin) {
                if (!std::getline(*fin, id)) break;
                if (!std::getline(*fin, seq)) break;
                if (!std::getline(*fin, plus)) break;
                if (!std::getline(*fin, qual)) break;
            } else {
                id   = gz_readline(fin_gz);
                seq  = gz_readline(fin_gz);
                plus = gz_readline(fin_gz);
                qual = gz_readline(fin_gz);
                if (id.empty()) break;
            }
        } catch (...) { break; }

        total_reads++;
        if (replace_non_atcg(seq, rng)) changed_reads++;

        // Write
        if (fout) {
            *fout << id << "\n" << seq << "\n" << plus << "\n" << qual << "\n";
        } else {
            gz_writeline(fout_gz, id);
            gz_writeline(fout_gz, seq);
            gz_writeline(fout_gz, plus);
            gz_writeline(fout_gz, qual);
        }
    }

    if (fin_gz) gzclose(fin_gz);
    if (fout_gz) gzclose(fout_gz);
    if (fin_file.is_open()) fin_file.close();
    if (fout_file.is_open()) fout_file.close();

    std::cout << "Total reads processed: " << total_reads << "\n";
    std::cout << "Reads changed (non-ATCG replaced): " << changed_reads << "\n";
    std::cout << "Fraction changed: " << (double)changed_reads / total_reads * 100 << "%\n";

    return 0;
}
