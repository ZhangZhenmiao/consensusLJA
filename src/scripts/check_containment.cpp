#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <filesystem>
#include <cctype>

namespace fs = std::filesystem;

// Function to reverse complement a DNA sequence
std::string reverseComplement(const std::string& sequence) {
    static const std::unordered_map<char, char> complement = {
        {'A', 'T'}, {'T', 'A'},
        {'C', 'G'}, {'G', 'C'},
        {'N', 'N'}, {'a', 't'},
        {'t', 'a'}, {'c', 'g'},
        {'g', 'c'}, {'n', 'n'}
    };

    std::string rc;
    rc.reserve(sequence.size());

    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        char base = *it;
        if (complement.find(base) != complement.end()) {
            rc += complement.at(base);
        }
        else {
            rc += 'N';  // Handle ambiguous bases
        }
    }

    return rc;
}

/**
 * @brief Check if contigs in a FASTA file are contained within other contigs using Mash,
 *        considering both original and reverse-complemented sequences
 *
 * @param fastaPath Path to the input FASTA file
 * @param kmerSize K-mer size for Mash (default 21)
 * @param sketchSize Sketch size for Mash (default 1000)
 * @param threshold Similarity threshold for containment (default 0.99)
 * @return std::unordered_map<std::string, std::vector<std::pair<std::string, bool>>>
 *         Map of contig to list of pairs (contained contig, isReverseComplement)
 */
std::unordered_map<std::string, std::vector<std::pair<std::string, bool>>>
checkContigContainmentWithRC(const std::string& fastaPath,
    int kmerSize = 21,
    int sketchSize = 1000,
    double threshold = 0.8) {

    std::unordered_map<std::string, std::vector<std::pair<std::string, bool>>> containmentMap;

    // Step 1: Create temporary directory for Mash files
    fs::path tempDir = fs::temp_directory_path() / "mash_containment_check_rc";
    if (!fs::exists(tempDir)) {
        fs::create_directory(tempDir);
    }

    // Step 2: Split FASTA into individual contig files (original and reverse complement)
    std::ifstream fastaFile(fastaPath);
    if (!fastaFile.is_open()) {
        throw std::runtime_error("Could not open FASTA file: " + fastaPath);
    }

    std::vector<std::pair<std::string, bool>> contigFiles;  // (path, isReverseComplement)
    std::string line, currentSequence, currentContig;
    std::ofstream currentOut, currentRcOut;

    while (getline(fastaFile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Process previous contig if any
            if (!currentContig.empty()) {
                // Write original sequence
                currentOut.close();

                // Create and write reverse complement
                std::string rcSequence = reverseComplement(currentSequence);
                fs::path rcContigPath = tempDir / (currentContig + "_rc.fa");
                contigFiles.emplace_back(rcContigPath.string(), true);
                currentRcOut.open(rcContigPath);
                if (!currentRcOut.is_open()) {
                    throw std::runtime_error("Could not create RC contig file: " + rcContigPath.string());
                }
                currentRcOut << ">" << currentContig << "_rc\n" << rcSequence << "\n";
                currentRcOut.close();

                currentSequence.clear();
            }

            // Start new contig
            size_t spacePos = line.find(' ');
            std::string contigName = (spacePos == std::string::npos) ?
                line.substr(1) :
                line.substr(1, spacePos - 1);

            // Create original contig file
            fs::path contigPath = tempDir / (contigName + ".fa");
            contigFiles.emplace_back(contigPath.string(), false);
            currentOut.open(contigPath);
            if (!currentOut.is_open()) {
                throw std::runtime_error("Could not create contig file: " + contigPath.string());
            }
            currentOut << line << "\n";
            currentContig = contigName;
        }
        else if (!currentContig.empty()) {
            // Remove any whitespace from sequence line
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            currentSequence += line;
            currentOut << line << "\n";
        }
    }

    // Process the last contig
    if (!currentContig.empty()) {
        currentOut.close();

        // Create and write reverse complement for last contig
        std::string rcSequence = reverseComplement(currentSequence);
        fs::path rcContigPath = tempDir / (currentContig + "_rc.fa");
        contigFiles.emplace_back(rcContigPath.string(), true);
        currentRcOut.open(rcContigPath);
        if (currentRcOut.is_open()) {
            currentRcOut << ">" << currentContig << "_rc\n" << rcSequence << "\n";
            currentRcOut.close();
        }
    }

    fastaFile.close();

    // Step 3: Create Mash sketches for each contig (original and RC)
    for (const auto& [contigFile, isRc] : contigFiles) {
        std::string mashCmd = "mash sketch -k " + std::to_string(kmerSize) +
            " -s " + std::to_string(sketchSize) +
            " -o " + contigFile + ".msh " + contigFile;
        int ret = system(mashCmd.c_str());
        if (ret != 0) {
            std::cerr << "Warning: Mash sketch failed for " << contigFile << std::endl;
        }
    }

    // Step 4: Compare all pairs of contigs using Mash
    for (size_t i = 0; i < contigFiles.size(); ++i) {
        const auto& [contigFile1, isRc1] = contigFiles[i];
        std::string contig1 = fs::path(contigFile1).stem().string();

        // Remove "_rc" suffix if present
        bool queryIsRc = contig1.size() > 3 && contig1.substr(contig1.size() - 3) == "_rc";
        std::string contig1Base = queryIsRc ? contig1.substr(0, contig1.size() - 3) : contig1;

        for (size_t j = 0; j < contigFiles.size(); ++j) {
            if (i == j) continue;

            const auto& [contigFile2, isRc2] = contigFiles[j];
            std::string contig2 = fs::path(contigFile2).stem().string();

            // Skip comparing a contig with its own RC
            std::string contig2Base = (contig2.size() > 3 && contig2.substr(contig2.size() - 3) == "_rc") ?
                contig2.substr(0, contig2.size() - 3) : contig2;
            if (contig1Base == contig2Base) continue;

            // Run Mash dist
            std::string mashDistCmd = "mash dist " + contigFile1 + ".msh " +
                contigFile2 + ".msh";

            // Capture output
            FILE* pipe = popen(mashDistCmd.c_str(), "r");
            if (!pipe) continue;

            char buffer[128];
            std::string result;
            while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                result += buffer;
            }
            pclose(pipe);

            // Parse Mash output
            std::istringstream iss(result);
            std::string ref1, ref2;
            double dist, pvalue;
            int sharedHashes;

            if (iss >> ref1 >> ref2 >> dist >> pvalue >> sharedHashes) {
                // Calculate containment
                double containment = static_cast<double>(sharedHashes) / sketchSize;

                if (containment >= threshold) {
                    // Determine if this is a reverse complement match
                    bool isReverseMatch = (contig2.size() > 3 && contig2.substr(contig2.size() - 3) == "_rc");

                    // Store the relationship (using base contig names without _rc suffix)
                    containmentMap[contig1Base].emplace_back(contig2Base, isReverseMatch);
                }
            }
        }
    }

    // Step 5: Clean up temporary files (optional)
    // fs::remove_all(tempDir);

    return containmentMap;
}

// Example usage:
int main() {
    try {
        auto results = checkContigContainmentWithRC("graph.final.fasta");

        for (const auto& [contig, contained] : results) {
            if (!contained.empty()) {
                std::cout << "Contig " << contig << " contains:\n";
                for (const auto& [c, isRc] : contained) {
                    std::cout << "  - " << c << (isRc ? " (reverse complement)" : "") << "\n";
                }
            }
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}