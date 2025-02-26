#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <unistd.h>

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

int main(int argc, char* argv[]) {
    std::string prefix = argv[1];
    // std::string ref_seq = "/Poppy/zmzhang/Rust_fungi/genome/reference.compressed.only_chrs.fasta";
    std::string ref_seq = "/Poppy/zmzhang/cLJA_Project/Bonobo/genome/mPanPan1.compressed.fasta";
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
    if (system((exeDir + "/get_reference.py -o " + prefix + ".ref.bam.stats " + prefix + ".ref.bam " + prefix + ".fasta").c_str()) != 0)
        exit(1);
}