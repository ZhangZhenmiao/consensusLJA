#include<string>

class DBGRunner {
public:
    std::string dot, fasta, aln, output;
    double low = 0;
    DBGRunner(std::string dot, std::string fasta, std::string aln, std::string output) {
        dot = dot;
        fasta = fasta;
        aln = aln;
        output = output;
        cleanDBG();
    }
    void cleanDBG();
};