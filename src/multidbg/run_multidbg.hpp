#include<string>

class MDBGRunner {
public:
    std::string dot, fasta, graph_dbg, paths_dbg, nodes, output;
    MDBGRunner(std::string dot, std::string fasta, std::string graph_dbg, std::string paths_dbg, std::string nodes, std::string output) {
        dot = dot;
        fasta = fasta;
        graph_dbg = graph_dbg;
        paths_dbg = paths_dbg;
        nodes = nodes;
        output = output;
        simplifyMDBG();
    }
    void simplifyMDBG();
};