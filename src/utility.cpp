#include "dot_graph.hpp"
#include <sstream>
#include <iomanip>

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