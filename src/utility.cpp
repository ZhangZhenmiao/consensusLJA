#include "dot_graph.hpp"
#include <sstream>
#include <iomanip>

std::string Graph::doubleToString(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << value;
    return stream.str();
}

template<typename T>
void Node::mergeMaps(std::unordered_map<std::string, std::vector<T>>& map1, const std::unordered_map<std::string, std::vector<T>>& map2) {
    for (const auto& pair : map2) {
        if (map1.find(pair.first) != map1.end()) {
            map1[pair.first].insert(map1[pair.first].end(), pair.second.begin(), pair.second.end());
        }
        else {
            map1[pair.first] = pair.second;
        }
    }
}

template void Graph::merge_vecs<std::string>(std::vector<std::string>&, std::vector<std::string>&);
template void Graph::merge_vecs<Edge>(std::vector<Edge>&, std::vector<Edge>&);

template<typename T>
void Graph::merge_vecs(std::vector<T>& e1, std::vector<T>& e2) {
    e1.insert(e1.end(), e2.begin(), e2.end());
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