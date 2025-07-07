#include "dot_graph.hpp"
#include <string>
#include <iostream>
#include <sys/stat.h>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <filesystem>
#include <map>
#include "src/dbg/dot_graph.hpp"
#include "src/dbg/read_aln.hpp"
namespace fs = std::filesystem;

using namespace multidbg;

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

void Graph::read_graph(std::string& output, std::string& graph_dot, const std::string& graph_fasta, const std::string& nodes_fasta, const std::string& graph_dbg, const std::string& paths_dbg) {
    fs::create_directory(output);
    if (graph.empty())
        read_from_dot(graph_dot, graph_fasta, nodes_fasta, graph_dbg, paths_dbg);
}

// Create histogram from edge multiplicities (bin size = 1)
std::map<int, int> Graph::create_histogram(const std::vector<double>& edges) {
    std::map<int, int> hist;
    for (double val : edges) {
        int bin = static_cast<int>(val);
        hist[bin]++;
    }
    return hist;
}

// Find peaks and valleys in histogram
void Graph::analyze_histogram(const std::map<int, int>& hist) {
    int window_size = 5;
    if (window_size < 3 || window_size % 2 == 0) {
        std::cerr << "Window size must be odd and >= 3" << std::endl;
        return;
    }

    int k = window_size / 2;

    if (hist.size() < window_size) {
        std::cerr << "Not enough data points for the chosen window size" << std::endl;
        return;
    }

    std::vector<std::pair<int, int>> peaks;
    std::vector<std::pair<int, int>> valleys;

    // Convert map to sorted vector for index-based access
    std::vector<std::pair<int, int>> vec(hist.begin(), hist.end());

    int max_center = 0;
    for (size_t i = k; i < vec.size() - k; ++i) {
        int center = vec[i].second;
        bool is_peak = true;
        bool is_valley = true;

        if (center > max_center) {
            mean_cov = vec[i].first;
            max_center = center;
        }

        for (int j = 1; j <= k; ++j) {
            if (center <= vec[i - j].second || center <= vec[i + j].second) {
                is_peak = false;
            }
            if (center >= vec[i - j].second || center >= vec[i + j].second) {
                is_valley = false;
            }
            // Early break for performance
            if (!is_peak && !is_valley) break;
        }

        if (is_peak) {
            peaks.emplace_back(vec[i].first, center);
        }
        else if (is_valley) {
            valleys.emplace_back(vec[i].first, center);
        }
    }

    // Output results
    if (!peaks.empty()) {
        // std::cout << "First peak at multiplicity: " << peaks[0].first
        //     << " (count: " << peaks[0].second << ")\n";

        error_peak = peaks[0].first;

        // std::cout << "Subsequent local minima: ";
        for (const auto& [valley, count] : valleys) {
            if (valley > peaks[0].first) {
                // std::cout << valley << " (count: " << count << ")\n";
                first_minima = valley;
                break;
            }
        }
        // std::cout << "Average coverage: " << mean_cov << std::endl;
    }
    else {
        // std::cout << "No peaks found in histogram" << std::endl;
    }
}

// Read graph from DOT file
void Graph::read_from_dot(const std::string& graph_dot, const std::string& graph_fasta, const std::string& nodes_fasta, const std::string& graph_dbg, const std::string& paths_dbg) {
    // load fasta sequence
    std::string e_id;
    std::string line, node1, node2, length, multiplicity, start_base;
    std::unordered_map<std::string, std::string> edge2sequence, node2sequence, nodenew2sequence;
    std::unordered_map<std::string, double> edgedbg2multi, edge2multi;
    std::unordered_map<std::string, double> edgedbg2len;
    std::unordered_map<std::string, std::pair<std::string, std::string>> edgedbg2nodes;

    // load edge multiplicities of dbg
    std::ifstream graph_dbg_file(graph_dbg);
    while (getline(graph_dbg_file, line)) {
        if (line.find("->") != std::string::npos) {
            size_t pos1 = line.find("->");
            std::string start_name = line.substr(1, pos1 - 3);
            size_t pos2 = line.find('[');
            std::string end_name = line.substr(pos1 + 4, pos2 - pos1 - 6);

            // parse edge label: starting base, length, multiplicity
            pos1 = line.find(")\" ");
            std::string label_all = line.substr(pos2 + 8, pos1 - pos2 - 7);
            std::string edge_label = label_all.substr(0, label_all.find(' '));
            label_all = label_all.substr(label_all.find(' ') + 1);
            pos1 = label_all.find('(');
            assert(pos1 != std::string::npos);
            char start_base = label_all[0];
            unsigned length = std::atoi(label_all.substr(2, pos1 - 2).c_str());
            double multiplicity = std::atof(label_all.substr(pos1 + 1, label_all.size() - pos1 - 2).c_str());
            edgedbg2multi[edge_label] = multiplicity;
            edgedbg2len[edge_label] = length;
            edgedbg2nodes[edge_label] = { start_name, end_name };
            // std::cout << "Multi for " << edge_label << " is " << multiplicity << ", len " << length << std::endl;
        }
    }
    graph_dbg_file.close();

    // read edge sequences of multidbg
    std::ifstream fasta_file(graph_fasta);
    while (getline(fasta_file, line)) {
        // line is contig name
        if (line.at(0) == '>') {
            e_id = line.substr(1);
        }
        // line is a contig
        else {
            std::string edge_sequence = line;
            assert(edge2sequence.find(e_id) == edge2sequence.end());
            edge2sequence[e_id] = edge_sequence;
            // std::cout << "Read " << e_id << " from " << graph_fasta << std::endl;
        }
    }

    // read node sequences of multidbg
    std::ifstream nodes_file(nodes_fasta);
    std::string node;
    while (getline(nodes_file, line)) {
        // line is contig name
        if (line.at(0) == '>') {
            node = line.substr(1);
        }
        // line is a contig
        else {
            std::string node_sequence = line;
            assert(node2sequence.find(node) == node2sequence.end());
            node2sequence[node] = node_sequence;
            // std::cout << "Read " << node << " from " << nodes_fasta << std::endl;
        }
    }

    // construct new node names of multidbg with reverse complementary information
    std::unordered_map<std::string, std::string> idMapping;
    std::unordered_set<std::string> nodes_without_rc;
    for (auto it = node2sequence.begin(); it != node2sequence.end(); ) {
        std::string currentId = it->first;
        std::string sequence = it->second;
        std::string revComp = reverse_complementary(sequence);

        // Find reverse complement in the map
        auto revIt = std::find_if(node2sequence.begin(), node2sequence.end(),
            [&revComp](const std::pair<std::string, std::string>& pair) {
                return pair.second == revComp;
            });

        if (revIt != node2sequence.end() && revIt->first != currentId) {
            // Rename current ID and reverse complement ID
            idMapping[currentId] = currentId;
            idMapping[revIt->first] = "-" + currentId;
            nodeid2Rev[currentId] = "-" + currentId;
            nodeid2Rev["-" + currentId] = currentId;
            nodenew2sequence[currentId] = sequence;
            nodenew2sequence["-" + currentId] = revComp;

            std::cout << "Node mapping: " << currentId << " to " << currentId << std::endl;
            std::cout << "Node mapping: " << revIt->first << " to " << "-" + currentId << std::endl;

            // Remove reverse complement from map
            node2sequence.erase(revIt);
        }
        else {
            // If no reverse complement is found, keep the ID as is
            if (revIt == node2sequence.end()) {
                std::cout << "Node cannot find reverse complement: " << currentId << std::endl;
                idMapping[currentId] = currentId;
                nodeid2Rev[currentId] = "-" + currentId;
                nodeid2Rev["-" + currentId] = currentId;
                nodenew2sequence[currentId] = sequence;
                nodenew2sequence["-" + currentId] = revComp;
                // this container stores initial node ids for nodes without rc; the new node names should be in idMapping
                nodes_without_rc.insert(currentId);
            }
            else {
                std::cout << "Node is self-reverse: " << currentId << std::endl;
                idMapping[currentId] = currentId;
                nodeid2Rev[currentId] = currentId;
                nodenew2sequence[currentId] = sequence;
            }
        }

        // Remove current ID from map
        it = node2sequence.erase(it);
    }

    //load multidbg edge paths in dbg, and construct coverages for multidbg edges
    std::ifstream paths_dbg_file(paths_dbg);
    std::string edge_name;
    int max_initial_edge_id = 0;
    while (getline(paths_dbg_file, line)) {
        if (line.at(0) == '>') {
            edge_name = line.substr(1);
            int edge_id = std::atoi(edge_name.c_str());
            if (edge_id > max_initial_edge_id)
                max_initial_edge_id = edge_id;
        }
        else {
            // std::cout << edge_name << ": Start";
            std::vector<double> multis;
            std::vector<double> multis_long;
            // std::cout << "Path for " << edge_name << ": ";
            while (line.find(' ') != std::string::npos) {
                std::string edge_id = line.substr(0, line.find(' '));
                assert(edgedbg2multi.find(edge_id) != edgedbg2multi.end());
                multis.push_back(edgedbg2multi[edge_id]);
                if (edgedbg2len[edge_id] >= 100000)
                    multis_long.push_back(edgedbg2multi[edge_id]);
                line = line.substr(line.find(' ') + 1);

                assert(edgedbg2nodes.find(edge_id) != edgedbg2nodes.end());
                if (edge_initial_to_path_in_dbg[edge_name].empty()) {
                    edge_initial_to_path_in_dbg[edge_name].push_back(edgedbg2nodes[edge_id].first);
                    // std::cout << edgedbg2nodes[edge_id].first;
                }
                assert(edgedbg2nodes[edge_id].first == edge_initial_to_path_in_dbg[edge_name].at(edge_initial_to_path_in_dbg[edge_name].size() - 1));
                edge_initial_to_path_in_dbg[edge_name].push_back(edgedbg2nodes[edge_id].second);

                // std::cout << " (" << edge_id << ")";
                // std::cout << " " << edgedbg2nodes[edge_id].second;
            }
            // std::cout << " path len " << edge_initial_to_DBG_reads[edge_name].dbg_path_sequence.size() << " mdbg edge len " << edge2sequence[edge_name].size() << std::endl;

            if (multis_long.size()) {
                double min_multi = 0;
                std::sort(multis_long.begin(), multis_long.end());
                for (auto&& m : multis_long) {
                    if (m > 0) {
                        min_multi = m;
                        break;
                    }
                }
                if (min_multi == 0)
                    min_multi = multis_long.at(multis_long.size() - 1);

                // std::cout << "; Min multi for " << edge_name << " is " << min_multi << std::endl;
                edge2multi[edge_name] = min_multi;
            }
            else {
                double min_multi = 0;
                std::sort(multis.begin(), multis.end());
                for (auto&& m : multis) {
                    if (m > 0) {
                        min_multi = m;
                        break;
                    }
                }
                if (min_multi == 0)
                    min_multi = multis.at(multis.size() - 1);

                // std::cout << "; Min multi for " << edge_name << " is " << min_multi << std::endl;
                edge2multi[edge_name] = min_multi;
            }
        }
    }
    paths_dbg_file.close();

    // load multidbg graph
    long cnt_edge = 0;
    std::vector<double> edge_multis;
    std::ifstream dot_file(graph_dot);
    while (getline(dot_file, line)) {
        // line is neithor a node nor an edge
        if (line.find('[') == std::string::npos)
            continue;
        // line is an edge, in this case, all the node should already be loaded
        if (line.find("->") != std::string::npos) {
            // extract node name
            size_t pos1 = line.find("->");
            std::string start_name = idMapping.at(line.substr(1, pos1 - 1));
            size_t pos2 = line.find('[');
            std::string end_name = idMapping.at(line.substr(pos1 + 2, pos2 - pos1 - 2));
            assert(this->graph.find(start_name) != this->graph.end() && this->graph.find(end_name) != this->graph.end());

            // parse edge label: starting base, length, multiplicity and sequence
            pos1 = line.find("size");
            std::string edge_label = line.substr(pos2 + 14, pos1 - pos2 - 16);

            // std::cout << "Read edge " << start_name << " -> " << end_name << " from " << graph_dot << std::endl;

            assert(edge2sequence.find(edge_label) != edge2sequence.end());
            assert(edge2multi.find(edge_label) != edge2multi.end());
            int length = edge2sequence[edge_label].size();
            char start_base = edge2sequence[edge_label].at(graph[start_name].sequence.size());
            Edge edge = Edge(start_base, length, edge2sequence[edge_label], edge2multi[edge_label]);
            assert(edge2sequence[edge_label].substr(0, graph[start_name].sequence.size()) == graph[start_name].sequence);
            edge.path_nodes_in_original_graph.push_back(start_name);
            edge.path_nodes_in_original_graph.push_back(end_name);
            edge.path_edges_in_original_graph.push_back(edge_label);
            // don't set edge label as initial labels, as we use a different naming rule for edge labels
            // edge.label = edge_label;

            this->graph[start_name].outgoing_edges[end_name].push_back(edge);
            this->graph[end_name].incoming_edges[start_name].push_back(edge);
            cnt_edge += 1;
            edge_multis.push_back(edge2multi[edge_label]);

            // normally reverse complements edges are stored as two lines, but for edges on nodes without rc, there is only one line, so we should add a dummy edge to keep the graph symmetric
            if (nodes_without_rc.find(start_name) != nodes_without_rc.end() || nodes_without_rc.find(end_name) != nodes_without_rc.end()) {
                std::string start_r = reverse_complementary_node(end_name);
                std::string end_r = reverse_complementary_node(start_name);

                std::string seq_r = reverse_complementary(edge2sequence[edge_label]);
                int length = seq_r.size();
                char start_base = seq_r.at(graph[start_r].sequence.size());
                Edge edge = Edge(start_base, length, seq_r, edge2multi[edge_label]);

                std::string new_edge_initial_label = std::to_string(++max_initial_edge_id);
                edge.path_nodes_in_original_graph.push_back(start_r);
                edge.path_nodes_in_original_graph.push_back(end_r);
                edge.path_edges_in_original_graph.push_back(new_edge_initial_label);

                // need to construct the dbg path for the newly added edge initial id, this is for fianl containment checking
                if (edge_initial_to_path_in_dbg[edge_label].size()) {
                    for (int i = int(edge_initial_to_path_in_dbg[edge_label].size()) - 1; i >= 0; --i) {
                        std::string node_in_dbg = edge_initial_to_path_in_dbg[edge_label].at(i);
                        std::string rc_node_in_dbg = node_in_dbg.at(0) != '-' ? "-" + node_in_dbg : node_in_dbg.substr(1);
                        edge_initial_to_path_in_dbg[new_edge_initial_label].push_back(rc_node_in_dbg);
                    }
                }

                this->graph[start_r].outgoing_edges[end_r].push_back(edge);
                this->graph[end_r].incoming_edges[start_r].push_back(edge);

                std::cout << "Add complementary edge " << start_r << " -> " << end_r << std::endl;

                cnt_edge += 1;
            }
        }
        // line is an node
        else {
            std::string node_name = line.substr(1, line.find('[') - 1);
            Node node;
            assert(idMapping.find(node_name) != idMapping.end());
            node.sequence = nodenew2sequence.at(idMapping.at(node_name));
            this->graph[idMapping.at(node_name)] = node;

            // normally reverse complements nodes are stored as two lines, but for nodes without rc, we should add a dummy node
            if (nodes_without_rc.find(node_name) != nodes_without_rc.end()) {
                Node node;
                node.sequence = nodenew2sequence.at(reverse_complementary_node(idMapping.at(node_name)));
                this->graph[reverse_complementary_node(idMapping.at(node_name))] = node;
                std::cout << "Add complementary node " << reverse_complementary_node(idMapping.at(node_name)) << std::endl;
            }
        }
    }
    dot_file.close();

    auto hist = create_histogram(edge_multis);
    analyze_histogram(hist);
    std::cout << "Read " << get_num_nodes() << " vertices, " << cnt_edge << " edges." << std::endl;
}

std::string Graph::get_unique_label(std::unordered_set<std::string>& labels) {
    int number = 10 + std::rand() % 90;
    while (true) {
        std::string numb = std::to_string(number);
        if (labels.find(numb) == labels.end())
            return numb;
        number++;
    }
}

std::string Graph::get_contracted_name(std::string node) {
    return node.substr(0, node.find_first_of('_')) + "+" + node.substr(node.find_last_of('_') + 1);
}

std::string Graph::get_contracted_label(std::string node) {
    std::string label = node.substr(0, node.find_first_of('_')) + "_" + node.substr(node.find_last_of('_') + 1) + "_N" + std::to_string(graph[node].number_of_contracted_edge) + "_L" + std::to_string(graph[node].length_of_contracted_edge);
    if (graph[node].circles >= 2) {
        label += ("\\nC" + std::to_string(graph[node].circles));
        label += ("_CL" + std::to_string(graph[node].total_length));
        label += ("_CM" + std::to_string(int(graph[node].median_length)));
    }
    return label;
}

void Graph::write_graph(const std::string& prefix, int thick, bool contracted, bool colored, std::unordered_set<std::string> nodes_retain) {
    std::string graph_dot = prefix + ".dot";
    std::string graph_fasta = prefix + ".fasta";
    std::string graph_path = prefix + ".path";

    std::cout << "Write graph " << graph_dot << ", fasta " << graph_fasta << std::endl;
    std::ofstream file_dot(graph_dot);
    std::ofstream file_fasta(graph_fasta);
    std::ofstream file_path(graph_path);

    int num_edges = 0;
    file_dot << "digraph {\nnodesep = 0.5;\n";
    int max_contracted = 0;
    std::string max_contracted_node;
    for (auto&& node : this->graph) {
        if (!nodes_retain.empty() && nodes_retain.find(node.first) == nodes_retain.end())
            continue;
        if (node.second.number_of_contracted_edge > 0) {
            std::string node_o = get_contracted_name(node.first);
            nodeid2Rev[node_o] = get_contracted_name(reverse_complementary_node(node.first));
            nodeid2Rev[get_contracted_name(reverse_complementary_node(node.first))] = get_contracted_name(node.first);
            file_dot << "\"" << node_o << "\" [style=filled fillcolor=\"white\" label=\"" << get_contracted_label(node.first) << "\"]\n";
        }
        else
            file_dot << node.first << " [style=filled fillcolor=\"white\" label=\"" << node.first + "_L" + std::to_string(graph[node.first].sequence.size()) << "\"]\n";
        if (node.second.number_of_contracted_edge > max_contracted) {
            max_contracted_node = get_contracted_name(node.first);
            max_contracted = node.second.number_of_contracted_edge;
        }
    }
    if (max_contracted > 0)
        std::cout << "Max contracted node: " << max_contracted_node << " has " << max_contracted << " edges." << std::endl;

    // construct labels for vertices
    std::unordered_map<std::string, std::unordered_set<std::string>> vertice2labels;
    for (auto&& node : this->graph) {
        std::string start_node = node.first;
        for (auto&& edges : node.second.outgoing_edges) {
            for (auto&& edge : edges.second) {
                if (!edge.label.empty())
                    vertice2labels[start_node].insert(edge.label.substr(edge.label.find('.') + 1));
            }
        }
    }

    std::unordered_map<std::string, std::string> color_map = {
        {"1A", "#325527"},
        {"1B", "#325527"},
        {"2A", "#628DCF"},
        {"2B", "#628DCF"},
        {"3A", "#41496B"},
        {"3B", "#41496B"},
        {"4A", "#12CCD6"},
        {"4B", "#12CCD6"},
        {"5A", "#3E16F3"},
        {"5B", "#3E16F3"},
        {"6A", "#E46C0A"},
        {"6B", "#E46C0A"},
        {"7A", "#446768"},
        {"7B", "#446768"},
        {"8A", "#FF0000"},
        {"8B", "#FF0000"},
        {"9A", "#3C06A6"},
        {"9B", "#3C06A6"},
        {"10A", "#6CB9AB"},
        {"10B", "#6CB9AB"},
        {"11A", "#988430"},
        {"11B", "#988430"},
        {"12A", "#4BAA54"},
        {"12B", "#4BAA54"},
        {"13A", "#154E54"},
        {"13B", "#154E54"},
        {"14A", "#A74C5D"},
        {"14B", "#A74C5D"},
        {"15A", "#528444"},
        {"15B", "#528444"},
        {"16A", "#B61664"},
        {"16B", "#B61664"},
        {"17A", "#8F3296"},
        {"17B", "#8F3296"},
        {"18A", "#E1A9E7"},
        {"18B", "#E1A9E7"},
        {"19A", "#54340D"},
        {"19B", "#54340D"},
        {"20A", "#316260"},
        {"20B", "#316260"},
        {"21A", "#8041AF"},
        {"21B", "#8041AF"},
        {"22A", "#5AB499"},
        {"22B", "#5AB499"},
        {"23A", "#952395"},
        {"23B", "#952395"},
        {"24A", "#70229F"},
        {"24B", "#70229F"},
        {"25A", "#4D4050"},
        {"25B", "#4D4050"},
        {"26A", "#969696"},
        {"26B", "#969696"},
        {"1M", "#325527"},
        {"1P", "#325527"},
        {"2M", "#628DCF"},
        {"2P", "#628DCF"},
        {"3M", "#41496B"},
        {"3P", "#41496B"},
        {"4M", "#12CCD6"},
        {"4P", "#12CCD6"},
        {"5M", "#3E16F3"},
        {"5P", "#3E16F3"},
        {"6M", "#E46C0A"},
        {"6P", "#E46C0A"},
        {"7M", "#446768"},
        {"7P", "#446768"},
        {"8M", "#FF0000"},
        {"8P", "#FF0000"},
        {"9M", "#3C06A6"},
        {"9P", "#3C06A6"},
        {"10M", "#6CB9AB"},
        {"10P", "#6CB9AB"},
        {"11M", "#988430"},
        {"11P", "#988430"},
        {"12M", "#4BAA54"},
        {"12P", "#4BAA54"},
        {"13M", "#154E54"},
        {"13P", "#154E54"},
        {"14M", "#A74C5D"},
        {"14P", "#A74C5D"},
        {"15M", "#528444"},
        {"15P", "#528444"},
        {"16M", "#B61664"},
        {"16P", "#B61664"},
        {"17M", "#8F3296"},
        {"17P", "#8F3296"},
        {"18M", "#E1A9E7"},
        {"18P", "#E1A9E7"},
        {"19M", "#54340D"},
        {"19P", "#54340D"},
        {"20M", "#316260"},
        {"20P", "#316260"},
        {"21M", "#8041AF"},
        {"21P", "#8041AF"},
        {"22M", "#5AB499"},
        {"22P", "#5AB499"},
        {"23M", "#952395"},
        {"23P", "#952395"},
        {"X", "#969696"},
        {"Y", "#969696"},
        {"mtDNA", "#FF0000"},
        {"Chr1", "#325527"},
        {"Chr2", "#628DCF"},
        {"Chr3", "#41496B"},
        {"Chr4", "#12CCD6"},
        {"Chr5", "#3E16F3"},
        {"Chr6", "#E46C0A"},
        {"Chr7", "#446768"},
        {"Chr8", "#FF0000"},
        {"Chr9", "#3C06A6"},
        {"Chr10", "#6CB9AB"},
        {"Chr11", "#988430"},
        {"Chr12", "#4BAA54"},
        {"Chr13", "#154E54"},
        {"Chr14", "#A74C5D"}
    };

    std::unordered_set <std::string> traversed_labels;
    for (auto&& node : this->graph) {
        std::string start_node = node.first;
        if (node.second.number_of_contracted_edge > 0) {
            start_node = get_contracted_name(node.first);
        }
        for (auto&& edges : node.second.outgoing_edges) {
            std::string end_node = edges.first;
            if (!nodes_retain.empty() && nodes_retain.find(node.first) == nodes_retain.end() && nodes_retain.find(end_node) == nodes_retain.end())
                continue;
            if (graph[edges.first].number_of_contracted_edge > 0) {
                end_node = get_contracted_name(edges.first);
            }
            for (auto&& edge : edges.second) {
                if (traversed_labels.find(edge.label) != traversed_labels.end()) {
                    num_edges += 1;
                    if (contracted || colored) {
                        if (edge.ref_ids.empty())
                            file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << static_cast<int>(round(edge.multiplicity)) << ")\" color=\"black\"]\n";
                        else {
                            file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << static_cast<int>(round(edge.multiplicity)) << ")";
                            for (int i = 0; i < edge.ref_ids.size(); ++i) {
                                file_dot << "\\n" << edge.ref_ids[i];
                            }
                            std::string chr = edge.ref_ids.at(0).substr(0, edge.ref_ids.at(0).find(' '));
                            if (chr.at(0) == '-') chr = chr.substr(1);
                            std::string color = color_map[chr];
                            file_dot << "\" color=\"" << color << "\"]\n";
                        }
                    }
                    else
                        file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << edge.multiplicity << ")\" color=\"black\"]\n";
                    continue;
                }
                if (!edge.label.empty()) {
                    traversed_labels.insert(edge.rc_label);
                }
                else {
                    std::unordered_set<std::string>& forward_labels = vertice2labels[start_node];
                    std::unordered_set<std::string>& reverse_labels = vertice2labels[reverse_complementary_node(end_node)];
                    std::string label_forward = start_node + '.' + get_unique_label(forward_labels), label_reverse = reverse_complementary_node(end_node) + '.' + get_unique_label(reverse_labels);
                    if (edge.sequence == reverse_complementary(edge.sequence)) {
                        assert(end_node == reverse_complementary_node(start_node));
                        label_reverse = label_forward;
                    }
                    edge.label = label_forward;
                    edge.rc_label = label_reverse;
                    // std::cout << "New label " << label_forward << " and " << label_reverse << std::endl;
                    bool flag = false;
                    for (auto&& edge_i : graph[edges.first].incoming_edges[node.first]) {
                        if (edge_i.sequence == edge.sequence) {
                            // the edge should appear only once
                            assert(flag == false);
                            flag = true;
                            edge_i.label = label_forward;
                            edge_i.rc_label = label_reverse;
                        }
                    }
                    assert(flag);
                    flag = false;
                    for (auto&& edge_r : graph[reverse_complementary_node(edges.first)].outgoing_edges[reverse_complementary_node(node.first)]) {
                        if (edge_r.sequence == reverse_complementary(edge.sequence) || graph[reverse_complementary_node(edges.first)].outgoing_edges[reverse_complementary_node(node.first)].size() == 1) {
                            assert(flag == false);
                            flag = true;
                            edge_r.rc_label = label_forward;
                            edge_r.label = label_reverse;
                        }
                    }
                    assert(flag);
                    flag = false;
                    for (auto&& edge_r_i : graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(edges.first)]) {
                        if (edge_r_i.sequence == reverse_complementary(edge.sequence) || graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(edges.first)].size() == 1) {
                            assert(flag == false);
                            flag = true;
                            edge_r_i.rc_label = label_forward;
                            edge_r_i.label = label_reverse;
                        }
                    }
                    assert(flag);

                    traversed_labels.insert(label_reverse);
                    vertice2labels[start_node].insert(label_forward.substr(label_forward.find('.') + 1));
                    vertice2labels[reverse_complementary_node(end_node)].insert(label_reverse.substr(label_reverse.find('.') + 1));
                }

                num_edges += 1;
                if (contracted || colored) {
                    if (edge.ref_ids.empty())
                        file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << static_cast<int>(round(edge.multiplicity)) << ")\" color=\"black\"]\n";
                    else {
                        file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << static_cast<int>(round(edge.multiplicity)) << ")";
                        for (int i = 0; i < edge.ref_ids.size(); ++i) {
                            file_dot << "\\n" << edge.ref_ids[i];
                        }
                        std::string chr = edge.ref_ids.at(0).substr(0, edge.ref_ids.at(0).find(' '));
                        if (chr.at(0) == '-') chr = chr.substr(1);
                        std::string color = color_map[chr];
                        file_dot << "\" color=\"" << color << "\"]\n";
                    }
                }
                else
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << edge.multiplicity << ")\" color=\"black\"]\n";
                file_fasta << ">" << edge.label << "_" << edge.rc_label << "\n";
                file_fasta << edge.sequence << "\n";
                file_path << ">" << edge.label << "_" << edge.rc_label << " " << edge.length << "\n";
                assert(edge.path_nodes_in_original_graph.size() == edge.path_edges_in_original_graph.size() + 1);
                if (edge.path_nodes_in_original_graph.size() >= 1)
                    file_path << edge.path_nodes_in_original_graph.at(0);
                for (int i = 1; i < edge.path_nodes_in_original_graph.size(); ++i)
                    file_path << "->(" << edge.path_edges_in_original_graph.at(i - 1) << ")->" << edge.path_nodes_in_original_graph.at(i);
                file_path << "\n";
            }
        }
    }
    file_dot << "}" << std::endl;
    file_dot.close();
    file_fasta.close();
    file_path.close();
    std::cout << "Total number of nodes: " << this->get_num_nodes() << std::endl;
    std::cout << "Total number of edges: " << num_edges << std::endl;
}

// assume no simple bubbles
void Graph::write_graph_contracted(const std::string& prefix, int min_length, bool simplify) {
    // std::cout << "----------Contracted visulization----------" << std::endl;
    std::unordered_map<std::string, Node> graph_vis = this->graph;
    struct Nodes_To_Contract
    {
        std::string node1;
        std::string node2;
        int edge_length;
        std::string edge_sequence;
        Nodes_To_Contract(std::string n1, std::string n2, int l, std::string seq) {
            node1 = n1;
            node2 = n2;
            edge_length = l;
            edge_sequence = seq;
        }
    };

    std::unordered_map<std::string, std::unordered_set<std::string>> node1_to_node2_scanned;
    std::vector<Nodes_To_Contract> nodes_to_contract;
    for (auto&& node : graph_vis) {
        // skip simple components
        std::unordered_set<std::string> connected_nodes;
        for (auto&& n : node.second.incoming_edges) {
            if (n.first != node.first)
                connected_nodes.insert(n.first);
        }
        for (auto&& n : node.second.outgoing_edges) {
            if (n.first != node.first)
                connected_nodes.insert(n.first);
        }
        if (connected_nodes.size() <= 1) {
            bool flag = true;
            for (auto&& c : connected_nodes) {
                for (auto&& n : graph_vis[c].outgoing_edges) {
                    if (n.first != node.first)
                        flag = false;
                }
                for (auto&& n : graph_vis[c].incoming_edges) {
                    if (n.first != node.first)
                        flag = false;
                }
            }
            if (flag)
                continue;
        }

        // filter edges to be contracted
        for (auto&& edge : node.second.outgoing_edges) {
            // do not deal with palindromic bulges
            if (node.first == reverse_complementary_node(edge.first) || node.first == edge.first)
                continue;
            if (node1_to_node2_scanned[node.first].find(edge.first) != node1_to_node2_scanned[node.first].end() || node1_to_node2_scanned[reverse_complementary_node(edge.first)].find(reverse_complementary_node(node.first)) != node1_to_node2_scanned[reverse_complementary_node(edge.first)].end())
                continue;
            for (int i = 0; i < edge.second.size();++i) {
                // the edge should be collapsed
                if (edge.second.at(i).length <= min_length) {
                    nodes_to_contract.emplace_back(Nodes_To_Contract(node.first, edge.first, edge.second.at(i).length, edge.second.at(i).sequence));
                    nodes_to_contract.emplace_back(Nodes_To_Contract(reverse_complementary_node(edge.first), reverse_complementary_node(node.first), edge.second.at(i).length, reverse_complementary(edge.second.at(i).sequence)));
                }
            }
            node1_to_node2_scanned[node.first].insert(edge.first);
            node1_to_node2_scanned[reverse_complementary_node(edge.first)].insert(reverse_complementary_node(node.first));
        }
    }

    auto process_contracting = [&](Nodes_To_Contract& node_pair, std::unordered_map<std::string, std::string>& node2contracted, bool add_seq = false) -> void
        {
            std::string node1 = node_pair.node1, node2 = node_pair.node2;
            while (node2contracted.find(node1) != node2contracted.end())
                node1 = node2contracted[node1];
            while (node2contracted.find(node2) != node2contracted.end())
                node2 = node2contracted[node2];

            if (node1 == node2) {
                graph_vis[node1].number_of_contracted_edge += 1;
                graph_vis[node1].length_of_contracted_edge += node_pair.edge_length;
                // std::cout << "Contract " << node1 << "(" << node_pair.node1 << ")" << " and " << node2 << "(" << node_pair.node2 << ") to " << node1 << ", number " << graph_vis[node1].number_of_contracted_edge << ", length " << graph_vis[node1].length_of_contracted_edge << std::endl;
                int i = 0;
                for (;i < graph_vis[node2].incoming_edges[node1].size();++i) {
                    if (graph_vis[node2].incoming_edges[node1][i].sequence == node_pair.edge_sequence || graph_vis[node2].incoming_edges[node1].size() == 1)
                        break;
                }
                assert(i != graph_vis[node2].incoming_edges[node1].size());
                graph_vis[node2].incoming_edges[node1].erase(graph_vis[node2].incoming_edges[node1].begin() + i);
                if (graph_vis[node2].incoming_edges[node1].empty())
                    graph_vis[node2].incoming_edges.erase(node1);

                int j = 0;
                for (;j < graph_vis[node1].outgoing_edges[node2].size();++j) {
                    if (graph_vis[node1].outgoing_edges[node2][j].sequence == node_pair.edge_sequence || graph_vis[node1].outgoing_edges[node2].size() == 1)
                        break;
                }
                assert(j != graph_vis[node1].outgoing_edges[node2].size());
                graph_vis[node1].outgoing_edges[node2].erase(graph_vis[node1].outgoing_edges[node2].begin() + j);
                if (graph_vis[node1].outgoing_edges[node2].empty())
                    graph_vis[node1].outgoing_edges.erase(node2);

                if (add_seq) {
                    if (graph_vis[node1].sequence.empty())
                        graph_vis[node1].sequence += ("NNNNNNNNNNNNNNNNNNNN" + node_pair.edge_sequence + "NNNNNNNNNNNNNNNNNNNN");
                    else
                        graph_vis[node1].sequence += (node_pair.edge_sequence + "NNNNNNNNNNNNNNNNNNNN");
                }

                return;
            }

            std::string node_contracted = node1 + "_" + node2;
            nodeid2Rev[node_contracted] = reverse_complementary_node(node2) + "_" + reverse_complementary_node(node1);
            nodeid2Rev[reverse_complementary_node(node2) + "_" + reverse_complementary_node(node1)] = node_contracted;
            node2contracted[node1] = node_contracted;
            node2contracted[node2] = node_contracted;

            graph_vis[node_contracted].number_of_contracted_edge = graph_vis[node1].number_of_contracted_edge + graph_vis[node2].number_of_contracted_edge + 1;
            graph_vis[node_contracted].length_of_contracted_edge = graph_vis[node1].length_of_contracted_edge + graph_vis[node2].length_of_contracted_edge + node_pair.edge_length;
            // std::cout << "Contract " << node1 << "(" << node_pair.node1 << ")" << " and " << node2 << "(" << node_pair.node2 << ") to " << node_contracted << ", number " << graph_vis[node_contracted].number_of_contracted_edge << ", length " << graph_vis[node_contracted].length_of_contracted_edge << std::endl;

            int i = 0;
            for (;i < graph_vis[node2].incoming_edges[node1].size();++i) {
                if (graph_vis[node2].incoming_edges[node1][i].sequence == node_pair.edge_sequence || graph_vis[node2].incoming_edges[node1].size() == 1)
                    break;
            }
            assert(i != graph_vis[node2].incoming_edges[node1].size());
            graph_vis[node2].incoming_edges[node1].erase(graph_vis[node2].incoming_edges[node1].begin() + i);
            if (graph_vis[node2].incoming_edges[node1].empty())
                graph_vis[node2].incoming_edges.erase(node1);

            int j = 0;
            for (;j < graph_vis[node1].outgoing_edges[node2].size();++j) {
                if (graph_vis[node1].outgoing_edges[node2][j].sequence == node_pair.edge_sequence || graph_vis[node1].outgoing_edges[node2].size() == 1)
                    break;
            }
            assert(j != graph_vis[node1].outgoing_edges[node2].size());
            graph_vis[node1].outgoing_edges[node2].erase(graph_vis[node1].outgoing_edges[node2].begin() + j);
            if (graph_vis[node1].outgoing_edges[node2].empty())
                graph_vis[node1].outgoing_edges.erase(node2);

            assert(i == j);

            if (graph_vis[node1].incoming_edges.find(node1) != graph_vis[node1].incoming_edges.end()) {
                merge_vecs(graph_vis[node1].incoming_edges[node_contracted], graph_vis[node1].incoming_edges[node1]);
                graph_vis[node1].incoming_edges.erase(node1);
                assert(graph_vis[node1].outgoing_edges.find(node1) != graph_vis[node1].outgoing_edges.end());
                merge_vecs(graph_vis[node1].outgoing_edges[node_contracted], graph_vis[node1].outgoing_edges[node1]);
                graph_vis[node1].outgoing_edges.erase(node1);
            }
            else
                assert(graph_vis[node1].outgoing_edges.find(node1) == graph_vis[node1].outgoing_edges.end());
            if (graph_vis[node1].incoming_edges.find(node2) != graph_vis[node1].incoming_edges.end()) {
                merge_vecs(graph_vis[node1].incoming_edges[node_contracted], graph_vis[node1].incoming_edges[node2]);
                graph_vis[node1].incoming_edges.erase(node2);
                assert(graph_vis[node2].outgoing_edges.find(node1) != graph_vis[node2].outgoing_edges.end());
                merge_vecs(graph_vis[node2].outgoing_edges[node_contracted], graph_vis[node2].outgoing_edges[node1]);
                graph_vis[node2].outgoing_edges.erase(node1);
            }
            else
                assert(graph_vis[node2].outgoing_edges.find(node1) == graph_vis[node2].outgoing_edges.end());
            if (graph_vis[node2].incoming_edges.find(node1) != graph_vis[node2].incoming_edges.end()) {
                merge_vecs(graph_vis[node2].incoming_edges[node_contracted], graph_vis[node2].incoming_edges[node1]);
                graph_vis[node2].incoming_edges.erase(node1);
                assert(graph_vis[node1].outgoing_edges.find(node2) != graph_vis[node1].outgoing_edges.end());
                merge_vecs(graph_vis[node1].outgoing_edges[node_contracted], graph_vis[node1].outgoing_edges[node2]);
                graph_vis[node1].outgoing_edges.erase(node2);
            }
            else
                assert(graph_vis[node1].outgoing_edges.find(node2) == graph_vis[node1].outgoing_edges.end());
            if (graph_vis[node2].incoming_edges.find(node2) != graph_vis[node2].incoming_edges.end()) {
                merge_vecs(graph_vis[node2].incoming_edges[node_contracted], graph_vis[node2].incoming_edges[node2]);
                graph_vis[node2].incoming_edges.erase(node2);
                assert(graph_vis[node2].outgoing_edges.find(node2) != graph_vis[node2].outgoing_edges.end());
                merge_vecs(graph_vis[node2].outgoing_edges[node_contracted], graph_vis[node2].outgoing_edges[node2]);
                graph_vis[node2].outgoing_edges.erase(node2);
            }
            else
                assert(graph_vis[node2].outgoing_edges.find(node2) == graph_vis[node2].outgoing_edges.end());

            graph_vis[node_contracted].mergeMaps(graph_vis[node_contracted].incoming_edges, graph_vis[node1].incoming_edges);
            graph_vis[node_contracted].mergeMaps(graph_vis[node_contracted].incoming_edges, graph_vis[node2].incoming_edges);
            graph_vis[node_contracted].mergeMaps(graph_vis[node_contracted].outgoing_edges, graph_vis[node1].outgoing_edges);
            graph_vis[node_contracted].mergeMaps(graph_vis[node_contracted].outgoing_edges, graph_vis[node2].outgoing_edges);

            for (auto&& node : graph_vis[node1].incoming_edges) {
                if (node.first == node_contracted)
                    continue;
                assert(node.first != node1 && node.first != node2);
                assert(graph_vis[node.first].outgoing_edges.find(node1) != graph_vis[node.first].outgoing_edges.end());
                merge_vecs(graph_vis[node.first].outgoing_edges[node_contracted], graph_vis[node.first].outgoing_edges[node1]);
                graph_vis[node.first].outgoing_edges.erase(node1);
            }
            for (auto&& node : graph_vis[node2].incoming_edges) {
                if (node.first == node_contracted)
                    continue;
                assert(node.first != node1 && node.first != node2);
                assert(graph_vis[node.first].outgoing_edges.find(node2) != graph_vis[node.first].outgoing_edges.end());
                merge_vecs(graph_vis[node.first].outgoing_edges[node_contracted], graph_vis[node.first].outgoing_edges[node2]);
                graph_vis[node.first].outgoing_edges.erase(node2);
            }
            for (auto&& node : graph_vis[node1].outgoing_edges) {
                if (node.first == node_contracted)
                    continue;
                assert(node.first != node1 && node.first != node2);
                assert(graph_vis[node.first].incoming_edges.find(node1) != graph_vis[node.first].incoming_edges.end());
                merge_vecs(graph_vis[node.first].incoming_edges[node_contracted], graph_vis[node.first].incoming_edges[node1]);
                graph_vis[node.first].incoming_edges.erase(node1);
            }
            for (auto&& node : graph_vis[node2].outgoing_edges) {
                if (node.first == node_contracted)
                    continue;
                assert(node.first != node1 && node.first != node2);
                assert(graph_vis[node.first].incoming_edges.find(node2) != graph_vis[node.first].incoming_edges.end());
                merge_vecs(graph_vis[node.first].incoming_edges[node_contracted], graph_vis[node.first].incoming_edges[node2]);
                graph_vis[node.first].incoming_edges.erase(node2);
            }

            graph_vis.erase(node1);
            graph_vis.erase(node2);

            for (auto&& node : graph_vis[node_contracted].outgoing_edges) {
                assert(graph_vis[node.first].incoming_edges[node_contracted].size() == graph_vis[node_contracted].outgoing_edges[node.first].size());
            }
            for (auto&& node : graph_vis[node_contracted].incoming_edges) {
                assert(graph_vis[node.first].outgoing_edges[node_contracted].size() == graph_vis[node_contracted].incoming_edges[node.first].size());
            }

            if (add_seq) {
                if (graph_vis[node_contracted].sequence.empty())
                    graph_vis[node_contracted].sequence += ("NNNNNNNNNNNNNNNNNNNN" + node_pair.edge_sequence + "NNNNNNNNNNNNNNNNNNNN");
                else
                    graph_vis[node_contracted].sequence += (node_pair.edge_sequence + "NNNNNNNNNNNNNNNNNNNN");
            }
        };

    std::unordered_map<std::string, std::string> node2contracted;
    for (auto&& node_pair : nodes_to_contract) {
        process_contracting(node_pair, node2contracted);
    }

    auto contract_self_loops = [&](std::string node, bool update_length = true)
        {
            if (graph_vis[node].outgoing_edges.find(node) != graph_vis[node].outgoing_edges.end() && graph_vis[node].outgoing_edges[node].size() >= 1 && graph_vis[node].number_of_contracted_edge >= 1) {
                long total_len = 0;
                std::vector<int> lens;

                std::vector<int> vector_to_remove;
                for (size_t i = 0; i < graph_vis[node].outgoing_edges[node].size();++i) {
                    if (graph_vis[node].outgoing_edges[node][i].sequence.size() >= 1000000)
                        continue;
                    vector_to_remove.push_back(int(i));
                    auto& e = graph_vis[node].outgoing_edges[node][i];
                    if (graph_vis[node].sequence.empty())
                        graph_vis[node].sequence += ("NNNNNNNNNNNNNNNNNNNN" + e.sequence + "NNNNNNNNNNNNNNNNNNNN");
                    else
                        graph_vis[node].sequence += (e.sequence + "NNNNNNNNNNNNNNNNNNNN");
                    total_len += e.length;
                    lens.push_back(e.length);
                    graph_vis[node].circles += 1;
                }
                if (!vector_to_remove.empty()) {
                    if (update_length) {
                        std::sort(lens.begin(), lens.end());
                        size_t size = lens.size();
                        if (size % 2 == 0) {
                            graph_vis[node].median_length = (lens[size / 2 - 1] + lens[size / 2]) / 2.0;
                        }
                        else {
                            graph_vis[node].median_length = lens[size / 2];
                        }

                        graph_vis[node].total_length = total_len;
                    }
                    else {
                        remove_items_from_vector(graph_vis[node].outgoing_edges[node], vector_to_remove);

                        vector_to_remove.clear();
                        for (size_t i = 0; i < graph_vis[node].incoming_edges[node].size();++i) {
                            if (graph_vis[node].incoming_edges[node][i].sequence.size() >= 1000000)
                                continue;
                            vector_to_remove.push_back(int(i));
                        }
                        remove_items_from_vector(graph_vis[node].incoming_edges[node], vector_to_remove);

                        if (graph_vis[node].outgoing_edges[node].empty()) {
                            graph_vis[node].outgoing_edges.erase(node);
                            graph_vis[node].incoming_edges.erase(node);
                        }
                    }
                }
            }
        };

    // remove self-loops, add statistics to the contracted node
    for (auto&& node : graph_vis) {
        contract_self_loops(node.first);
    }

    // contract edges of similar length to contracted nodes
    int contracted_edges = 1;
    while (contracted_edges) {
        contracted_edges = 0;
        std::vector<Nodes_To_Contract> nodes_to_contract_similar;
        node1_to_node2_scanned.clear();
        for (auto&& node : graph_vis) {
            //skip non-contracted nodes and contracted nodes with no circles
            if (node.second.number_of_contracted_edge == 0 || node.second.median_length == 0)
                continue;
            for (auto&& n_in : graph_vis[node.first].incoming_edges) {
                if (node.first == reverse_complementary_node(n_in.first) || node.first == n_in.first)
                    continue;
                if (node1_to_node2_scanned[n_in.first].find(node.first) != node1_to_node2_scanned[node.first].end() || node1_to_node2_scanned[reverse_complementary_node(node.first)].find(reverse_complementary_node(n_in.first)) != node1_to_node2_scanned[reverse_complementary_node(node.first)].end())
                    continue;
                for (auto&& e : n_in.second) {
                    if (e.length * 0.8 < node.second.median_length && e.length < 1000000) {
                        node1_to_node2_scanned[n_in.first].insert(node.first);
                        node1_to_node2_scanned[reverse_complementary_node(node.first)].insert(reverse_complementary_node(n_in.first));
                        nodes_to_contract_similar.emplace_back(Nodes_To_Contract(n_in.first, node.first, e.length, e.sequence));
                        nodes_to_contract_similar.emplace_back(Nodes_To_Contract(reverse_complementary_node(node.first), reverse_complementary_node(n_in.first), e.length, reverse_complementary(e.sequence)));
                        contracted_edges += 1;
                    }
                }
            }

            std::vector<std::string> n_outs;
            for (auto&& n_out : graph_vis[node.first].outgoing_edges) {
                if (node.first == reverse_complementary_node(n_out.first) || node.first == n_out.first)
                    continue;
                if (node1_to_node2_scanned[node.first].find(n_out.first) != node1_to_node2_scanned[node.first].end() || node1_to_node2_scanned[reverse_complementary_node(n_out.first)].find(reverse_complementary_node(node.first)) != node1_to_node2_scanned[reverse_complementary_node(n_out.first)].end())
                    continue;
                for (auto&& e : n_out.second) {
                    if (e.length * 0.8 < node.second.median_length && e.length < 1000000) {
                        node1_to_node2_scanned[node.first].insert(n_out.first);
                        node1_to_node2_scanned[reverse_complementary_node(n_out.first)].insert(reverse_complementary_node(node.first));
                        nodes_to_contract_similar.emplace_back(Nodes_To_Contract(node.first, n_out.first, e.length, e.sequence));
                        nodes_to_contract_similar.emplace_back(Nodes_To_Contract(reverse_complementary_node(n_out.first), reverse_complementary_node(node.first), e.length, reverse_complementary(e.sequence)));
                        contracted_edges += 1;
                    }
                }
            }
        }

        for (auto&& node_pair : nodes_to_contract_similar) {
            process_contracting(node_pair, node2contracted, true);
        }

        // remove self-loops
        for (auto&& node : graph_vis) {
            contract_self_loops(node.first, true);
        }
    }

    for (auto&& node : graph_vis) {
        contract_self_loops(node.first, false);
    }

    if (simplify) {
        // update sequences related to contracted node
        for (auto&& node : graph_vis) {
            //skip non-contracted nodes and contracted nodes with no circles
            if (node.second.number_of_contracted_edge == 0)
                continue;

            bool flag = false;
            if (node.second.incoming_edges.size() == 1 && node.second.outgoing_edges.size() == 1) {
                std::string in_node, out_node;
                for (auto&& e : node.second.incoming_edges)
                    in_node = e.first;
                for (auto&& e : node.second.outgoing_edges)
                    out_node = e.first;

                std::vector<std::string>& in_path = graph_vis[in_node].outgoing_edges.at(node.first).at(0).path_nodes_in_original_graph;
                std::vector<std::string>& out_path = graph_vis[node.first].outgoing_edges.at(out_node).at(0).path_nodes_in_original_graph;

                if (in_path.at(in_path.size() - 1) == out_path.at(0)) {
                    graph_vis[node.first].sequence = graph.at(in_path.at(in_path.size() - 1)).sequence;
                    flag = true;
                    for (auto&& e : node.second.incoming_edges) {
                        for (auto&& e_in : e.second) {
                            assert(e_in.length == e_in.sequence.size());
                            assert(e_in.sequence.substr(e_in.sequence.size() - graph_vis[node.first].sequence.size()) == graph_vis[node.first].sequence);
                        }
                        for (auto&& e_out : graph_vis[e.first].outgoing_edges[node.first]) {
                            assert(e_out.length == e_out.sequence.size());
                            assert(e_out.sequence.substr(e_out.sequence.size() - graph_vis[node.first].sequence.size()) == graph_vis[node.first].sequence);
                        }
                    }
                    for (auto&& e : node.second.outgoing_edges) {
                        for (auto&& e_out : e.second) {
                            assert(e_out.length == e_out.sequence.size());
                            assert(e_out.sequence.substr(0, graph_vis[node.first].sequence.size()) == graph_vis[node.first].sequence);
                        }
                        for (auto&& e_in : graph_vis[e.first].incoming_edges[node.first]) {
                            assert(e_in.length == e_in.sequence.size());
                            assert(e_in.sequence.substr(0, graph_vis[node.first].sequence.size()) == graph_vis[node.first].sequence);
                        }
                    }
                    std::cout << "Modify contracted node sequence for " << node.first << ", new length " << graph_vis[node.first].sequence.size() << std::endl;
                }
            }

            if (node.second.sequence.empty()) {
                node.second.sequence = "NNNNNNNNNNNNNNNNNNNN";
            }

            // add sequences of contracted node to edges for consistency, skip continuous path of incoming and outgoing edges
            if (!flag) {
                for (auto&& e : node.second.incoming_edges) {
                    for (auto&& e_in : e.second) {
                        e_in.sequence = e_in.sequence + node.second.sequence;
                        e_in.length += graph_vis[node.first].sequence.size();
                        assert(e_in.sequence.size() == e_in.length);
                    }
                    for (auto&& e_out : graph_vis[e.first].outgoing_edges[node.first]) {
                        e_out.sequence = e_out.sequence + node.second.sequence;
                        e_out.length += graph_vis[node.first].sequence.size();
                        assert(e_out.sequence.length() == e_out.length);
                    }
                }
                for (auto&& e : node.second.outgoing_edges) {
                    for (auto&& e_out : e.second) {
                        e_out.sequence = node.second.sequence + e_out.sequence;
                        e_out.length += graph_vis[node.first].sequence.size();
                        assert(e_out.sequence.length() == e_out.length);
                    }
                    for (auto&& e_in : graph_vis[e.first].incoming_edges[node.first]) {
                        e_in.sequence = node.second.sequence + e_in.sequence;
                        e_in.length += graph_vis[node.first].sequence.size();
                        assert(e_in.sequence.size() == e_in.length);
                    }
                }
            }

        }

        // auto graph_cp = graph;
        this->graph = graph_vis;
        // this->graph = graph_cp;
        unsigned removed_bulges = 1;
        while (removed_bulges) {
            merge_non_branching_paths(true);
            multi_bulge_removal(removed_bulges);
        }
        this->write_graph(prefix, 1000000, true);
    }
    else {
        auto graph_cp = graph;
        this->graph = graph_vis;
        this->write_graph(prefix, 1000000, true);
        this->graph = graph_cp;
    }
    return;
}

void Graph::write_graph_colored(const std::string& prefix, const std::string& genomes) {
    // read reference sequences from fasta file
    struct FastaEntry {
        std::string header;
        std::string sequence;
        void clear() { header.clear();  sequence.clear(); }
    };
    std::ifstream genome_seq(genomes);
    std::vector<FastaEntry> fastaEntries;
    std::string line;
    FastaEntry currentEntry;
    while (std::getline(genome_seq, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }
        if (line[0] == '>') {
            // This is a header line
            if (!currentEntry.header.empty() || !currentEntry.sequence.empty()) {
                if (currentEntry.sequence.size() > 1000)
                    fastaEntries.push_back(currentEntry);
                currentEntry.clear(); // Reset for the next entry
            }
            currentEntry.header = line.substr(1);
        }
        else {
            currentEntry.sequence += line;
        }
    }
    if (!currentEntry.header.empty() || !currentEntry.sequence.empty()) {
        if (currentEntry.sequence.size() > 1000)
            fastaEntries.push_back(currentEntry);
    }
    std::cout << fastaEntries.size() << std::endl;
    genome_seq.close();

    // for each edge, find a reference if aligned
    int cnt = 0;
    int cnt_aligned = 0, cnt_unaligned = 0;
    for (auto&& node : graph) {
        cnt += 1;
        if (cnt % 100 == 0)
            std::cout << "Aligned " << cnt << " nodes, " << cnt_aligned << " edges aligned" << std::endl;
        for (auto&& edge : node.second.outgoing_edges) {
            for (auto&& e : edge.second) {
#pragma omp parallel for
                for (int i = 0; i < fastaEntries.size(); ++i) {
                    // double PI = 1.0 * matches_by_edlib(e.sequence, fastaEntries[i].sequence) / e.sequence.size();
                    // std::cout << PI << std::endl;
                    // if (PI > 0.9)
                    //     e.ref_ids.push_back(fastaEntries[i].header);
                    if (fastaEntries[i].sequence.find(e.sequence) != std::string::npos || fastaEntries[i].sequence.find(reverse_complementary(e.sequence)) != std::string::npos) {
#pragma omp critical
                        e.ref_ids.push_back(fastaEntries[i].header);
                    }
                }
                if (!e.ref_ids.empty()) cnt_aligned++;
                else cnt_unaligned++;
            }
        }
        for (auto&& edge : node.second.incoming_edges) {
            for (auto&& e : edge.second) {
#pragma omp parallel for
                for (int i = 0; i < fastaEntries.size(); ++i) {
                    // double PI = 1.0 * matches_by_edlib(e.sequence, fastaEntries[i].sequence) / e.sequence.size();
                    // std::cout << PI << std::endl;
                    // if (PI > 0.9)
                    //     e.ref_ids.push_back(fastaEntries[i].header);
                    if (fastaEntries[i].sequence.find(e.sequence) != std::string::npos || fastaEntries[i].sequence.find(reverse_complementary(e.sequence)) != std::string::npos) {
#pragma omp critical
                        e.ref_ids.push_back(fastaEntries[i].header);
                    }
                }
                if (!e.ref_ids.empty()) cnt_aligned++;
                else cnt_unaligned++;
            }
        }
    }
    write_graph(prefix, 1000000, false, true);
}

void Graph::write_graph_colored_from_bam(const std::string& prefix, const std::string& bam_processed) {
    std::unordered_map<std::string, std::vector<std::string>> edge2refs;
    std::ifstream ref_mapping(bam_processed);
    std::string line;

    std::unordered_map<std::string, std::vector< std::string>> query_to_refs;
    while (std::getline(ref_mapping, line)) {
        std::istringstream iss(line);
        std::string query_name, ref_name, aligned_len, identity;
        iss >> query_name >> ref_name >> aligned_len >> identity;
        if (iss.fail()) continue;

        query_to_refs[query_name].push_back(ref_name + " " + aligned_len + " " + identity);
        // std::cout << "Read alignment for " << query_name << ": " << ref_name + " " + aligned_len + " " + identity << std::endl;
    }

    int cnt = 0, cnt_aligned = 0, cnt_unaligned = 0;
    for (auto&& node : graph) {
        cnt += 1;
        // if (cnt % 100 == 0)
        //     std::cout << "Processed " << cnt << " nodes, " << cnt_aligned << " edges aligned, " << cnt_unaligned << " edges unaligned." << std::endl;
        for (auto&& edge : node.second.outgoing_edges) {
            for (auto&& e : edge.second) {
                if (query_to_refs.find(e.label) != query_to_refs.end()) {
                    // std::cout << "Find ref for " << e.label << std::endl;
                    e.ref_ids = query_to_refs[e.label];
                }
                if (!e.ref_ids.empty()) cnt_aligned++;
                else cnt_unaligned++;
            }
        }
        for (auto&& edge : node.second.incoming_edges) {
            for (auto&& e : edge.second) {
                if (query_to_refs.find(e.label) != query_to_refs.end())
                    e.ref_ids = query_to_refs[e.label];
                if (!e.ref_ids.empty()) cnt_aligned++;
                else cnt_unaligned++;
            }
        }
    }
    // std::cout << "Processed " << cnt << " nodes, " << cnt_aligned << " edges aligned, " << cnt_unaligned << " edges unaligned." << std::endl;
    write_graph(prefix, 1000000, false, true);
}

void Graph::write_graph_gfa(const std::string& prefix) {
    std::string graph_gfa = prefix + ".gfa";
    std::cout << "Write graph " << graph_gfa << std::endl;
    std::ofstream file_gfa(graph_gfa);

    std::unordered_map<std::string, std::pair<std::string, char>> edgeid2gfa_node;
    file_gfa << "H\tVN:Z:1.0" << std::endl;
    int cnt = 10000;
    std::unordered_set<std::string> traversed_labels;
    for (auto&& n1 : graph) {
        for (auto&& n2 : n1.second.outgoing_edges) {
            for (auto&& e : n2.second) {
                if (traversed_labels.find(e.label) == traversed_labels.end()) {
                    traversed_labels.insert(e.label);
                    traversed_labels.insert(e.rc_label);std::string new_label = std::to_string(cnt++);
                    edgeid2gfa_node[e.label] = std::make_pair(new_label, '+');
                    edgeid2gfa_node[e.rc_label] = std::make_pair(new_label, '-');
                    if (n1.first != n2.first) {
                        // edgeid2gfa_node[e.label] = std::make_pair(e.label + "_" + e.rc_label, '+');
                        // edgeid2gfa_node[e.rc_label] = std::make_pair(e.label + "_" + e.rc_label, '-');
                        // file_gfa << "S\t" << e.label + "_" + e.rc_label << "\t" << e.sequence << '\n';
                        file_gfa << "S\t" << new_label << "\t" << e.sequence << '\n';
                        assert(e.sequence.size() > 0);
                    }
                    else {
                        if (e.sequence.size() >= n1.second.sequence.size() * 2)
                            file_gfa << "S\t" << new_label << "\t" << e.sequence.substr(0, e.sequence.size() - n1.second.sequence.size()) << '\n';
                        else
                            file_gfa << "S\t" << new_label << "\t" << e.sequence << '\n';
                        assert(e.sequence.size() > 0);
                    }
                }
            }
        }
        if (n1.second.outgoing_edges.empty() && n1.second.incoming_edges.empty()) {
            if (traversed_labels.find(n1.first) == traversed_labels.end()) {
                traversed_labels.insert(n1.first);
                traversed_labels.insert(reverse_complementary_node(n1.first));
                std::string new_label = std::to_string(cnt++);
                file_gfa << "S\t" << new_label << "\t" << n1.second.sequence << '\n';
            }
        }
    }

    // std::unordered_set<std::string> traversed_nodes;
    // for (auto&& n : graph) {
    //     if (traversed_nodes.find(n.first) != traversed_nodes.end())
    //         continue;
    //     traversed_nodes.insert(n.first);
    //     traversed_nodes.insert(reverse_complementary_node(n.first));
    //     for (auto&& n_o : n.second.outgoing_edges) {
    //         for (auto&& e_o : n_o.second) {
    //             for (auto&& n_i : n.second.incoming_edges) {
    //                 for (auto&& e_i : n_i.second) {
    //                     assert(edgeid2gfa_node.find(e_i.label) != edgeid2gfa_node.end());
    //                     assert(edgeid2gfa_node.find(e_o.label) != edgeid2gfa_node.end());
    //                     file_gfa << "L\t" << edgeid2gfa_node[e_i.label].first << "\t" << edgeid2gfa_node[e_i.label].second << "\t" << edgeid2gfa_node[e_o.label].first << "\t" << edgeid2gfa_node[e_o.label].second << "\t" << n.second.sequence.size() << "M\n";
    //                 }
    //             }
    //         }
    //     }
    // }

    file_gfa.close();
}