#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "dot_graph.hpp"
#include <cassert>
#include <algorithm>
#include <cmath>
#include "edlib.h"
#include <sstream>
#include <iomanip>
#include <sys/stat.h>

Edge::Edge() {}

Edge::Edge(char& base, const unsigned& length, const std::string& sequence, const double& multi) {
    this->start_base = base;
    this->length = length;
    this->sequence = sequence;
    this->multiplicity = multi;
}

void Edge::add_multi_from_edge_or_path(Edge& edge_or_path) {
    this->multiplicity += edge_or_path.multiplicity;
}

void Edge::add_multi_from_edge_or_path(Path& edge_or_path) {
    this->multiplicity += edge_or_path.min_multi;
}

void Edge::remove_multi_from_path(Path& edge_or_path) {
    if (this->multiplicity >= edge_or_path.min_multi)
        this->multiplicity -= edge_or_path.min_multi;
    else
        this->multiplicity = 0;
}

Node::Node() {}

Path::Path() {}

void Path::update_min_multi(Edge edge) {
    if (this->min_multi == 0)
        this->min_multi = edge.multiplicity;
    else
        this->min_multi = std::min(this->min_multi, edge.multiplicity);
}

Bulge::Bulge() {}

Bulge::Bulge(Path p1, Path p2, double identity) {
    this->leg1 = p1;
    this->leg2 = p2;
    this->max_identity = identity;
    this->num_edges = p1.nodes.size() + p2.nodes.size();
}

void Bulge::check_conflict(Bulge& bulge) {
    std::vector<Path> path_resolutions_1, path_resolutions_2;
    path_resolutions_1.insert(path_resolutions_1.end(), this->path_resolutions1.begin(), this->path_resolutions1.end());
    path_resolutions_1.insert(path_resolutions_1.end(), this->path_resolutions2.begin(), this->path_resolutions2.end());
    path_resolutions_2.insert(path_resolutions_2.end(), bulge.path_resolutions1.begin(), bulge.path_resolutions1.end());
    path_resolutions_2.insert(path_resolutions_2.end(), bulge.path_resolutions2.begin(), bulge.path_resolutions2.end());
    for (auto&& i : path_resolutions_1) {
        for (auto&& j : path_resolutions_2) {
            if (i.nodes.size() != j.nodes.size())
                continue;

            bool flag = false;
            for (int k = 1; k < i.nodes.size() - 1; ++k) {
                if (i.nodes[k] != j.nodes[k]) {
                    flag = true;
                    break;
                }
            }
            if (flag)
                continue;

            if (i.nodes[0] == j.nodes[0] && i.nodes[i.nodes.size() - 1] != j.nodes[j.nodes.size() - 1]) {
                this->is_confict = true;
                bulge.is_confict = true;
            }
            else if (i.nodes[0] != j.nodes[0] && i.nodes[i.nodes.size() - 1] == j.nodes[j.nodes.size() - 1]) {
                this->is_confict = true;
                bulge.is_confict = true;
            }
        }
    }
}

bool Graph::check_2_in_2_out(Path& leg1) {
    for (int i = 1; i < leg1.nodes.size() - 1; ++i) {
        Path path_resolution;
        if (this->graph[leg1.nodes[i]].incoming_edges.size() > 1) {
            if (this->graph[leg1.nodes[i]].incoming_edges.size() == 2 && this->graph[leg1.nodes[i]].incoming_edges.find(leg1.nodes[i]) != this->graph[leg1.nodes[i]].incoming_edges.end())
                continue;
            int j = i;
            this->add_node_to_path(path_resolution, leg1.nodes[i - 1]);
            while (j < leg1.nodes.size() - 1) {
                bool flag = false;
                if (graph[leg1.nodes[j]].outgoing_edges.size() > 2)
                    flag = true;
                if (graph[leg1.nodes[j]].outgoing_edges.size() == 2 && graph[leg1.nodes[j]].outgoing_edges.find(leg1.nodes[j]) == graph[leg1.nodes[j]].outgoing_edges.end())
                    flag = true;
                // find a 2-in-2-out path
                if (flag) {
                    Path path_2_2 = path_resolution;
                    this->add_node_to_path(path_2_2, leg1.nodes[j], leg1.bulge_legs[j - 1]);
                    this->add_node_to_path(path_2_2, leg1.nodes[j + 1], leg1.bulge_legs[j]);
                    return true;
                }

                this->add_node_to_path(path_resolution, leg1.nodes[j], leg1.bulge_legs[j - 1]);
                ++j;
            }
        }
    }
    return false;
}

void Graph::find_2_in_2_out(Bulge& bulge) {
    Path& leg1 = bulge.leg1, leg2 = bulge.leg2;
    for (int i = 1; i < leg1.nodes.size() - 1; ++i) {
        Path path_resolution;
        if (this->graph[leg1.nodes[i]].incoming_edges.size() > 1) {
            if (this->graph[leg1.nodes[i]].incoming_edges.size() == 2 && this->graph[leg1.nodes[i]].incoming_edges.find(leg1.nodes[i]) != this->graph[leg1.nodes[i]].incoming_edges.end())
                continue;
            int j = i;
            this->add_node_to_path(path_resolution, leg1.nodes[i - 1]);
            while (j < leg1.nodes.size() - 1) {
                bool flag = false;
                if (graph[leg1.nodes[j]].outgoing_edges.size() > 2)
                    flag = true;
                if (graph[leg1.nodes[j]].outgoing_edges.size() == 2 && graph[leg1.nodes[j]].outgoing_edges.find(leg1.nodes[j]) == graph[leg1.nodes[j]].outgoing_edges.end())
                    flag = true;
                // find a 2-in-2-out path
                if (flag) {
                    Path path_2_2 = path_resolution;
                    this->add_node_to_path(path_2_2, leg1.nodes[j], leg1.bulge_legs[j - 1]);
                    this->add_node_to_path(path_2_2, leg1.nodes[j + 1], leg1.bulge_legs[j]);
                    // std::cout << "2-in-2-out in bulge leg " << leg1.nodes[0];
                    // for (int k = 1; k < leg1.nodes.size();++k)
                    //     std::cout << "->" << leg1.nodes[k];
                    // std::cout << ": " << path_2_2.nodes[0];
                    // for (int k = 1; k < path_2_2.nodes.size();++k)
                    //     std::cout << "->" << path_2_2.nodes[k];
                    // std::cout << std::endl;
                    bulge.path_resolutions1.emplace_back(path_2_2);
                }

                this->add_node_to_path(path_resolution, leg1.nodes[j], leg1.bulge_legs[j - 1]);
                ++j;
            }
        }
    }

    for (int i = 1; i < leg2.nodes.size() - 1; ++i) {
        Path path_resolution;
        if (this->graph[leg2.nodes[i]].incoming_edges.size() > 1) {
            if (this->graph[leg2.nodes[i]].incoming_edges.size() == 2 && this->graph[leg2.nodes[i]].incoming_edges.find(leg2.nodes[i]) != this->graph[leg2.nodes[i]].incoming_edges.end())
                continue;
            int j = i;
            this->add_node_to_path(path_resolution, leg2.nodes[i - 1]);
            while (j < leg2.nodes.size() - 1) {
                bool flag = false;
                if (graph[leg2.nodes[j]].outgoing_edges.size() > 2)
                    flag = true;
                if (graph[leg2.nodes[j]].outgoing_edges.size() == 2 && graph[leg2.nodes[j]].outgoing_edges.find(leg2.nodes[j]) == graph[leg2.nodes[j]].outgoing_edges.end())
                    flag = true;
                // find a 2-in-2-out path
                if (flag) {
                    Path path_2_2 = path_resolution;
                    this->add_node_to_path(path_2_2, leg2.nodes[j], leg2.bulge_legs[j - 1]);
                    this->add_node_to_path(path_2_2, leg2.nodes[j + 1], leg2.bulge_legs[j]);
                    // std::cout << "2-in-2-out in bulge leg " << leg2.nodes[0];
                    // for (int k = 1; k < leg2.nodes.size();++k)
                    //     std::cout << "->" << leg2.nodes[k];
                    // std::cout << ": " << path_2_2.nodes[0];
                    // for (int k = 1; k < path_2_2.nodes.size();++k)
                    //     std::cout << "->" << path_2_2.nodes[k];
                    // std::cout << std::endl;
                    bulge.path_resolutions2.emplace_back(path_2_2);
                }

                this->add_node_to_path(path_resolution, leg2.nodes[j], leg2.bulge_legs[j - 1]);
                ++j;
            }
        }
    }
}

Genome::Genome() {}

void Genome::add_node_to_end(std::string node) {
    this->genome_path.push_back(node);
}

bool Genome::find_path(Path& path, std::vector<int>& pos) {
    pos.clear();
    if (path.nodes.empty())
        return false;
    std::string node = path.nodes.at(0);
    auto iter = std::find(this->genome_path.begin(), this->genome_path.end(), node);
    while (iter != this->genome_path.end()) {
        bool flag = true;
        for (int i = 1; i < path.nodes.size(); ++i) {
            if (iter + i == this->genome_path.end()) {
                flag = false;
                break;
            }
            if (path.nodes.at(i) != *(iter + i)) {
                flag = false;
                break;
            }
        }
        if (flag) {
            pos.push_back(iter - this->genome_path.begin());
            iter = std::find(iter + path.nodes.size() - 1, this->genome_path.end(), node);
        }
        else
            iter = std::find(iter + 1, this->genome_path.end(), node);
    }
    if (pos.empty())
        return false;
    else
        return true;
}

void Genome::construct_edge_orders() {
    this->node2order.clear();
    for (int i = 0; i + 1 < this->genome_path.size(); ++i) {
        this->node2order[this->genome_path.at(i)][this->genome_path.at(i + 1)].push_back(i);
    }
}

void Genome::erase_node(std::string node) {
    auto iter = std::find(this->genome_path.begin(), this->genome_path.end(), node);
    while (iter != this->genome_path.end()) {
        this->genome_path.erase(iter);
        iter = std::find(this->genome_path.begin(), this->genome_path.end(), node);
    }
}

Graph::Graph() {}

int Graph::get_num_nodes() {
    return this->graph.size();
}

void Graph::read_graph(std::string& output, std::string& restart_from, std::string& graph_dot, const std::string& graph_fasta, const std::string& nodes_fasta) {
    if (mkdir(output.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1) {
        if (errno == EEXIST) {
            if (!restart_from.empty())
                read_from_dot(output + "/" + restart_from + ".dot", output + "/" + restart_from + ".fasta", output + "/" + restart_from + ".vertices.fasta");
            else {
                std::cout << "Output directory already exists. Please delete it or set restart_from." << std::endl;
                exit(0);
            }
        }
        else {
            std::cout << "Create output directory faliled." << std::endl;
            exit(0);
        }
    }

    if (graph.empty())
        read_from_dot(graph_dot, graph_fasta, nodes_fasta);
}

// Read graph from DOT file
void Graph::read_from_dot(const std::string& graph_dot, const std::string& graph_fasta, const std::string& nodes_fasta) {
    // load fasta sequence
    std::string e_id;
    std::string line, node1, node2, length, multiplicity, start_base;
    std::unordered_map<std::string, std::string> edge2sequence, node2sequence;

    // read input fasta
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
            std::cout << "Read " << e_id << " from " << graph_fasta << std::endl;
        }
    }

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
            std::cout << "Read " << node << " from " << nodes_fasta << std::endl;
        }
    }

    // load graph
    long cnt_edge = 0;
    std::ifstream dot_file(graph_dot);
    while (getline(dot_file, line)) {
        // line is neithor a node nor an edge
        if (line.find('[') == std::string::npos)
            continue;
        // line is an edge, in this case, all the node should already be loaded
        if (line.find("->") != std::string::npos) {
            // extract node name
            size_t pos1 = line.find("->");
            std::string start_name = line.substr(1, pos1 - 1);
            size_t pos2 = line.find('[');
            std::string end_name = line.substr(pos1 + 2, pos2 - pos1 - 2);
            assert(this->graph.find(start_name) != this->graph.end() && this->graph.find(end_name) != this->graph.end());

            // parse edge label: starting base, length, multiplicity and sequence
            pos1 = line.find("size");
            std::string edge_label = line.substr(pos2 + 14, pos1 - pos2 - 16);

            std::cout << "Read edge " << start_name << " -> " << end_name << ", " << edge_label << " from " << graph_dot << std::endl;

            assert(edge2sequence.find(edge_label) != edge2sequence.end());
            int length = edge2sequence[edge_label].size();
            double multiplicity = 0;
            char start_base;
            Edge edge = Edge(start_base, length, edge2sequence[edge_label], multiplicity);
            assert(edge2sequence[edge_label].substr(0, graph[start_name].sequence.size()) == graph[start_name].sequence);
            edge.path_nodes_in_original_graph.push_back(start_name);
            edge.path_nodes_in_original_graph.push_back(end_name);
            edge.path_edges_in_original_graph.push_back(edge_label);
            edge.label = edge_label;

            this->graph[start_name].outgoing_edges[end_name].push_back(edge);
            this->graph[end_name].incoming_edges[start_name].push_back(edge);
            cnt_edge += 1;
        }
        // line is an node
        else {
            std::string node_name = line.substr(1, line.find('[') - 1);
            Node node;
            assert(node2sequence.find(node_name) != node2sequence.end());
            node.sequence = node2sequence[node_name];
            this->graph[node_name] = node;
        }
    }
    dot_file.close();
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

std::string Graph::doubleToString(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << value;
    return stream.str();
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
        if (contracted && node.second.number_of_contracted_edge > 0)
            file_dot << "\"" << node.first.substr(0, node.first.find_first_of('_')) + "_" + node.first.substr(node.first.find_last_of('_') + 1) + "_N" + std::to_string(node.second.number_of_contracted_edge) + "_L" + std::to_string(node.second.length_of_contracted_edge) << "\" [style=filled fillcolor=\"white\"]\n";
        else
            file_dot << node.first << " [style=filled fillcolor=\"white\"]\n";
        if (node.second.number_of_contracted_edge > max_contracted) {
            max_contracted_node = node.first.substr(0, node.first.find_first_of('_')) + "_" + node.first.substr(node.first.find_last_of('_') + 1) + "_N" + std::to_string(node.second.number_of_contracted_edge) + "_L" + std::to_string(node.second.length_of_contracted_edge);
            max_contracted = node.second.number_of_contracted_edge;
        }
    }
    if (max_contracted > 0)
        std::cout << "Max contracted node: " << max_contracted_node << " has " << max_contracted << " edges." << std::endl;

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
        {"Y", "#969696"}
    };

    std::unordered_set<std::string> edgelabels;
    for (auto&& node : this->graph) {
        for (auto&& edges : node.second.outgoing_edges) {
            for (auto&& edge : edges.second) {
                edgelabels.insert(edge.label);
            }
        }
    }
    for (auto&& node : this->graph) {
        std::string start_node = node.first;
        if (contracted && node.second.number_of_contracted_edge > 0)
            start_node = node.first.substr(0, node.first.find_first_of('_')) + "_" + node.first.substr(node.first.find_last_of('_') + 1) + "_N" + std::to_string(node.second.number_of_contracted_edge) + "_L" + std::to_string(node.second.length_of_contracted_edge);
        for (auto&& edges : node.second.outgoing_edges) {
            std::string end_node = edges.first;
            if (!nodes_retain.empty() && nodes_retain.find(node.first) == nodes_retain.end() && nodes_retain.find(end_node) == nodes_retain.end())
                continue;
            if (contracted && graph[edges.first].number_of_contracted_edge > 0)
                end_node = edges.first.substr(0, edges.first.find_first_of('_')) + '_' + edges.first.substr(edges.first.find_last_of('_') + 1) + "_N" + std::to_string(graph[edges.first].number_of_contracted_edge) + "_L" + std::to_string(graph[edges.first].length_of_contracted_edge);
            for (auto&& edge : edges.second) {
                num_edges += 1;
                if (edge.label.empty()) {
                    edge.label = get_unique_label(edgelabels);
                    edgelabels.insert(edge.label);
                }
                if (contracted || colored) {
                    if (edge.ref_ids.empty())
                        file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.length << "\" color=\"black\"]\n";
                    else {
                        file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.length;
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
                    file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.length << "\" color=\"black\"]\n";

                file_fasta << ">" << edge.label << "\n";
                file_fasta << edge.sequence << "\n";
                file_path << ">" << edge.label << " " << edge.length << "\n";
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

void Graph::write_graph_contracted(const std::string& prefix, int min_length) {
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
        for (auto&& edge : node.second.outgoing_edges) {
            // do not deal with palindromic bulges
            if (node.first == edge.first)
                continue;
            if (node1_to_node2_scanned[node.first].find(edge.first) != node1_to_node2_scanned[node.first].end())
                continue;
            for (int i = 0; i < edge.second.size();++i) {
                // the edge should be collapsed
                if (edge.second.at(i).length <= min_length) {
                    nodes_to_contract.emplace_back(Nodes_To_Contract(node.first, edge.first, edge.second.at(i).length, edge.second.at(i).sequence));
                }
            }
            node1_to_node2_scanned[node.first].insert(edge.first);
        }
    }

    std::unordered_map<std::string, std::string> node2contracted;
    for (auto&& node_pair : nodes_to_contract) {
        std::string node1 = node_pair.node1, node2 = node_pair.node2;
        while (node2contracted.find(node1) != node2contracted.end())
            node1 = node2contracted[node1];
        while (node2contracted.find(node2) != node2contracted.end())
            node2 = node2contracted[node2];

        if (node1 == node2) {
            double sum_multi = 0;
            for (auto&& n : graph[node_pair.node1].incoming_edges) {
                for (auto&& e : n.second)
                    sum_multi += e.multiplicity;
            }
            if (sum_multi > graph_vis[node1].max_in_multi) {
                graph_vis[node1].max_in_multi = sum_multi;
                graph_vis[node1].max_in_node = node_pair.node1;
            }

            sum_multi = 0;
            for (auto&& n : graph[node_pair.node2].incoming_edges) {
                for (auto&& e : n.second)
                    sum_multi += e.multiplicity;
            }
            if (sum_multi > graph_vis[node1].max_in_multi) {
                graph_vis[node1].max_in_multi = sum_multi;
                graph_vis[node1].max_in_node = node_pair.node2;
            }

            graph_vis[node1].number_of_contracted_edge += 1;
            graph_vis[node1].length_of_contracted_edge += node_pair.edge_length;
            // std::cout << "Contract " << node1 << "(" << node_pair.node1 << ")" << " and " << node2 << "(" << node_pair.node2 << ") to " << node1 << ", number " << graph_vis[node1].number_of_contracted_edge << ", length " << graph_vis[node1].length_of_contracted_edge << std::endl;
            int i = 0;
            for (;i < graph_vis[node2].incoming_edges[node1].size();++i) {
                if (graph_vis[node2].incoming_edges[node1][i].sequence == node_pair.edge_sequence)
                    break;
            }
            assert(i != graph_vis[node2].incoming_edges[node1].size());
            graph_vis[node2].incoming_edges[node1].erase(graph_vis[node2].incoming_edges[node1].begin() + i);
            if (graph_vis[node2].incoming_edges[node1].empty())
                graph_vis[node2].incoming_edges.erase(node1);

            int j = 0;
            for (;j < graph_vis[node1].outgoing_edges[node2].size();++j) {
                if (graph_vis[node1].outgoing_edges[node2][j].sequence == node_pair.edge_sequence)
                    break;
            }
            assert(j != graph_vis[node1].outgoing_edges[node2].size());
            graph_vis[node1].outgoing_edges[node2].erase(graph_vis[node1].outgoing_edges[node2].begin() + j);
            if (graph_vis[node1].outgoing_edges[node2].empty())
                graph_vis[node1].outgoing_edges.erase(node2);

            continue;
        }

        std::string node_contracted = node1 + "_" + node2;
        node2contracted[node1] = node_contracted;
        node2contracted[node2] = node_contracted;

        double sum_multi = 0;
        for (auto&& n : graph[node_pair.node1].incoming_edges) {
            for (auto&& e : n.second)
                sum_multi += e.multiplicity;
        }
        if (sum_multi > graph_vis[node_contracted].max_in_multi) {
            graph_vis[node_contracted].max_in_multi = sum_multi;
            graph_vis[node_contracted].max_in_node = node_pair.node1;
        }

        sum_multi = 0;
        for (auto&& n : graph[node_pair.node2].incoming_edges) {
            for (auto&& e : n.second)
                sum_multi += e.multiplicity;
        }
        if (sum_multi > graph_vis[node_contracted].max_in_multi) {
            graph_vis[node_contracted].max_in_multi = sum_multi;
            graph_vis[node_contracted].max_in_node = node_pair.node2;
        }
        graph_vis[node_contracted].number_of_contracted_edge = graph_vis[node1].number_of_contracted_edge + graph_vis[node2].number_of_contracted_edge + 1;
        graph_vis[node_contracted].length_of_contracted_edge = graph_vis[node1].length_of_contracted_edge + graph_vis[node2].length_of_contracted_edge + node_pair.edge_length;
        // std::cout << "Contract " << node1 << "(" << node_pair.node1 << ")" << " and " << node2 << "(" << node_pair.node2 << ") to " << node_contracted << ", number " << graph_vis[node_contracted].number_of_contracted_edge << ", length " << graph_vis[node_contracted].length_of_contracted_edge << std::endl;

        int i = 0;
        for (;i < graph_vis[node2].incoming_edges[node1].size();++i) {
            if (graph_vis[node2].incoming_edges[node1][i].sequence == node_pair.edge_sequence)
                break;
        }
        assert(i != graph_vis[node2].incoming_edges[node1].size());
        graph_vis[node2].incoming_edges[node1].erase(graph_vis[node2].incoming_edges[node1].begin() + i);
        if (graph_vis[node2].incoming_edges[node1].empty())
            graph_vis[node2].incoming_edges.erase(node1);

        int j = 0;
        for (;j < graph_vis[node1].outgoing_edges[node2].size();++j) {
            if (graph_vis[node1].outgoing_edges[node2][j].sequence == node_pair.edge_sequence)
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
    }
    auto graph_cp = graph;
    this->graph = graph_vis;
    this->write_graph(prefix, 1000000, true);
    this->graph = graph_cp;
    return;
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
    if (node.at(0) == '-')
        return node.substr(1);
    else
        return '-' + node;
}