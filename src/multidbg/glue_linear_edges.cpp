#include <string>
#include <iostream>
#include <cassert>
#include <queue>
#include <algorithm>
#include <set>
#include "dot_graph.hpp"
#include "edlib.h"
#include <cmath>
#include <filesystem>
#include <unistd.h>
#include <cstdlib>
#include <utils.hpp>
#include <fstream>

using namespace multidbg;
namespace fs = std::filesystem;

void Graph::merge_tips_L(unsigned& num_tips, std::unordered_map<std::string, std::vector<std::string>> nodes2bc) {
    num_tips = 0;
    std::set<std::string> nodes_to_remove;

    for (auto&& node : this->graph) {
        std::vector<std::string> outgoing_tips;
        for (auto&& n : graph[node.first].outgoing_edges) {
            if (graph[n.first].outgoing_edges.size() == 0) {
                assert(graph[reverse_complementary_node(n.first)].incoming_edges.size() == 0);
                if (n.second.size() == 1)
                    outgoing_tips.push_back(n.first);
            }
        }
        if (outgoing_tips.size() >= 2) {
            // move edges of all tips to the first tip
            int max_index_bc = -1;
            long max_length_bc = 0;
            int max_index_all = -1;
            long max_length_all = 0;
            for (int i = 0; i < outgoing_tips.size();++i) {
                if (nodes2bc.find(graph[outgoing_tips[i]].sequence) != nodes2bc.end() && graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length > max_length_bc) {
                    max_length_bc = graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length;
                    max_index_bc = i;
                }
                if (graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length > max_length_all) {
                    max_length_all = graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length;
                    max_index_all = i;
                }
            }
            int max_index = -1;
            if (max_index_bc != -1)
                max_index = max_index_bc;
            else
                max_index = max_index_all;

            for (int i = 0; i < outgoing_tips.size();++i) {
                if (i == max_index)
                    continue;
                if (nodes2bc.find(graph[outgoing_tips[i]].sequence) != nodes2bc.end())
                    continue;
                int prefix_size = 100000;
                std::string prefix_tip_target = graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).sequence.substr(0, prefix_size);
                std::string prefix_tip_to_merge = graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).sequence.substr(0, prefix_size);
                double sim = matches_by_edlib(prefix_tip_target, prefix_tip_to_merge);
                std::cout << "[MergeTip] Check tip " << outgoing_tips.at(i) << " length " << graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length << " to tip " << outgoing_tips.at(max_index) << " length " << graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length << " sim " << sim << std::endl;

                if (sim < 0.2)
                    continue;
                // if (sim < 0.9) {
                //     double sim_len = 1.0 * std::min(graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length, graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length)
                //         / std::max(graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length, graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length);
                //     if (sim_len < 0.8)
                //         continue;
                // }
                std::cout << "[MergeTip] Merge tip " << outgoing_tips.at(i) << " length " << graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length << " to tip " << outgoing_tips.at(max_index) << " length " << graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length << " sim " << sim << std::endl;
                merge_vecs(graph[node.first].outgoing_edges[outgoing_tips[max_index]], graph[node.first].outgoing_edges[outgoing_tips[i]]);
                merge_vecs(graph[outgoing_tips[max_index]].incoming_edges[node.first], graph[outgoing_tips[i]].incoming_edges[node.first]);
                merge_vecs(graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(outgoing_tips[max_index])], graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(outgoing_tips[i])]);
                merge_vecs(graph[reverse_complementary_node(outgoing_tips[max_index])].outgoing_edges[reverse_complementary_node(node.first)], graph[reverse_complementary_node(outgoing_tips[i])].outgoing_edges[reverse_complementary_node(node.first)]);
                graph[node.first].outgoing_edges.erase(outgoing_tips[i]);
                graph[outgoing_tips[i]].incoming_edges.erase(node.first);
                graph[reverse_complementary_node(node.first)].incoming_edges.erase(reverse_complementary_node(outgoing_tips[i]));
                graph[reverse_complementary_node(outgoing_tips[i])].outgoing_edges.erase(reverse_complementary_node(node.first));
                if (graph[outgoing_tips[i]].incoming_edges.size() == 0) {
                    assert(graph[reverse_complementary_node(outgoing_tips[i])].outgoing_edges.size() == 0);
                    nodes_to_remove.insert(outgoing_tips[i]);
                    nodes_to_remove.insert(reverse_complementary_node(outgoing_tips[i]));
                }
                num_tips++;
            }
        }
    }
    for (auto&& n : nodes_to_remove) {
        this->graph.erase(n);
    }
    unsigned bulges = 1;
    while (bulges) {
        multi_bulge_removal(bulges);
    }

    merge_non_branching_paths(true);
}

void Graph::merge_tips_into_edges_L(unsigned& num_tips, double ratio, bool only_tips, bool merge_long_tips, std::unordered_map<std::string, std::vector<std::string>> nodes2bc) {
    if (!only_tips) {
        unsigned removed_paths = 1;
        while (removed_paths) {
            resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.6, true, 2);
        }
        merge_non_branching_paths(true);
    }

    num_tips = 0;
    // traverse all nodes
    std::set<std::string> nodes_to_remove;
    for (auto&& node : graph) {
        if (nodes_to_remove.find(node.first) != nodes_to_remove.end())
            continue;
        std::vector<std::string> non_tips, tips;
        for (auto&& node_out : node.second.outgoing_edges) {
            if (graph[node_out.first].outgoing_edges.empty() && graph[node_out.first].incoming_edges.size() == 1) {
                if (node_out.second.size() == 1 && nodes2bc.find(graph[node_out.first].sequence) == nodes2bc.end())
                    tips.push_back(node_out.first);
            }
            else {
                if (node_out.second.size() == 1)
                    non_tips.push_back(node_out.first);
            }
        }
        if (tips.empty())
            continue;

        // check whether the prefix of the tip is similar to the prefix to the other edges; check whether the tip is short
        for (auto&& t : tips) {
            double max_sim = 0;
            std::string max_edge;
            // std::cout << "Check tip " << node.first << "->" << t << " multi " << graph[node.first].outgoing_edges[t].at(0).multiplicity << std::endl;
            for (auto&& e : non_tips) {
                // the tip length * ratio should be shorter than the edge length
                if (graph[node.first].outgoing_edges[t].at(0).length * ratio > graph[node.first].outgoing_edges[e].at(0).length)
                    continue;
                if (!merge_long_tips) {
                    // the tip multi * ratio should be lower than the edge multi
                    if (graph[node.first].outgoing_edges[t].at(0).multiplicity * ratio > graph[node.first].outgoing_edges[e].at(0).multiplicity)
                        continue;
                }

                if (merge_long_tips) {
                    // long tips should be merged only if the tip length * 0.8 is shorter than the edge length
                    if (1.0 * graph[node.first].outgoing_edges[e].at(0).length / graph[node.first].outgoing_edges[t].at(0).length < 0.8 && graph[node.first].outgoing_edges[t].at(0).length >= 1000000)
                        continue;
                }
                else {
                    // long tips should be merged only at the final stage
                    if (graph[node.first].outgoing_edges[t].at(0).length >= 1000000)
                        continue;
                }

                // should be very careful merging extra-long tips to edges, because we only mapps the 1Mb prefix
                if (graph[node.first].outgoing_edges[t].at(0).length >= 20000000)
                    continue;

                // calculate similarity
                int prefix_len = 1000000;
                std::string prefix_tip = graph[node.first].outgoing_edges[t].at(0).sequence.substr(0, prefix_len);
                std::string prefix_edge = graph[node.first].outgoing_edges[e].at(0).sequence.substr(0, prefix_len);
                double sim = matches_by_edlib(prefix_tip, prefix_edge);
                if (sim > max_sim) {
                    max_sim = sim;
                    max_edge = e;
                }

                // skip if node.first is in a palindromic bulge - dangerous
                bool flag = false;
                std::string curr_n = node.first;
                std::unordered_set<std::string> traversed_nodes;
                while (true) {
                    traversed_nodes.insert(curr_n);
                    if (graph[curr_n].incoming_edges.size() != 1) {
                        if (curr_n == node.first)
                            break;
                        if (curr_n != reverse_complementary_node(node.first))
                            break;
                        else {
                            flag = true;
                            break;
                        }
                    }
                    std::string n_incoming;
                    for (auto&& n : graph[curr_n].incoming_edges) {
                        n_incoming = n.first;
                    }
                    if (traversed_nodes.find(n_incoming) != traversed_nodes.end())
                        break;
                    curr_n = n_incoming;
                }
                if (flag)
                    continue;
            }
            if (max_edge.empty())
                continue;

            std::cout << "[RepairTip] Check tip " << node.first << "->" << t << " multi " << graph[node.first].outgoing_edges[t].at(0).multiplicity << " len " << graph[node.first].outgoing_edges[t].at(0).length << " to edge " << node.first << "->" << max_edge << " multi " << graph[node.first].outgoing_edges[max_edge].at(0).multiplicity << " len " << graph[node.first].outgoing_edges[max_edge].at(0).length << " with sim " << max_sim << std::endl;

            if (max_sim < 0.8)
                continue;

            graph[node.first].outgoing_edges[max_edge].at(0).multiplicity += 1.0 * graph[node.first].outgoing_edges[t].at(0).multiplicity * graph[node.first].outgoing_edges[t].at(0).length / graph[node.first].outgoing_edges[max_edge].at(0).length;
            graph[max_edge].incoming_edges[node.first].at(0).multiplicity = graph[node.first].outgoing_edges[max_edge].at(0).multiplicity;

            graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).multiplicity +=
                1.0 * graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(t)].at(0).multiplicity *
                graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(t)].at(0).length /
                graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).length;
            graph[reverse_complementary_node(max_edge)].outgoing_edges[reverse_complementary_node(node.first)].at(0).multiplicity = graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).multiplicity;

            nodes_to_remove.insert(t);
            nodes_to_remove.insert(reverse_complementary_node(t));
            num_tips += 2;
            std::cout << "[RepairTip] Tip " << node.first << "->" << t << " multi " << graph[node.first].outgoing_edges[t].at(0).multiplicity << " is merged to edge " << node.first << "->" << max_edge << " with sim " << max_sim << std::endl;
            std::cout << "[RepairTip] Tip " << reverse_complementary_node(t) << "->" << reverse_complementary_node(node.first) << " multi " << graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(t)].at(0).multiplicity << " is merged to edge " << reverse_complementary_node(max_edge) << "->" << reverse_complementary_node(node.first) << " with sim " << max_sim << " multi " << graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).multiplicity << std::endl;
            graph[node.first].outgoing_edges.erase(t);
            graph[t].incoming_edges.erase(node.first);
            graph[reverse_complementary_node(node.first)].incoming_edges.erase(reverse_complementary_node(t));
            graph[reverse_complementary_node(t)].outgoing_edges.erase(reverse_complementary_node(node.first));
        }
    }
    for (auto&& n : nodes_to_remove)
        graph.erase(n);
}

void Graph::merge_tips_into_edges_further_L(unsigned& num_tips, double ratio, std::unordered_map<std::string, std::vector<std::string>> nodes2bc) {
    num_tips = 0;
    // traverse all nodes
    std::set<std::string> nodes_to_remove;
    for (auto&& node : graph) {
        if (nodes_to_remove.find(node.first) != nodes_to_remove.end())
            continue;
        std::vector<std::string> non_tips, tips;
        for (auto&& node_out : node.second.outgoing_edges) {
            if (graph[node_out.first].outgoing_edges.empty() && graph[node_out.first].incoming_edges.size() == 1) {
                if (node_out.second.size() == 1 && nodes2bc.find(graph[node_out.first].sequence) == nodes2bc.end())
                    tips.push_back(node_out.first);
            }
            else {
                if (node_out.second.size() == 1)
                    non_tips.push_back(node_out.first);
            }
        }
        if (tips.empty())
            continue;

        // only consider the special case
        if (tips.size() != 1 || non_tips.size() != 1)
            continue;

        std::string tip = tips.at(0);
        std::string non_tip = non_tips.at(0);
        Edge tip_edge = graph[node.first].outgoing_edges[tip][0];
        Edge non_tip_edge = graph[node.first].outgoing_edges[non_tip][0];

        // if (tip_edge.length >= 1000000)
        //     continue; // skip long tips

        // this case should be resolved by merge_tips_into_edges 
        if (tip_edge.length * 0.8 <= non_tip_edge.length)
            continue;

        // filter out cases that are not suitable for merging
        if (graph[non_tip].outgoing_edges.size() != 1)
            continue;
        std::string out_node_non_tip;
        for (auto&& n_o : graph[non_tip].outgoing_edges)
            out_node_non_tip = n_o.first;
        if (graph[non_tip].outgoing_edges[out_node_non_tip].size() != 1)
            continue;

        Path path;
        add_node_to_path(path, node.first);
        add_node_to_path(path, non_tip);
        add_node_to_path(path, out_node_non_tip);

        std::cout << "[RepairTip] Check tip " << node.first << "->" << tip << " len " << tip_edge.sequence.size() << " and " << node.first << "->" << non_tip << "->" << out_node_non_tip << " len " << path.sequence.size() << std::endl;

        // if the path length is shorter than the tip length, do not merge
        if (tip_edge.length * 0.6 > path.length)
            continue;

        int prefix_size = std::max(1000000, int(non_tip_edge.sequence.size()) + 100000);

        // global alignment of long sequences will be very time consuming
        if (prefix_size > 2000000)
            continue;

        std::string prefix_tip = tip_edge.sequence.substr(0, prefix_size);
        std::string prefix_edge = path.sequence.substr(0, prefix_size);

        double sim = matches_by_edlib(prefix_tip, prefix_edge);
        std::cout << "[RepairTip] Check tip " << node.first << "->" << tip << " multi " << tip_edge.multiplicity << " len " << tip_edge.length << " to " << node.first << "->" << non_tip << "->" << out_node_non_tip << " multi " << path.multiplicity << " len " << path.length << ": sim " << sim << std::endl;
        if (sim < ratio)
            continue;

        nodes_to_remove.insert(tip);
        nodes_to_remove.insert(reverse_complementary_node(tip));
        num_tips += 2;
        std::cout << "[RepairTip] Tip " << node.first << "->" << tip << " multi " << tip_edge.multiplicity << " is merged to path " << node.first << "->" << non_tip << "->" << out_node_non_tip << " with sim " << sim << std::endl;
        std::cout << "[RepairTip] Tip " << reverse_complementary_node(tip) << "->" << reverse_complementary_node(node.first) << " multi " << tip_edge.multiplicity << " is merged to path " << reverse_complementary_node(out_node_non_tip) << "->" << reverse_complementary_node(non_tip) << "->" << reverse_complementary_node(node.first) << " with sim " << sim << std::endl;
        graph[node.first].outgoing_edges.erase(tip);
        graph[tip].incoming_edges.erase(node.first);
        graph[reverse_complementary_node(node.first)].incoming_edges.erase(reverse_complementary_node(tip));
        graph[reverse_complementary_node(tip)].outgoing_edges.erase(reverse_complementary_node(node.first));
    }
    for (auto&& n : nodes_to_remove)
        graph.erase(n);
}

void Graph::remove_deadend_edges_L(unsigned& removed_edges, std::unordered_map<std::string, std::vector<std::string>> nodes2bc) {
    removed_edges = 0;
    std::set<std::string> nodes_to_remove;
    for (auto&& node : graph) {
        std::string deadend_out;
        bool remove = false;
        double sim = 0;

        for (auto&& n_o : node.second.outgoing_edges) {
            deadend_out = n_o.first;

            // No loop deadends
            if (deadend_out == node.first || deadend_out == reverse_complementary_node(node.first))
                continue;

            // allow no bulges between node.first and deadend_out
            if (node.second.outgoing_edges[deadend_out].size() > 1)
                continue;

            Edge& e_deadend = node.second.outgoing_edges[deadend_out][0];

            // node.first --> n_out is a valid deadend edge
            if (graph[deadend_out].outgoing_edges.empty() && nodes2bc.find(node.second.sequence) == nodes2bc.end() && nodes2bc.find(graph[deadend_out].sequence) == nodes2bc.end()) {
                remove = true;
                break;
            }
        }
        if (remove) {
            std::cout << "[RepairTip] Deadend edge " << node.first << " -> " << deadend_out << " is removed (sim=" << sim << ")" << std::endl;
            std::cout << "[RepairTip] Deadend edge " << reverse_complementary_node(deadend_out) << " -> " << reverse_complementary_node(node.first) << " is removed (sim=" << sim << ")" << std::endl;
            graph[node.first].outgoing_edges.erase(deadend_out);
            graph[deadend_out].incoming_edges.erase(node.first);
            graph[reverse_complementary_node(deadend_out)].outgoing_edges.erase(reverse_complementary_node(node.first));
            graph[reverse_complementary_node(node.first)].incoming_edges.erase(reverse_complementary_node(deadend_out));

            removed_edges += 2;

            if (graph[node.first].outgoing_edges.empty() && graph[node.first].incoming_edges.empty()) {
                nodes_to_remove.insert(node.first);
                nodes_to_remove.insert(reverse_complementary_node(node.first));
            }

            if (graph[deadend_out].outgoing_edges.empty() && graph[deadend_out].incoming_edges.empty()) {
                nodes_to_remove.insert(deadend_out);
                nodes_to_remove.insert(reverse_complementary_node(deadend_out));
            }
        }
    }

    for (auto&& n : nodes_to_remove)
        graph.erase(n);

    merge_non_branching_paths(true);

}

void Graph::write_graph_L(const std::string& prefix, int thick, bool contracted, bool colored, std::unordered_set<std::string> nodes_retain, std::unordered_map<std::string, std::vector<std::string>> nodes2bc) {
    std::string graph_dot = prefix + ".dot";
    std::string graph_fasta = prefix + ".fasta";
    std::string graph_path = prefix + ".path";

    std::cout << "[WriteGraph] Write graph " << graph_dot << ", fasta " << graph_fasta << std::endl;
    std::ofstream file_dot(graph_dot);
    std::ofstream file_fasta(graph_fasta);
    std::ofstream file_path(graph_path);

    // extract single edges
    std::unordered_set<std::string> nodes2remove;
    for (auto&& n : graph) {
        if (n.second.outgoing_edges.size() == 1 && n.second.incoming_edges.size() == 1 && n.second.outgoing_edges.find(n.first) != n.second.outgoing_edges.end())
            nodes2remove.insert(n.first);
        if (n.second.outgoing_edges.size() != 1)
            continue;
        if (n.second.incoming_edges.size() != 0)
            continue;

        std::string next_node;
        for (auto&& n_o : n.second.outgoing_edges) {
            next_node = n_o.first;
        }
        if (graph[next_node].outgoing_edges.size() != 0)
            continue;
        if (graph[next_node].incoming_edges.size() != 1)
            continue;
        if (n.second.outgoing_edges[next_node].size() != 1)
            continue;

        if (nodes2bc.find(n.second.sequence) == nodes2bc.end() || nodes2bc.find(graph[next_node].sequence) == nodes2bc.end()) {
            nodes2remove.insert(n.first);
            nodes2remove.insert(next_node);
        }
    }

    for (auto&& n : nodes2remove) {
        graph.erase(n);
    }

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
        else if (node.first.find("+") != std::string::npos) {
            file_dot << "\"" << node.first << "\" [style=filled fillcolor=\"white\" label=\"" << node.first + "_L" + std::to_string(graph[node.first].sequence.size()) << "\"]\n";
        }
        else {
            if (nodes2bc.find(node.second.sequence) != nodes2bc.end()) {
                std::string bc = "";
                for (auto&& b : nodes2bc[node.second.sequence])
                    bc += ("\\nBC: " + b);
                file_dot << node.first << " [style=filled fillcolor=\"white\" label=\"" << node.first + "_L" + std::to_string(graph[node.first].sequence.size()) << bc << "\"]\n";
            }
            else
                file_dot << node.first << " [style=filled fillcolor=\"white\" label=\"" << node.first + "_L" + std::to_string(graph[node.first].sequence.size()) << "\"]\n";
        }
        if (node.second.number_of_contracted_edge > max_contracted) {
            max_contracted_node = get_contracted_name(node.first);
            max_contracted = node.second.number_of_contracted_edge;
        }
    }
    if (max_contracted > 0)
        std::cout << "[WriteGraph] Max contracted node: " << max_contracted_node << " has " << max_contracted << " edges." << std::endl;

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
                            file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << edge.multiplicity << ")\" color=\"black\"]\n";
                        else {
                            file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << edge.multiplicity << ")";
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
                        file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << edge.multiplicity << ")\" color=\"black\"]\n";
                    else {
                        file_dot << "\"" << start_node << "\" -> \"" << end_node << "\" [label=\"" << edge.label << " " << edge.start_base << " " << edge.length << "(" << edge.multiplicity << ")";
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
    std::cout << "[WriteGraph] Total number of nodes: " << this->get_num_nodes() << std::endl;
    std::cout << "[WriteGraph] Total number of edges: " << num_edges << std::endl;
}

void Graph::write_graph_contracted_L(const std::string& prefix, int min_length, bool simplify, std::unordered_map<std::string, std::vector<std::string>> nodes2bc) {
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

            if (nodes2bc.find(graph[node.first].sequence) != nodes2bc.end() || nodes2bc.find(graph[edge.first].sequence) != nodes2bc.end())
                continue;

            // skip tips
            if (node.second.outgoing_edges.empty() || node.second.incoming_edges.empty())
                continue;
            if (graph_vis[edge.first].outgoing_edges.empty() || graph_vis[edge.first].incoming_edges.empty())
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
                    // update_length means update length only, and does not change the graph
                    if (!update_length) {
                        if (graph_vis[node].sequence.empty())
                            graph_vis[node].sequence += ("NNNNNNNNNNNNNNNNNNNN" + e.sequence + "NNNNNNNNNNNNNNNNNNNN");
                        else
                            graph_vis[node].sequence += (e.sequence + "NNNNNNNNNNNNNNNNNNNN");
                    }
                    total_len += e.length;
                    lens.push_back(e.length);
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
                        graph_vis[node].circles = lens.size();
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
        if (nodes2bc.find(graph[node.first].sequence) != nodes2bc.end())
            continue;
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

                if (nodes2bc.find(graph[n_in.first].sequence) != nodes2bc.end())
                    continue;

                // skip tips
                if (graph_vis[n_in.first].outgoing_edges.empty() || graph_vis[n_in.first].incoming_edges.empty())
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

                if (nodes2bc.find(graph[n_out.first].sequence) != nodes2bc.end())
                    continue;

                // skip tips
                if (graph_vis[n_out.first].outgoing_edges.empty() || graph_vis[n_out.first].incoming_edges.empty())
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
            if (nodes2bc.find(graph[node.first].sequence) != nodes2bc.end())
                continue;
            contract_self_loops(node.first, true);
        }
    }

    for (auto&& node : graph_vis) {
        if (nodes2bc.find(graph[node.first].sequence) != nodes2bc.end())
            continue;
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
                    std::cout << "[WriteGraph] Modify contracted node sequence for " << node.first << ", new length " << graph_vis[node.first].sequence.size() << std::endl;
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
        this->write_graph_L(prefix, 1000000, true, true, std::unordered_set<std::string>(), nodes2bc);
    }
    else {
        auto graph_cp = graph;
        this->graph = graph_vis;
        this->write_graph_L(prefix, 1000000, true, true, std::unordered_set<std::string>(), nodes2bc);
        this->graph = graph_cp;
    }
    return;
}

void Graph::write_prefix_siffux_linear_edges(std::string output, std::string jumbodbg, int threads, int k_mer, unsigned& num_glued) {
    num_glued = 0;
    execute_command("mkdir -p " + output);
    std::string output_prefix_suffix = output + "/graph_linear_edges.fa";
    std::unordered_set<std::string> traversed_nodes;
    std::ofstream outfile(output_prefix_suffix);
    if (!outfile.is_open()) {
        throw std::runtime_error("Failed to open file: " + output_prefix_suffix);
    }

    std::cout << "[Connect] Write prefixes (and suffixes) of linear edges (and tips) to fasta" << std::endl;
    int bc_id = 1;
    std::unordered_map<std::string, std::vector<std::string>> kmer2bc;
    std::unordered_map<std::string, std::vector<std::string>> kmer2nodes;

    for (auto&& node : graph) {
        if (traversed_nodes.find(node.first) != traversed_nodes.end())
            continue;
        if (node.second.outgoing_edges.size() == 1) {
            std::string n_out;
            for (auto&& n : node.second.outgoing_edges)
                n_out = n.first;

            if (graph[n_out].incoming_edges.size() != 1 || !graph[n_out].outgoing_edges.empty())
                continue;

            if (node.second.outgoing_edges[n_out].size() != 1)
                continue;

            traversed_nodes.insert(node.first);
            traversed_nodes.insert(reverse_complementary_node(n_out));

            std::string edge_label_forward = "Edge: " + node.second.outgoing_edges[n_out].at(0).label + " Length: " + format_with_commas(node.second.outgoing_edges[n_out].at(0).sequence.size());
            if (!node.second.outgoing_edges[n_out].at(0).ref_ids.empty()) {
                for (auto&& r : node.second.outgoing_edges[n_out].at(0).ref_ids)
                    edge_label_forward += ("\\n" + r);
            }

            std::string edge_label_reverse = "Edge: " + graph[reverse_complementary_node(n_out)].outgoing_edges[reverse_complementary_node(node.first)].at(0).label + " Length: " + format_with_commas(node.second.outgoing_edges[n_out].at(0).sequence.size());
            if (!graph[reverse_complementary_node(n_out)].outgoing_edges[reverse_complementary_node(node.first)].at(0).ref_ids.empty()) {
                for (auto&& r : graph[reverse_complementary_node(n_out)].outgoing_edges[reverse_complementary_node(node.first)].at(0).ref_ids)
                    edge_label_reverse += ("\\n" + r);
            }

            // find a valid linear edge
            if (node.second.incoming_edges.empty()) {
                std::string seq1 = replaceNsWithRandomBases(node.second.sequence.substr(0, 5000));
                outfile << ">" << node.first << "\n";
                outfile << seq1 << "\n";

                std::string seq2 = replaceNsWithRandomBases(graph[n_out].sequence.substr(graph[n_out].sequence.size() - 5000));
                outfile << ">" << n_out << "\n";
                outfile << seq2 << "\n";

                std::string suffix_start = seq1.substr(seq1.size() - k_mer);
                std::string prefix_end = seq2.substr(0, k_mer);
                kmer2bc[suffix_start].push_back(std::to_string(bc_id) + " " + edge_label_forward);
                kmer2bc[prefix_end].push_back(std::to_string(-bc_id) + " " + edge_label_forward);
                kmer2nodes[suffix_start].push_back(node.first);
                kmer2nodes[prefix_end].push_back(n_out);

                bc_id++;
                kmer2bc[reverse_complementary(suffix_start)].push_back(std::to_string(-bc_id) + " " + edge_label_reverse);
                kmer2bc[reverse_complementary(prefix_end)].push_back(std::to_string(bc_id) + " " + edge_label_reverse);
                kmer2nodes[reverse_complementary(suffix_start)].push_back(reverse_complementary_node(node.first));
                kmer2nodes[reverse_complementary(prefix_end)].push_back(reverse_complementary_node(n_out));

                bc_id++;
            }
            else {
                std::string seq2 = replaceNsWithRandomBases(graph[n_out].sequence.substr(graph[n_out].sequence.size() - 5000));
                outfile << ">" << n_out << "\n";
                outfile << seq2 << "\n";

                std::string prefix_end = seq2.substr(0, k_mer);
                kmer2bc[prefix_end].push_back(std::to_string(-bc_id) + " " + edge_label_forward);
                kmer2nodes[prefix_end].push_back(n_out);

                bc_id++;
                kmer2bc[reverse_complementary(prefix_end)].push_back(std::to_string(bc_id) + " " + edge_label_reverse);
                kmer2nodes[reverse_complementary(prefix_end)].push_back(reverse_complementary_node(n_out));

                bc_id++;
            }
        }
    }
    outfile.close();

    // for (auto&& bc : bc2edge) {
    //     std::cout << bc.first << ": " << bc.second << std::endl;
    // }

    if (!fs::is_directory(output + "/graph_linear_edges.dbg"))
        execute_command(jumbodbg + " --reads " + output + "/graph_linear_edges.fa" + " -t " + std::to_string(threads) + " --coverage -k " + std::to_string(k_mer) + " -o " + output + "/graph_linear_edges.dbg");

    unsigned removed_paths = 1;
    unsigned removed_whirls = 1;
    unsigned removed_bulges = 1;
    int decoupled = 1;
    unsigned cnt_rounds = 0;
    int total_removed = 0;
    unsigned total_whirls = 0;
    unsigned removed_edges = 1;
    unsigned removed_tips = 1;

    Graph graph_linear_edges;
    graph_linear_edges.restart_from_dot(output + "/graph_linear_edges.dbg/graph.dot", output + "/graph_linear_edges.dbg/graph.fasta", k_mer);
    // graph_linear_edges.write_graph(output + "/graph_linear_edges.ori");
    // graph_linear_edges.get_annotation(output + "/graph_linear_edges.ori");
    // graph_linear_edges.write_graph(output + "/graph_linear_edges.ori.color", 1000000, false, true, std::unordered_set<std::string>(), kmer2bc);

    removed_bulges = 1;
    total_removed = 0;
    while (removed_bulges) {
        graph_linear_edges.multi_bulge_removal(removed_bulges);
        total_removed += removed_bulges;
    }

    removed_whirls = 1;
    total_removed = 0;
    while (removed_whirls) {
        graph_linear_edges.general_whirl_removal(removed_whirls);
        total_removed += removed_whirls;
    }

    removed_paths = 1;
    total_removed = 0;
    while (removed_paths) {
        graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 3, 0.8, true, 3);
        total_removed += removed_paths;
    }
    removed_paths = 1;
    while (removed_paths) {
        graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 4, 0.8, true, 2);
        total_removed += removed_paths;
    }
    removed_paths = 1;
    while (removed_paths) {
        graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.8, true, 2);
        total_removed += removed_paths;
    }

    removed_paths = 1;
    while (removed_paths) {
        graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0.6, true, 2);
        total_removed += removed_paths;
    }

    // ensure all below outputting graph_linear_edges have no simple bulges, or the program will fail
    decoupled = 1;
    total_removed = 0;
    while (decoupled) {
        graph_linear_edges.resolve_edges_in_reverse_complement(decoupled, true);
        total_removed += decoupled;
    }

    total_removed = 0;
    removed_tips = 1;
    while (removed_tips) {
        graph_linear_edges.merge_tips_into_edges_L(removed_tips, 0.8, false, false, kmer2bc);
        total_removed += removed_tips;
        if (removed_tips > 0)
            std::cout << "Merged " << removed_tips << " tips to edges" << std::endl;
    }

    while (true) {
        bool flag = true;
        removed_paths = 1;
        while (removed_paths) {
            graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0.6, true, 2);
            if (removed_paths)
                flag = false;
            if (removed_paths > 0)
                std::cout << "[Connect] Detoured " << removed_paths << " paths" << std::endl;
        }

        decoupled = 1;
        while (decoupled) {
            graph_linear_edges.resolve_edges_in_reverse_complement(decoupled, true);
            if (decoupled)
                flag = false;
            if (decoupled > 0)
                std::cout << "[Connect] Decoupled " << decoupled << " strands" << std::endl;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph_linear_edges.merge_tips_into_edges_L(removed_tips, 0.8, false, false, kmer2bc);
            if (removed_tips)
                flag = false;
            if (removed_tips > 0)
                std::cout << "[Connect] Merged " << removed_tips << " tips to edges" << std::endl;
        }

        removed_whirls = 1;
        while (removed_whirls) {
            graph_linear_edges.general_whirl_removal(removed_whirls);
            graph_linear_edges.merge_non_branching_paths(true);
            if (removed_whirls)
                flag = false;
            if (removed_whirls > 0)
                std::cout << "[Connect] Removed " << removed_whirls << " whirls" << std::endl;
        }

        removed_bulges = 1;
        while (removed_bulges) {
            graph_linear_edges.multi_bulge_removal(removed_bulges);
            graph_linear_edges.merge_non_branching_paths(true);
            if (removed_bulges)
                flag = false;
            if (removed_bulges > 0)
                std::cout << "[Connect] Removed " << removed_bulges << " bulges" << std::endl;
        }

        if (flag)
            break;
    }

    while (true)
    {
        bool flag = true;
        removed_paths = 1;
        while (removed_paths) {
            graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0.6, true, 2, true);
            if (removed_paths)
                flag = false;
            if (removed_paths > 0)
                std::cout << "[Connect] Detoured " << removed_paths << " paths" << std::endl;
        }

        decoupled = 1;
        while (decoupled) {
            graph_linear_edges.resolve_edges_in_reverse_complement(decoupled, true);
            if (decoupled)
                flag = false;
            if (decoupled > 0)
                std::cout << "[Connect] Decoupled " << decoupled << " strands" << std::endl;
        }

        if (flag)
            break;
    }

    while (true) {
        bool flag = true;
        removed_paths = 1;
        while (removed_paths) {
            graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0, true, 2, true);
            if (removed_paths)
                flag = false;
        }
        decoupled = 1;
        while (decoupled) {
            graph_linear_edges.resolve_edges_in_reverse_complement(decoupled, false);
            if (decoupled)
                flag = false;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph_linear_edges.merge_tips_into_edges_L(removed_tips, 0.8, false, false, kmer2bc);
            if (removed_tips)
                flag = false;
        }

        removed_whirls = 1;
        while (removed_whirls) {
            graph_linear_edges.general_whirl_removal(removed_whirls);
            graph_linear_edges.merge_non_branching_paths(true);
            if (removed_whirls)
                flag = false;
        }

        removed_bulges = 1;
        while (removed_bulges) {
            graph_linear_edges.multi_bulge_removal(removed_bulges);
            graph_linear_edges.merge_non_branching_paths(true);
            if (removed_bulges)
                flag = false;
        }

        if (flag)
            break;
    }

    while (true) {
        bool flag = true;

        removed_paths = 1;
        while (removed_paths) {
            graph_linear_edges.resolving_bulge_with_two_multi_edge_paths(removed_paths, 8, 0, true, 2, true);
            if (removed_paths)
                flag = false;
        }

        decoupled = 1;
        while (decoupled) {
            graph_linear_edges.resolve_edges_in_reverse_complement(decoupled);
            if (decoupled)
                flag = false;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph_linear_edges.merge_tips_into_edges_L(removed_tips, 0.8, false, true, kmer2bc);
            if (removed_tips)
                flag = false;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph_linear_edges.merge_tips_L(removed_tips, kmer2bc);
            if (removed_tips)
                flag = false;
        }

        removed_tips = 1;
        while (removed_tips) {
            graph_linear_edges.merge_tips_into_edges_further_L(removed_tips, 0.8, kmer2bc);
            if (removed_tips)
                flag = false;
        }

        removed_whirls = 1;
        while (removed_whirls) {
            graph_linear_edges.general_whirl_removal(removed_whirls, false, true);
            graph_linear_edges.merge_non_branching_paths(true);
            if (removed_whirls)
                flag = false;
        }

        removed_bulges = 1;
        while (removed_bulges) {
            graph_linear_edges.multi_bulge_removal(removed_bulges);
            graph_linear_edges.merge_non_branching_paths(true);
            if (removed_bulges)
                flag = false;
        }

        removed_edges = 1;
        while (removed_edges) {
            graph_linear_edges.remove_deadend_edges_L(removed_edges, kmer2bc);
            if (removed_edges)
                flag = false;
        }

        if (flag)
            break;
    }

    graph_linear_edges.write_graph_L(output + "/graph_linear_edges.final", 1000000, false, true, std::unordered_set<std::string>(), kmer2bc);
    graph_linear_edges.write_graph_contracted_L(output + "/graph_linear_edges.final.contracted.600", 600, false, kmer2bc);

    for (auto&& n : graph_linear_edges.graph) {
        if (n.second.incoming_edges.empty() && n.second.outgoing_edges.size() == 1) {
            std::string n_out;
            for (auto&& n_o : n.second.outgoing_edges)
                n_out = n_o.first;

            if (n.second.outgoing_edges[n_out].size() != 1)
                continue;

            if (kmer2nodes.find(n.second.sequence) == kmer2nodes.end())
                continue;

            if (kmer2nodes.find(graph_linear_edges.graph[n_out].sequence) == kmer2nodes.end())
                continue;

            int max_i = -1;
            int max_len_node1 = 0;
            for (int i = 0; i < kmer2nodes.at(n.second.sequence).size(); ++i) {
                std::string node1 = kmer2nodes.at(n.second.sequence).at(i);
                std::string node1_in;
                if (!graph[node1].outgoing_edges.empty())
                    continue;
                for (auto&& n1_i : graph[node1].incoming_edges)
                    node1_in = n1_i.first;

                int len = graph[node1_in].outgoing_edges[node1].at(0).sequence.size();
                if (len > max_len_node1) {
                    max_i = i;
                    max_len_node1 = len;
                }

            }

            int max_j = -1;
            int max_len_node2 = 0;
            for (int j = 0; j < kmer2nodes.at(graph_linear_edges.graph[n_out].sequence).size(); ++j) {
                std::string node2 = kmer2nodes.at(graph_linear_edges.graph[n_out].sequence).at(j);
                std::string node2_out;
                if (!graph[node2].incoming_edges.empty())
                    continue;
                for (auto&& n2_o : graph[node2].outgoing_edges)
                    node2_out = n2_o.first;

                int len = graph[node2].outgoing_edges[node2_out].at(0).sequence.size();
                if (len > max_len_node2) {
                    max_j = j;
                    max_len_node2 = len;
                }

            }

            if (max_i == -1 || max_j == -1)
                continue;

            std::string node1 = kmer2nodes.at(n.second.sequence).at(max_i);
            std::string node2 = kmer2nodes.at(graph_linear_edges.graph[n_out].sequence).at(max_j);

            if (node1 == reverse_complementary_node(node2))
                continue;

            // node1 and node 2 can be glued
            if (graph[node1].outgoing_edges.empty() && graph[node2].incoming_edges.empty() && graph[node2].outgoing_edges.find(node1) == graph[node2].outgoing_edges.end()) {
                std::string node1_in, node2_out;
                assert(graph[node1].incoming_edges.size() == 1);
                for (auto&& n1_i : graph[node1].incoming_edges)
                    node1_in = n1_i.first;
                assert(graph[node2].outgoing_edges.size() == 1);
                for (auto n2_o : graph[node2].outgoing_edges)
                    node2_out = n2_o.first;

                if (node1 == reverse_complementary_node(node2_out) || node2 == reverse_complementary_node(node1_in))
                    continue;

                std::cout << "[Connect] Glue node " << node1 << " and " << node2 << " based on edge " << n.first << " -> " << n_out << std::endl;

                // delete last 5000 bp for node1 and first 5000 bp for node2
                graph[node1].sequence = graph[node1].sequence.substr(0, graph[node1].sequence.size() - 5000);
                assert(graph[node1_in].outgoing_edges[node1].size() == 1);
                graph[node1_in].outgoing_edges[node1].at(0).sequence = graph[node1_in].outgoing_edges[node1].at(0).sequence.substr(0, graph[node1_in].outgoing_edges[node1].at(0).sequence.size() - 5000);
                graph[node1_in].outgoing_edges[node1].at(0).length = graph[node1_in].outgoing_edges[node1].at(0).sequence.size();
                graph[node1].incoming_edges[node1_in].at(0).sequence = graph[node1_in].outgoing_edges[node1].at(0).sequence.substr(0, graph[node1_in].outgoing_edges[node1].at(0).sequence.size() - 5000);
                graph[node1].incoming_edges[node1_in].at(0).length = graph[node1].incoming_edges[node1_in].at(0).sequence.size();

                graph[node2].sequence = graph[node2].sequence.substr(5000);
                assert(graph[node2].outgoing_edges[node2_out].size() == 1);
                graph[node2].outgoing_edges[node2_out].at(0).sequence = graph[node2].outgoing_edges[node2_out].at(0).sequence.substr(5000);
                graph[node2].outgoing_edges[node2_out].at(0).length = graph[node2].outgoing_edges[node2_out].at(0).sequence.size();
                graph[node2_out].incoming_edges[node2].at(0).sequence = graph[node2].outgoing_edges[node2_out].at(0).sequence.substr(5000);
                graph[node2_out].incoming_edges[node2].at(0).length = graph[node2_out].incoming_edges[node2].at(0).sequence.size();

                // determine the edge sequence
                std::string seq = n.second.outgoing_edges[n_out].at(0).sequence;
                std::string new_edge_seq = graph[node1].sequence + seq + graph[node2].sequence;

                Edge edge(new_edge_seq.at(graph[node1].sequence.size()), new_edge_seq.size(), new_edge_seq, mean_cov);
                edge.path_edges_in_original_graph.push_back(node1 + "_" + node2);
                edge.path_nodes_in_original_graph.push_back(node1);
                edge.path_nodes_in_original_graph.push_back(node2);

                graph[node1].outgoing_edges[node2].push_back(edge);
                graph[node2].incoming_edges[node1].push_back(edge);
                num_glued += 1;
            }
        }
    }

    merge_non_branching_paths(true);
}

void Graph::remove_contained_contigs_minimap(std::string output, int threads, unsigned& num_contained) {
    num_contained = 0;
    execute_command("mkdir -p " + output);
    std::string output_prefix_suffix = output + "/linear_prefix_and_suffix.fasta";
    std::ofstream outfile_extracted(output_prefix_suffix);
    if (!outfile_extracted.is_open()) {
        throw std::runtime_error("Failed to open file: " + output_prefix_suffix);
    }

    std::string output_all = output + "/linear_all.fasta";
    std::ofstream outfile_all(output_all);
    if (!outfile_all.is_open()) {
        throw std::runtime_error("Failed to open file: " + output_all);
    }

    std::unordered_map<std::string, std::string> print2node;

    std::unordered_set<std::string> traversed_nodes;
    std::cout << "[Deduplicate] Write prefixes and suffixes of linear edges to fasta" << std::endl;
    int extract_length = 1000000;
    for (auto&& node : graph) {
        print2node[get_contracted_name(node.first)] = node.first;
        if (traversed_nodes.find(node.first) != traversed_nodes.end())
            continue;
        if (node.second.outgoing_edges.size() == 1 && node.second.incoming_edges.size() == 1 && node.second.outgoing_edges.find(node.first) != node.second.outgoing_edges.end()) {
            if (node.second.outgoing_edges[node.first].size() > 1)
                continue;

            traversed_nodes.insert(node.first);
            traversed_nodes.insert(reverse_complementary_node(node.first));

            std::string seq = node.second.outgoing_edges[node.first].at(0).sequence;
            outfile_extracted << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(node.first) << "_0" << "\n";
            outfile_extracted << seq << "\n";

            outfile_all << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(node.first) << "\n";
            outfile_all << seq << "\n";
        }
        if (node.second.outgoing_edges.size() == 1) {
            std::string n_out;
            for (auto&& n : node.second.outgoing_edges)
                n_out = n.first;

            if (graph[n_out].incoming_edges.size() != 1 || !graph[n_out].outgoing_edges.empty())
                continue;

            if (node.second.outgoing_edges[n_out].size() != 1)
                continue;

            traversed_nodes.insert(node.first);
            traversed_nodes.insert(reverse_complementary_node(n_out));

            // find a valid linear edge
            if (node.second.incoming_edges.empty()) {
                std::string seq = node.second.outgoing_edges[n_out].at(0).sequence;
                if (seq.size() <= extract_length) {
                    std::string seq1 = seq.substr(0, extract_length);
                    outfile_extracted << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(n_out) << "_0" << "\n";
                    outfile_extracted << seq1 << "\n";
                }
                else {
                    std::string seq1 = seq.substr(0, extract_length);
                    outfile_extracted << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(n_out) << "_1" << "\n";
                    outfile_extracted << seq1 << "\n";

                    std::string seq2 = seq.substr(seq.size() - extract_length);
                    outfile_extracted << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(n_out) << "_2" << "\n";
                    outfile_extracted << seq2 << "\n";
                }

                outfile_all << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(n_out) << "\n";
                outfile_all << seq << "\n";
            }
        }
    }
    outfile_extracted.close();
    outfile_all.close();

    std::string out_bam_prefix = output + "/align.prefix_suffix.all";;

    if (!fs::is_regular_file(out_bam_prefix + ".bam")) {
        execute_command(("minimap2 -ax asm20 --eqx -Y -p 0.1 " + output_all + " " + output_prefix_suffix + " -t " + std::to_string(threads) + " 2>/dev/null | grep -v '^@' > " + out_bam_prefix + ".sam").c_str());
        if (!std::filesystem::exists(output_all + ".fai")) {
            execute_command(("samtools faidx " + output_all).c_str());
        }
        execute_command(("cut -f1,2 " + output_all + ".fai | awk " + R"('{print "@SQ\tSN:"$1"\tLN:"$2}')" + " > " + out_bam_prefix + ".header.sam").c_str());
        execute_command(("cat " + out_bam_prefix + ".header.sam " + out_bam_prefix + ".sam | samtools sort -@ " + std::to_string(threads) + " -o " + out_bam_prefix + ".bam").c_str());
    }

    std::string exeDir = getExecutablePath();
    execute_command((exeDir + "/../src/scripts/remove_cognate.py -o " + out_bam_prefix + ".results " + out_bam_prefix + ".bam").c_str());
    std::ifstream infile(out_bam_prefix + ".results");

    // Check if file opened successfully
    if (!infile.is_open()) {
        throw std::runtime_error("Error opening file: " + out_bam_prefix + ".results");
    }

    std::string line;
    while (std::getline(infile, line)) {
        size_t underscore_pos = line.find('_');
        if (underscore_pos != std::string::npos) {
            std::string node1 = print2node.at(line.substr(0, underscore_pos));
            std::string node2 = print2node.at(line.substr(underscore_pos + 1));
            std::cout << "[Deduplicate] Remove contained edge " << node1 << " -> " << node2 << std::endl;
            std::cout << "[Deduplicate] Remove contained edge " << reverse_complementary_node(node2) << " -> " << reverse_complementary_node(node1) << std::endl;

            graph.erase(node1);
            graph.erase(node2);
            graph.erase(reverse_complementary_node(node2));
            graph.erase(reverse_complementary_node(node1));
            num_contained += 2;
        }
    }

    infile.close();
}

void Graph::connect_linear_and_tips_using_spanning_reads(std::string output, int threads, std::string reads, double identity) {
    execute_command("mkdir -p " + output);

    std::string output_all = output + "/linear_all.fasta";
    std::ofstream outfile_all(output_all);
    if (!outfile_all.is_open()) {
        throw std::runtime_error("Failed to open file: " + output_all);
    }

    std::unordered_map<std::string, std::string> print2node;
    std::unordered_set<std::string> traversed_nodes;
    for (auto&& node : graph) {
        print2node[get_contracted_name(node.first)] = node.first;
        if (traversed_nodes.find(node.first) != traversed_nodes.end())
            continue;
        if (node.second.outgoing_edges.size() == 1 && node.second.incoming_edges.size() == 1 && node.second.outgoing_edges.find(node.first) != node.second.outgoing_edges.end()) {
            if (node.second.outgoing_edges[node.first].size() > 1)
                continue;

            traversed_nodes.insert(node.first);
            traversed_nodes.insert(reverse_complementary_node(node.first));

            std::string seq = node.second.outgoing_edges[node.first].at(0).sequence;
            outfile_all << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(node.first) << "\n";
            outfile_all << seq << "\n";
        }
        if (node.second.outgoing_edges.size() == 1) {
            std::string n_out;
            for (auto&& n : node.second.outgoing_edges)
                n_out = n.first;

            if (graph[n_out].incoming_edges.size() != 1 || !graph[n_out].outgoing_edges.empty())
                continue;

            if (node.second.outgoing_edges[n_out].size() != 1)
                continue;

            traversed_nodes.insert(node.first);
            traversed_nodes.insert(reverse_complementary_node(n_out));

            // find a valid linear edge
            if (node.second.incoming_edges.empty()) {
                std::string seq = node.second.outgoing_edges[n_out].at(0).sequence;
                outfile_all << ">" << get_contracted_name(node.first) << "_" << get_contracted_name(n_out) << "\n";
                outfile_all << seq << "\n";
            }
        }
    }
    outfile_all.close();

    std::string out_bam = output + "/linear_all.flanks.bam";
    std::string out_result = output + "/linear_all.flanks.results";

    std::string exeDir = getExecutablePath();
    std::string compress = exeDir + "/../lib/LJA/bin/compress";
    std::string spanning_reads_script = exeDir + "/../src/scripts/spanning_reads.py";

    if (!fs::is_regular_file(out_result)) {
        if (execute_command((spanning_reads_script + " " + output_all + " " + reads + " " + out_bam + " " + out_result + " -c " + compress + " -t " + std::to_string(threads) + " -i " + std::to_string(identity)).c_str()) != 0)
            throw std::runtime_error("Failed to execute: " + spanning_reads_script + " " + output_all + " " + reads + " " + out_bam + " " + out_result + " -c " + compress + " -t " + std::to_string(threads));
    }

    std::ifstream infile(out_result);

    // Check if file opened successfully
    if (!infile.is_open()) {
        throw std::runtime_error("Error opening file: " + out_result);
    }

    std::string line;
    struct Connection {
        std::string node1;
        std::string node2;
        std::string node1_in;
        std::string node2_out;
        std::string connecting_seq;
        int len_both_edges;
    };
    std::vector<Connection> connections;
    while (std::getline(infile, line)) {
        size_t underscore_pos = line.find('_');
        size_t space_pos = line.find('\t');
        if (underscore_pos != std::string::npos) {
            std::string node1 = print2node.at(line.substr(0, underscore_pos));
            std::string node2 = print2node.at(line.substr(underscore_pos + 1, space_pos - underscore_pos - 1));
            std::string connecting_seq = line.substr(space_pos + 1);

            if (node1 == reverse_complementary_node(node2))
                continue;

            if (graph.find(node1) == graph.end() || graph[node1].incoming_edges.size() != 1 || graph[node1].outgoing_edges.size() != 0)
                continue;
            if (graph.find(node2) == graph.end() || graph[node2].incoming_edges.size() != 0 || graph[node2].outgoing_edges.size() != 1)
                continue;

            std::string node1_in, node2_out;
            for (auto&& n1_i : graph[node1].incoming_edges)
                node1_in = n1_i.first;
            for (auto n2_o : graph[node2].outgoing_edges)
                node2_out = n2_o.first;

            if (node1 == reverse_complementary_node(node2_out) || node2 == reverse_complementary_node(node1_in))
                continue;

            assert(graph[node1_in].outgoing_edges[node1].size() == 1);
            assert(graph[node2].outgoing_edges[node2_out].size() == 1);

            connections.push_back({ node1, node2, node1_in, node2_out, connecting_seq,
                                   (int)graph[node1_in].outgoing_edges[node1].at(0).sequence.size() + (int)graph[node2].outgoing_edges[node2_out].at(0).sequence.size() });
        }
    }
    infile.close();

    // sort connections by the length of both edges descending
    std::sort(connections.begin(), connections.end(), [](const Connection& a, const Connection& b) {
        return a.len_both_edges > b.len_both_edges;
        });

    for (auto&& conn : connections) {
        if (graph.find(conn.node1) == graph.end() || graph.find(conn.node2) == graph.end())
            continue;

        std::string node1 = conn.node1;
        std::string node2 = conn.node2;

        if (graph.find(node1) == graph.end() || graph[node1].incoming_edges.size() != 1 || graph[node1].outgoing_edges.size() != 0)
            continue;
        if (graph.find(node2) == graph.end() || graph[node2].incoming_edges.size() != 0 || graph[node2].outgoing_edges.size() != 1)
            continue;

        std::string node1_in, node2_out;
        for (auto&& n1_i : graph[node1].incoming_edges)
            node1_in = n1_i.first;
        for (auto n2_o : graph[node2].outgoing_edges)
            node2_out = n2_o.first;

        if (node1 == reverse_complementary_node(node2_out) || node2 == reverse_complementary_node(node1_in))
            continue;

        if (node1 == reverse_complementary_node(node1_in) || node2 == reverse_complementary_node(node2_out))
            continue;

        assert(graph[node1_in].outgoing_edges[node1].size() == 1);
        assert(graph[node2].outgoing_edges[node2_out].size() == 1);

        std::string connecting_seq = conn.connecting_seq;

        Edge edge_in = graph[node1_in].outgoing_edges[node1].at(0);
        Edge edge_out = graph[node2].outgoing_edges[node2_out].at(0);

        std::string edge_in_unique = edge_in.sequence.size() >= 20000 ? edge_in.sequence.substr(0, edge_in.sequence.size() - 20000) : "";
        std::string edge_out_unique = edge_out.sequence.size() >= 20000 ? edge_out.sequence.substr(20000) : "";

        // for forward strand
        std::cout << "[Connect] Connect " << node1 << " and " << node2 << " using spanning reads" << std::endl;
        std::string new_edge_seq = edge_in_unique + connecting_seq + edge_out_unique;
        Edge new_edge(new_edge_seq.at(graph[node1_in].sequence.size()), new_edge_seq.size(), new_edge_seq, std::max(edge_in.multiplicity, edge_out.multiplicity));
        new_edge.path_edges_in_original_graph.push_back(node1 + "_" + node2);
        new_edge.path_nodes_in_original_graph.push_back(node1);
        new_edge.path_nodes_in_original_graph.push_back(node2);
        graph[node1_in].outgoing_edges.erase(node1);
        graph[node1].incoming_edges.erase(node1_in);
        graph[node2].outgoing_edges.erase(node2_out);
        graph[node2_out].incoming_edges.erase(node2);
        graph[node1_in].outgoing_edges[node2_out].push_back(new_edge);
        graph[node2_out].incoming_edges[node1_in].push_back(new_edge);
        graph.erase(node1);
        graph.erase(node2);

        // for reverse strand
        std::cout << "[Connect] Connect " << reverse_complementary_node(node2) << " and " << reverse_complementary_node(node1) << " using spanning reads" << std::endl;
        std::string new_edge_seq_rc = reverse_complementary(new_edge_seq);
        Edge new_edge_rc(new_edge_seq_rc.at(graph[node2_out].sequence.size()), new_edge_seq_rc.size(), new_edge_seq_rc, std::max(edge_in.multiplicity, edge_out.multiplicity));
        new_edge_rc.path_edges_in_original_graph.push_back(reverse_complementary_node(node2) + "_" + reverse_complementary_node(node1));
        new_edge_rc.path_nodes_in_original_graph.push_back(reverse_complementary_node(node2));
        new_edge_rc.path_nodes_in_original_graph.push_back(reverse_complementary_node(node1));
        graph[reverse_complementary_node(node2_out)].outgoing_edges.erase(reverse_complementary_node(node2));
        graph[reverse_complementary_node(node2)].incoming_edges.erase(reverse_complementary_node(node2_out));
        graph[reverse_complementary_node(node1)].outgoing_edges.erase(reverse_complementary_node(node1_in));
        graph[reverse_complementary_node(node1_in)].incoming_edges.erase(reverse_complementary_node(node1));
        graph[reverse_complementary_node(node2_out)].outgoing_edges[reverse_complementary_node(node1_in)].push_back(new_edge_rc);
        graph[reverse_complementary_node(node1_in)].incoming_edges[reverse_complementary_node(node2_out)].push_back(new_edge_rc);
        graph.erase(reverse_complementary_node(node1));
        graph.erase(reverse_complementary_node(node2));
    }
}