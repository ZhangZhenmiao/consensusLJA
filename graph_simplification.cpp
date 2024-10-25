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

bool Graph::check_non_branching(std::string node, bool merge_self_loop) {
    // in any case, 1-in-1-out are non-branching
    if (this->graph[node].incoming_edges.size() == 1 && this->graph[node].outgoing_edges.size() == 1) {
        bool flag = true;
        for (auto&& n : this->graph[node].incoming_edges) {
            if (n.second.size() > 1)
                flag = false;
        }
        for (auto&& n : this->graph[node].outgoing_edges) {
            if (n.second.size() > 1)
                flag = false;
        }
        return flag;
    }
    // 1-in-1-out plus a self-loop
    else if (merge_self_loop && this->graph[node].incoming_edges.find(node) != this->graph[node].incoming_edges.end() && this->graph[node].incoming_edges.size() <= 2 && this->graph[node].outgoing_edges.size() <= 2) {
        bool flag = true;
        for (auto&& n : this->graph[node].incoming_edges) {
            if (n.second.size() > 1)
                flag = false;
        }
        for (auto&& n : this->graph[node].outgoing_edges) {
            if (n.second.size() > 1)
                flag = false;
        }
        return flag;
    }
    else
        return false;
}

bool Graph::add_node_to_path(Path& path, std::string node, int bulge_leg, bool reduce_reverse) {
    if (path.nodes.empty()) {
        // the node should exist
        if (graph.find(node) == graph.end())
            return false;
        path.nodes.push_back(node);
        path.path_nodes_in_original_graph.push_back(node);
        return true;
    }
    std::string prev_node = path.nodes.at(path.nodes.size() - 1);

    // the edge should exist
    if (this->graph[prev_node].outgoing_edges.find(node) == this->graph[prev_node].outgoing_edges.end())
        return false;

    path.nodes.push_back(node);
    path.bulge_legs.push_back(bulge_leg);
    path.path_nodes_in_original_graph.pop_back();
    merge_vecs(path.path_nodes_in_original_graph, graph[prev_node].outgoing_edges[node].at(bulge_leg).path_nodes_in_original_graph);
    merge_vecs(path.path_edges_in_original_graph, graph[prev_node].outgoing_edges[node].at(bulge_leg).path_edges_in_original_graph);
    if (path.path_nodes_in_original_graph.size() != path.path_edges_in_original_graph.size() + 1)
        write_graph("debug_add_node" + node);
    assert(path.path_nodes_in_original_graph.size() == path.path_edges_in_original_graph.size() + 1);
    Edge& edge = this->graph[prev_node].outgoing_edges[node].at(bulge_leg);
    if (!reduce_reverse || prev_node != reverse_complementary_node(node)) {
        path.update_min_multi(edge);
        path.multiplicity = (path.multiplicity * path.length + edge.multiplicity * edge.length) / (path.length + edge.length);
    }
    else {
        if (path.min_multi == 0)
            path.min_multi = edge.multiplicity / 2;
        else
            path.min_multi = std::min(path.min_multi, edge.multiplicity / 2);
        path.multiplicity = (path.multiplicity * path.length + edge.multiplicity * edge.length / 2) / (path.length + edge.length);
    }

    path.safe_to_extract = true;
    for (int i = 1; i < path.nodes.size() - 1; ++i) {
        int bulge_in = path.bulge_legs[i - 1];
        int bulge_out = path.bulge_legs[i];
        bool flag_in = false, flag_out = false;
        if (graph[path.nodes[i]].incoming_edges[path.nodes[i - 1]].at(bulge_in).multiplicity - MIN_MULTI < path.min_multi && graph[path.nodes[i]].incoming_edges[path.nodes[i - 1]].size() <= 1 && graph[path.nodes[i]].incoming_edges.size() <= 1) {
            flag_in = true;
        }
        if (graph[path.nodes[i]].outgoing_edges[path.nodes[i + 1]].at(bulge_out).multiplicity - MIN_MULTI < path.min_multi && graph[path.nodes[i]].outgoing_edges[path.nodes[i + 1]].size() <= 1 && graph[path.nodes[i]].outgoing_edges.size() <= 1) {
            flag_out = true;
        }
        if (flag_in || flag_out)
            path.safe_to_extract = false;

        if (graph[path.nodes[i]].outgoing_edges.find(path.nodes[i]) != graph[path.nodes[i]].outgoing_edges.end())
            path.safe_to_extract = false;
    }

    // update path sequence
    path.length += edge.length;
    if (path.sequence.empty())
        path.sequence = edge.sequence;
    else
        path.sequence += edge.sequence.substr(this->k);

    if (path.nodes.size() >= 3 && this->graph[prev_node].outgoing_edges.size() > 1)
        path.is_unambiguous = false;

    if (prev_node == reverse_complementary_node(node)) {
        path.index2palindromic_seq[path.nodes.size() - 2] = this->graph[prev_node].outgoing_edges[node].at(bulge_leg).sequence;
    }

    // if (!path.index2palindromic_seq.empty())
    //     path.safe_to_extract = false;
    return true;
}

std::string Graph::merge_edges(std::string node, bool merge_self_loop) {
    std::string node1, node2;
    if (graph[node].incoming_edges.size() == 2) {
        for (auto&& i : graph[node].incoming_edges) {
            if (i.first != node) {
                node1 = i.first;
                if (i.second.size() != 1)
                    return "";
            }
        }
    }
    else {
        for (auto&& i : graph[node].incoming_edges) {
            node1 = i.first;
            if (i.second.size() != 1)
                return "";
        }
    }
    if (this->graph[node].outgoing_edges.size() == 2) {
        for (auto&& j : graph[node].outgoing_edges) {
            if (j.first != node) {
                if (j.second.size() != 1)
                    return "";
                node2 = j.first;
            }
        }
    }
    else {
        for (auto&& j : graph[node].outgoing_edges) {
            if (j.second.size() != 1)
                return "";
            node2 = j.first;
        }
    }

    if (node1 == node && node == node2)
        return "";
    Path new_path;
    add_node_to_path(new_path, node1, 0, false);
    add_node_to_path(new_path, node, 0, false);
    add_node_to_path(new_path, node2, 0, false);
    if (merge_self_loop && graph[node].incoming_edges.size() == 2 && graph[node].outgoing_edges.size() == 2) {
        assert(graph[node].outgoing_edges.find(node) != graph[node].outgoing_edges.end());
        int times = std::round(graph[node].outgoing_edges[node].at(0).multiplicity / new_path.multiplicity);
        Path path_with_loop;
        add_node_to_path(path_with_loop, node1, 0, false);
        add_node_to_path(path_with_loop, node, 0, false);
        for (int i = 0; i < times; ++i) {
            add_node_to_path(path_with_loop, node, 0, false);
        }
        add_node_to_path(path_with_loop, node2, 0, false);
        new_path.length = path_with_loop.length;
        new_path.sequence = path_with_loop.sequence;
        new_path.path_edges_in_original_graph = path_with_loop.path_edges_in_original_graph;
        new_path.path_nodes_in_original_graph = path_with_loop.path_nodes_in_original_graph;
    }
    else if (merge_self_loop && node1 == node) {
        assert(graph[node].incoming_edges.size() == 1);
        int times = std::round(graph[node].outgoing_edges[node].at(0).multiplicity / graph[node].outgoing_edges[node2].at(0).multiplicity);
        Path path_with_loop;
        add_node_to_path(path_with_loop, node1, 0, false);
        for (int i = 0; i < times; ++i) {
            add_node_to_path(path_with_loop, node, 0, false);
        }
        add_node_to_path(path_with_loop, node2, 0, false);
        new_path.length = path_with_loop.length;
        new_path.sequence = path_with_loop.sequence;
        new_path.multiplicity = graph[node].outgoing_edges[node2].at(0).multiplicity;
        new_path.path_edges_in_original_graph = path_with_loop.path_edges_in_original_graph;
        new_path.path_nodes_in_original_graph = path_with_loop.path_nodes_in_original_graph;
    }
    else if (merge_self_loop && node2 == node) {
        assert(graph[node].outgoing_edges.size() == 1);
        int times = std::round(graph[node].outgoing_edges[node].at(0).multiplicity / graph[node1].outgoing_edges[node].at(0).multiplicity);
        Path path_with_loop;
        add_node_to_path(path_with_loop, node1, 0, false);
        add_node_to_path(path_with_loop, node2, 0, false);
        for (int i = 0; i < times; ++i) {
            add_node_to_path(path_with_loop, node, 0, false);
        }
        new_path.length = path_with_loop.length;
        new_path.sequence = path_with_loop.sequence;
        new_path.multiplicity = graph[node1].outgoing_edges[node].at(0).multiplicity;
        new_path.path_edges_in_original_graph = path_with_loop.path_edges_in_original_graph;
        new_path.path_nodes_in_original_graph = path_with_loop.path_nodes_in_original_graph;
    }

    Edge new_edge(new_path.sequence.at(k), new_path.length, new_path.sequence, new_path.multiplicity);
    new_edge.path_nodes_in_original_graph = new_path.path_nodes_in_original_graph;
    new_edge.path_edges_in_original_graph = new_path.path_edges_in_original_graph;
    assert(new_edge.path_nodes_in_original_graph.size() == new_edge.path_edges_in_original_graph.size() + 1);

    this->graph[node1].outgoing_edges.erase(node);
    this->graph[node].outgoing_edges.erase(node2);
    this->graph[node1].outgoing_edges[node2].push_back(new_edge);

    this->graph[node2].incoming_edges.erase(node);
    this->graph[node].incoming_edges.erase(node1);
    this->graph[node2].incoming_edges[node1].push_back(new_edge);
    if (node1 != node && node2 != node)
        this->graph.erase(node);
    // std::cout << "Merging non-branching " << node1 << "->" << node << "->" << node2 << " multi " << new_edge.multiplicity << std::endl;
    return new_edge.sequence;
}

void Graph::merge_non_branching_paths(bool merge_self_loop) {
    std::set<std::string> nodes_to_remove;
    for (auto&& node : this->graph) {
        if (nodes_to_remove.find(node.first) != nodes_to_remove.end() || nodes_to_remove.find(reverse_complementary_node(node.first)) != nodes_to_remove.end())
            continue;
        if (check_non_branching(node.first, merge_self_loop) && check_non_branching(reverse_complementary_node(node.first), merge_self_loop))
            nodes_to_remove.insert(node.first);
    }

    for (auto&& node : nodes_to_remove) {
        std::string seq1, seq2;
        if (check_non_branching(node, merge_self_loop))
            seq1 = merge_edges(node, merge_self_loop);
        if (check_non_branching(reverse_complementary_node(node), merge_self_loop))
            seq2 = merge_edges(reverse_complementary_node(node), merge_self_loop);
        // if (!seq1.empty() && !seq2.empty() && seq1 != reverse_complementary(seq2)) {
        //     std::cout << "Merging non-branching " << node << " and " << reverse_complementary_node(node) << " not reverse complementary (may exist self-reverse complementary edges)." << std::endl;
        //     // write_graph("debug_" + node);
        // }
    }
}

std::string Graph::collapse_bulge(std::string node1, std::string node2, unsigned& removed_bulges, double similarity) {
    auto& edges = this->graph[node1].outgoing_edges[node2];
    auto& edges_in = graph[node2].incoming_edges[node1];

    while (true) {
        int leg1 = -1, leg2 = -1;
        double sim = 0;
        for (int i = 0; i < edges.size() - 1; ++i) {
            for (int j = i + 1; j < edges.size(); ++j) {
                if (edges.at(j).length > edges.at(i).length)
                    sim = 1.0 * edges.at(i).length / edges.at(j).length;
                else
                    sim = 1.0 * edges.at(j).length / edges.at(i).length;
                if (sim >= similarity) {
                    leg1 = i;
                    leg2 = j;
                    break;
                }
            }
            if (leg1 != -1 && leg2 != -1)
                break;
        }
        if (leg1 == -1 || leg2 == -1)
            break;

        // keep the longer edge; if there is a tie, keep the edge with higher multiplicity
        if ((edges.at(leg2).length > edges.at(leg1).length) || (edges.at(leg2).length == edges.at(leg1).length && edges.at(leg2).multiplicity > edges.at(leg1).multiplicity)) {
            // edges.at(leg1).multiplicity = 1.0 * edges.at(leg1).multiplicity * edges.at(leg1).length / edges.at(leg2).length;
            edges.at(leg2).add_multi_from_edge_or_path(edges.at(leg1));
            edges.erase(edges.begin() + leg1);

            // edges_in.at(leg1).multiplicity = 1.0 * edges_in.at(leg1).multiplicity * edges_in.at(leg1).length / edges_in.at(leg2).length;
            edges_in.at(leg2).add_multi_from_edge_or_path(edges_in.at(leg1));
            edges_in.erase(edges_in.begin() + leg1);
        }
        else {
            // edges.at(leg2).multiplicity = 1.0 * edges.at(leg2).multiplicity * edges.at(leg2).length / edges.at(leg1).length;
            edges.at(leg1).add_multi_from_edge_or_path(edges.at(leg2));
            edges.erase(edges.begin() + leg2);

            // edges_in.at(leg2).multiplicity = 1.0 * edges_in.at(leg2).multiplicity * edges_in.at(leg2).length / edges_in.at(leg1).length;
            edges_in.at(leg1).add_multi_from_edge_or_path(edges_in.at(leg2));
            edges_in.erase(edges_in.begin() + leg2);
        }
        std::cout << "Simple bulge: " << node1 << "->" << node2 << ", sim (length) " << sim << " , muti " << graph[node1].outgoing_edges[node2].at(0).multiplicity << std::endl;
        removed_bulges += 1;
    }
    return edges.at(0).sequence;
}

void Graph::merge_tips(unsigned& num_tips, bool restrict_length, bool only_edges) {
    num_tips = 0;
    std::set<std::string> nodes_to_remove;
    for (auto&& node : graph) {
        if (nodes_to_remove.find(node.first) != nodes_to_remove.end())
            continue;
        std::vector<std::string> non_tips, tips;
        for (auto&& node_out : node.second.outgoing_edges) {
            if (graph[node_out.first].outgoing_edges.empty() && graph[node_out.first].incoming_edges.size() == 1) {
                if (node_out.second.size() == 1)
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
            std::string prefix_tip = graph[node.first].outgoing_edges[t].at(0).sequence.substr(k, 100000);
            double max_sim = 0;
            std::string max_edge;
            for (auto&& e : non_tips) {
                if (restrict_length == true && graph[node.first].outgoing_edges[t].at(0).length > graph[node.first].outgoing_edges[e].at(0).length && 1.0 * graph[node.first].outgoing_edges[e].at(0).length / graph[node.first].outgoing_edges[t].at(0).length < 0.9)
                    continue;
                else if (!restrict_length && graph[node.first].outgoing_edges[t].at(0).length > 10000000)
                    continue;

                // calculate similarity
                std::string prefix_edge = graph[node.first].outgoing_edges[e].at(0).sequence.substr(k, 100000);
                double sim = 1.0 * matches_by_edlib(prefix_tip, prefix_edge) / std::min(prefix_tip.size(), prefix_edge.size());
                if (sim > max_sim) {
                    max_sim = sim;
                    max_edge = e;
                }
            }
            if (max_edge.empty())
                continue;

            if (max_sim < 0.8)
                continue;

            if (restrict_length == true) {
                graph[node.first].outgoing_edges[max_edge].at(0).multiplicity += 1.0 * graph[node.first].outgoing_edges[t].at(0).multiplicity * graph[node.first].outgoing_edges[t].at(0).length / graph[node.first].outgoing_edges[max_edge].at(0).length;
                graph[max_edge].incoming_edges[node.first].at(0).multiplicity = graph[node.first].outgoing_edges[max_edge].at(0).multiplicity;

                graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).multiplicity +=
                    1.0 * graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(t)].at(0).multiplicity *
                    graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(t)].at(0).length /
                    graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).length;
                graph[reverse_complementary_node(max_edge)].outgoing_edges[reverse_complementary_node(node.first)].at(0).multiplicity = graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).multiplicity;
            }

            nodes_to_remove.insert(t);
            nodes_to_remove.insert(reverse_complementary_node(t));
            num_tips += 2;
            std::cout << "Tip " << node.first << "->" << t << " multi " << graph[node.first].outgoing_edges[t].at(0).multiplicity << " is merged to edge " << node.first << "->" << max_edge << " with sim " << max_sim << " multi " << graph[node.first].outgoing_edges[max_edge].at(0).multiplicity << std::endl;
            std::cout << "Tip " << reverse_complementary_node(t) << "->" << reverse_complementary_node(node.first) << " multi " << graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(t)].at(0).multiplicity << " is merged to edge " << reverse_complementary_node(max_edge) << "->" << reverse_complementary_node(node.first) << " with sim " << max_sim << " multi " << graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).multiplicity << std::endl;
            graph[node.first].outgoing_edges.erase(t);
            graph[t].incoming_edges.erase(node.first);
            graph[reverse_complementary_node(node.first)].incoming_edges.erase(reverse_complementary_node(t));
            graph[reverse_complementary_node(t)].outgoing_edges.erase(reverse_complementary_node(node.first));
        }
    }
    for (auto&& n : nodes_to_remove)
        graph.erase(n);

    if (only_edges == true) {
        unsigned bulges = 1;
        while (bulges) {
            multi_bulge_removal(bulges);
        }

        merge_non_branching_paths(true);
        return;
    }

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
            int max_index = 0;
            long max_length = 0;
            for (int i = 0; i < outgoing_tips.size();++i) {
                if (graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length > max_length) {
                    max_length = graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length;
                    max_index = i;
                }
            }
            for (int i = 0; i < outgoing_tips.size();++i) {
                if (i == max_index)
                    continue;
                std::string prefix_tip_target = graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).sequence.substr(k, 500000);
                std::string prefix_tip_to_merge = graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).sequence.substr(k, 500000);
                double sim = 1.0 * matches_by_edlib(prefix_tip_target, prefix_tip_to_merge) / std::min(prefix_tip_target.size(), prefix_tip_to_merge.size());
                if (sim < 0.8)
                    continue;
                if (sim < 0.9) {
                    double sim_len = 1.0 * std::min(graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length, graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length)
                        / std::max(graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length, graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length);
                    if (sim_len < 0.8)
                        continue;
                }
                std::cout << "Merge tip " << outgoing_tips.at(i) << " length " << graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length << " to tip " << outgoing_tips.at(max_index) << " length " << graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length << " sim (prefix) " << sim << std::endl;
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

// remove simple bulges, allow arbitrary multiplicity of legs
void Graph::multi_bulge_removal(unsigned& removed_bulges, bool skip_rc_bulges) {
    std::vector<std::string> nodes_to_remove;
    // for (auto&& node : this->graph) {
    //     std::vector<std::string> outgoing_tips;
    //     for (auto&& n : graph[node.first].outgoing_edges) {
    //         if (graph[n.first].outgoing_edges.size() == 0) {
    //             assert(graph[reverse_complementary_node(n.first)].incoming_edges.size() == 0);
    //             if (n.second.size() == 1)
    //                 outgoing_tips.push_back(n.first);
    //         }
    //     }
    //     if (outgoing_tips.size() >= 2) {
    //         // move edges of all tips to the first tip
    //         int max_index = 0;
    //         long max_length = 0;
    //         for (int i = 0; i < outgoing_tips.size();++i) {
    //             if (graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length > max_length) {
    //                 max_length = graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length;
    //                 max_index = i;
    //             }
    //         }
    //         for (int i = 0; i < outgoing_tips.size();++i) {
    //             if (i == max_index)
    //                 continue;
    //             std::string prefix_tip_target = graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).sequence.substr(k, 100000);
    //             std::string prefix_tip_to_merge = graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).sequence.substr(k, 100000);
    //             double sim = 1.0 * matches_by_edlib(prefix_tip_target, prefix_tip_to_merge) / std::min(prefix_tip_target.size(), prefix_tip_to_merge.size());
    //             if (sim < 0.8)
    //                 continue;
    //             std::cout << "Merge tip " << outgoing_tips.at(i) << " length " << graph[node.first].outgoing_edges[outgoing_tips[i]].at(0).length << " to tip " << outgoing_tips.at(max_index) << " length " << graph[node.first].outgoing_edges[outgoing_tips[max_index]].at(0).length << " sim (prefix) " << sim << std::endl;
    //             merge_vecs(graph[node.first].outgoing_edges[outgoing_tips[max_index]], graph[node.first].outgoing_edges[outgoing_tips[i]]);
    //             merge_vecs(graph[outgoing_tips[max_index]].incoming_edges[node.first], graph[outgoing_tips[i]].incoming_edges[node.first]);
    //             merge_vecs(graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(outgoing_tips[max_index])], graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(outgoing_tips[i])]);
    //             merge_vecs(graph[reverse_complementary_node(outgoing_tips[max_index])].outgoing_edges[reverse_complementary_node(node.first)], graph[reverse_complementary_node(outgoing_tips[i])].outgoing_edges[reverse_complementary_node(node.first)]);
    //             graph[node.first].outgoing_edges.erase(outgoing_tips[i]);
    //             graph[outgoing_tips[i]].incoming_edges.erase(node.first);
    //             graph[reverse_complementary_node(node.first)].incoming_edges.erase(reverse_complementary_node(outgoing_tips[i]));
    //             graph[reverse_complementary_node(outgoing_tips[i])].outgoing_edges.erase(reverse_complementary_node(node.first));
    //             if (graph[outgoing_tips[i]].incoming_edges.size() == 0) {
    //                 assert(graph[reverse_complementary_node(outgoing_tips[i])].outgoing_edges.size() == 0);
    //                 nodes_to_remove.push_back(outgoing_tips[i]);
    //                 nodes_to_remove.push_back(reverse_complementary_node(outgoing_tips[i]));
    //             }
    //         }
    //     }
    // }
    // for (auto&& n : nodes_to_remove) {
    //     this->graph.erase(n);
    // }
    removed_bulges = 0;
    for (auto&& node : this->graph) {
        for (auto&& i : node.second.outgoing_edges) {
            if (i.second.size() >= 2) {
                std::string node1 = node.first;
                std::string node2 = i.first;

                if (skip_rc_bulges) {
                    if (node1 == reverse_complementary_node(node2))
                        continue;
                }
                // collapse forward bulge
                std::string seq1 = this->collapse_bulge(node1, node2, removed_bulges);

                // if node1 and node2 are reverse complementary nodes, skip collapsing reverse bulge
                // if (reverse_complementary_node(node2) != node1) {
                //     std::string seq2 = this->collapse_bulge(reverse_complementary_node(node2), reverse_complementary_node(node1), removed_bulges);
                //     if (seq1 != reverse_complementary(seq2))
                //         std::cout << "Problematic bulge collapsing: resulting sequences not reverse complementary" << std::endl;
                // }
            }
        }
    }
    this->merge_non_branching_paths();
}

void Graph::gluing_broken_bulges(unsigned& removed_bulges) {
    std::vector<std::string> nodes_to_remove;
    for (auto&& node : this->graph) {
        for (auto&& node2 : node.second.outgoing_edges) {
            if (node2.second.size() >= 2)
                continue;
            // find the first incoming tip in the end node
            std::string incoming_tip;
            Edge incoming_edge;
            for (auto&& in_node : graph[node2.first].incoming_edges) {
                if (in_node.second.size() == 1 && graph[in_node.first].incoming_edges.empty()) {
                    incoming_tip = in_node.first;
                    incoming_edge = in_node.second.at(0);
                }
            }
            // find the first outgoing tip in the start node
            std::string outgoing_tip;
            Edge outgoing_edge;
            for (auto&& out_node : graph[node.first].outgoing_edges) {
                if (out_node.second.size() == 1 && graph[out_node.first].outgoing_edges.empty()) {
                    outgoing_tip = out_node.first;
                    outgoing_edge = out_node.second.at(0);
                }
            }

            if (!incoming_tip.empty() && !outgoing_tip.empty() && outgoing_tip != node2.first && incoming_tip != node.first) {
                std::string seq_broken = outgoing_edge.sequence + incoming_edge.sequence.substr(k);
                // double similarity = 1.0 * matches_by_edlib(seq_broken, node2.second.at(0).sequence) / std::max(seq_broken.size(), node2.second.at(0).sequence.size());
                double similarity = 1.0 * std::min(outgoing_edge.length + incoming_edge.length, node2.second.at(0).length) / std::max(outgoing_edge.length + incoming_edge.length, node2.second.at(0).length);
                if (similarity < 0.9)
                    continue;
                std::cout << "Broken bulge (and reverse complementary): " << node.first << "->" << node2.first << ", " << node.first << "->" << outgoing_tip << " " << incoming_tip << "->" << node2.first << ", sim " << similarity << std::endl;
                std::vector<std::string> nodes;
                std::string out_tip_seq = outgoing_edge.sequence.substr(outgoing_edge.sequence.size() - k);
                for (auto&& n : graph[incoming_tip].outgoing_edges) {
                    nodes.push_back(n.first);
                    for (auto&& e : n.second) {
                        e.sequence = out_tip_seq + e.sequence.substr(k);
                    }
                    for (auto&& e : graph[n.first].incoming_edges[incoming_tip]) {
                        e.sequence = out_tip_seq + e.sequence.substr(k);
                    }
                }
                for (auto&& n : graph[reverse_complementary_node(incoming_tip)].incoming_edges) {
                    for (auto&& e : n.second) {
                        e.sequence = e.sequence.substr(0, e.sequence.size() - k) + reverse_complementary(out_tip_seq);
                    }
                    for (auto&& e : graph[n.first].outgoing_edges[reverse_complementary_node(incoming_tip)]) {
                        e.sequence = e.sequence.substr(0, e.sequence.size() - k) + reverse_complementary(out_tip_seq);
                    }
                }
                for (auto&& n : nodes) {
                    merge_vecs(graph[n].incoming_edges[outgoing_tip], graph[n].incoming_edges[incoming_tip]);
                    merge_vecs(graph[outgoing_tip].outgoing_edges[n], graph[incoming_tip].outgoing_edges[n]);
                    graph[n].incoming_edges.erase(incoming_tip);
                    graph[incoming_tip].outgoing_edges.erase(n);

                    merge_vecs(graph[reverse_complementary_node(n)].outgoing_edges[reverse_complementary_node(outgoing_tip)], graph[reverse_complementary_node(n)].outgoing_edges[reverse_complementary_node(incoming_tip)]);
                    merge_vecs(graph[reverse_complementary_node(outgoing_tip)].incoming_edges[reverse_complementary_node(n)], graph[reverse_complementary_node(incoming_tip)].incoming_edges[reverse_complementary_node(n)]);
                    graph[reverse_complementary_node(n)].outgoing_edges.erase(reverse_complementary_node(incoming_tip));
                    graph[reverse_complementary_node(incoming_tip)].incoming_edges.erase(reverse_complementary_node(n));
                }
                assert(graph[incoming_tip].incoming_edges.empty() && graph[incoming_tip].outgoing_edges.empty());
                assert(graph[reverse_complementary_node(incoming_tip)].incoming_edges.empty() && graph[reverse_complementary_node(incoming_tip)].outgoing_edges.empty());
                nodes_to_remove.push_back(incoming_tip);
                nodes_to_remove.push_back(reverse_complementary_node(incoming_tip));
            }
        }
    }
    for (auto&& n : nodes_to_remove)
        this->graph.erase(n);
    this->merge_non_branching_paths();
    this->multi_bulge_removal(removed_bulges);
}

void Graph::merge_tips_into_edges(unsigned& num_tips) {
    unsigned removed_paths = 1;
    while (removed_paths) {
        resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 2);
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
                if (node_out.second.size() == 1)
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
            std::string prefix_tip = graph[node.first].outgoing_edges[t].at(0).sequence.substr(k, 100000);
            double max_sim = 0;
            std::string max_edge;
            for (auto&& e : non_tips) {
                // the tip should be short: at least shorter than the edge
                // if (graph[node.first].outgoing_edges[t].at(0).length > graph[node.first].outgoing_edges[e].at(0).length)
                if (graph[node.first].outgoing_edges[t].at(0).length > graph[node.first].outgoing_edges[e].at(0).length && 1.0 * graph[node.first].outgoing_edges[e].at(0).length / graph[node.first].outgoing_edges[t].at(0).length < 0.9)
                    continue;

                if (graph[node.first].outgoing_edges[t].at(0).multiplicity > graph[node.first].outgoing_edges[e].at(0).multiplicity * 0.8)
                    continue;

                // calculate similarity
                std::string prefix_edge = graph[node.first].outgoing_edges[e].at(0).sequence.substr(k, 100000);
                double sim = 1.0 * matches_by_edlib(prefix_tip, prefix_edge) / std::min(prefix_tip.size(), prefix_edge.size());
                if (sim > max_sim) {
                    max_sim = sim;
                    max_edge = e;
                }
            }
            if (max_edge.empty())
                continue;

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
            std::cout << "Tip " << node.first << "->" << t << " multi " << graph[node.first].outgoing_edges[t].at(0).multiplicity << " is merged to edge " << node.first << "->" << max_edge << " with sim " << max_sim << " multi " << graph[node.first].outgoing_edges[max_edge].at(0).multiplicity << std::endl;
            std::cout << "Tip " << reverse_complementary_node(t) << "->" << reverse_complementary_node(node.first) << " multi " << graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(t)].at(0).multiplicity << " is merged to edge " << reverse_complementary_node(max_edge) << "->" << reverse_complementary_node(node.first) << " with sim " << max_sim << " multi " << graph[reverse_complementary_node(node.first)].incoming_edges[reverse_complementary_node(max_edge)].at(0).multiplicity << std::endl;
            graph[node.first].outgoing_edges.erase(t);
            graph[t].incoming_edges.erase(node.first);
            graph[reverse_complementary_node(node.first)].incoming_edges.erase(reverse_complementary_node(t));
            graph[reverse_complementary_node(t)].outgoing_edges.erase(reverse_complementary_node(node.first));
        }
    }
    for (auto&& n : nodes_to_remove)
        graph.erase(n);
    merge_non_branching_paths(true);
}

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

void Graph::resolve_2_in_2_out_mdbg(std::string multidbg, std::string output) {
    unsigned removed_bulges = 1;
    while (removed_bulges) {
        this->multi_bulge_removal(removed_bulges);
    }
    std::vector<std::string> nodes_to_remove;
    std::vector<std::string> source_nodes, sink_nodes;
    std::unordered_set<std::string> retained_nodes;
    //search for 2-in-2-out edge
    for (auto&& node : graph) {
        if (node.second.incoming_edges.size() >= 2 && node.second.outgoing_edges.size() == 1 && node.second.incoming_edges.find(node.first) == node.second.incoming_edges.end()) {
            std::string node_sink = node.first;
            bool flag = false;
            for (auto&& n : graph[node_sink].outgoing_edges) {
                if (n.first != node_sink) {
                    node_sink = n.first;
                    flag = true;
                }
            }
            if (flag == false)
                continue;

            // find a 2-in-2-out component
            if (this->graph[node_sink].incoming_edges.size() == 1 && this->graph[node_sink].outgoing_edges.size() >= 2 && graph[node_sink].outgoing_edges.find(node_sink) == graph[node_sink].outgoing_edges.end() && std::find(source_nodes.begin(), source_nodes.end(), node.first) == source_nodes.end()) {
                source_nodes.push_back(node.first);
                sink_nodes.push_back(node_sink);
                source_nodes.push_back(reverse_complementary_node(node_sink));
                sink_nodes.push_back(reverse_complementary_node(node.first));
                retained_nodes.insert(node.first);
                retained_nodes.insert(reverse_complementary_node(node.first));
                retained_nodes.insert(node_sink);
                retained_nodes.insert(reverse_complementary_node(node_sink));
            }
        }
        else if (node.second.incoming_edges.size() == 2 && node.second.outgoing_edges.size() == 2 && node.second.outgoing_edges.find(node.first) == node.second.outgoing_edges.end() && std::find(source_nodes.begin(), source_nodes.end(), node.first) == source_nodes.end()) {
            source_nodes.push_back(node.first);
            sink_nodes.push_back(node.first);
            source_nodes.push_back(reverse_complementary_node(node.first));
            sink_nodes.push_back(reverse_complementary_node(node.first));
            retained_nodes.insert(node.first);
            retained_nodes.insert(reverse_complementary_node(node.first));
        }
    }

    std::cout << "Nodes in 2-in-2-out: " << retained_nodes.size() << " " << source_nodes.size() << std::endl;

    write_graph(output + ".2in_2out", 10000, false, false, retained_nodes);
    // write_graph(output + ".2in_2out", 10000, false, false);
    // if (system(("get_prefix_suffix.py " + output + ".2in_2out.fasta " + output + ".2in_2out.pre.fasta").c_str()) != 0) {
    //     exit(1);
    // }
    if (system(("minimap2 -ax map-pb " + multidbg + " " + output + ".2in_2out.fasta -t 100 | grep -v '^@' > " + output + ".2in_2out.sam").c_str()) != 0) {
        exit(1);
    }
    if (!std::filesystem::exists(multidbg + ".fai")) {
        if (system(("samtools faidx " + multidbg).c_str()) != 0)
            exit(1);
    }
    if (system(("cut -f1,2 " + multidbg + ".fai | awk " + R"('{print "@SQ\tSN:"$1"\tLN:"$2}')" + " > " + output + ".2in_2out.header.sam").c_str()) != 0)
        exit(1);
    if (system(("cat " + output + ".2in_2out.header.sam " + output + ".2in_2out.sam | samtools sort -@ 50 -o " + output + ".2in_2out.bam").c_str()) != 0) {
        exit(1);
    }
    // if (system(("minimap2 -ax map-pb " + multidbg + " " + output + ".2in_2out.fasta -t 100 | samtools sort -@ 50 -o " + output + ".2in_2out.bam").c_str()) != 0) {
    //     exit(1);
    // }
    std::string exeDir = getExecutablePath();
    if (system((exeDir + "/get_reference.py -o " + output + ".2in_2out.bam.stats " + output + ".2in_2out.bam " + output + ".2in_2out.fasta").c_str()) != 0)
        exit(1);
    write_graph_colored_from_bam(output + ".2in_2out.multi", output + ".2in_2out.bam.stats");

    for (int i = 0; i < source_nodes.size();++i)
        this->resolve_2_in_2_out(source_nodes[i], sink_nodes[i], nodes_to_remove);

    for (auto&& node : nodes_to_remove) {
        this->graph.erase(node);
    }
    nodes_to_remove.clear();
    this->merge_non_branching_paths(true);
}

void Graph::resolve_2_in_2_out(std::string node1, std::string node2, std::vector<std::string>& nodes_to_remove) {
    std::vector<std::string> incoming_nodes, outgoing_nodes;

    for (auto&& e : this->graph[node1].incoming_edges) {
        incoming_nodes.push_back(e.first);
        if (e.second.size() != 1)
            return;
    }
    for (auto&& e : this->graph[node2].outgoing_edges) {
        outgoing_nodes.push_back(e.first);
        if (e.second.size() != 1)
            return;
    }

    if (node1 == node2) {
        std::unordered_set<std::string> tmp_ids, supporting11, supporting12, supporting21, supporting22;
        // check incoming 1
        tmp_ids.clear();
        for (auto&& i : graph[incoming_nodes[0]].outgoing_edges[node1].at(0).ref_ids) {
            tmp_ids.insert(i.substr(0, i.find(' ')));
        }
        // out 1
        for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[0]].at(0).ref_ids) {
            std::string edge_id = i.substr(0, i.find(' '));
            if (tmp_ids.find(edge_id) != tmp_ids.end())
                supporting11.insert(edge_id);
        }
        // out 2
        for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[1]].at(0).ref_ids) {
            std::string edge_id = i.substr(0, i.find(' '));
            if (tmp_ids.find(edge_id) != tmp_ids.end())
                supporting12.insert(edge_id);
        }

        // check incoming 2
        tmp_ids.clear();
        for (auto&& i : graph[incoming_nodes[1]].outgoing_edges[node1].at(0).ref_ids) {
            tmp_ids.insert(i.substr(0, i.find(' ')));
        }
        // out 1
        for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[0]].at(0).ref_ids) {
            std::string edge_id = i.substr(0, i.find(' '));
            if (tmp_ids.find(edge_id) != tmp_ids.end())
                supporting21.insert(edge_id);
        }
        // out 2
        for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[1]].at(0).ref_ids) {
            std::string edge_id = i.substr(0, i.find(' '));
            if (tmp_ids.find(edge_id) != tmp_ids.end())
                supporting22.insert(edge_id);
        }

        // check confict
        bool flag11 = false, flag12 = false;
        if (supporting11.size() != 0 && supporting12.size() == 0 && supporting21.size() == 0) {
            flag11 = true;
        }
        if (supporting22.size() != 0 && supporting21.size() == 0 && supporting12.size() == 0) {
            flag11 = true;
        }
        if (supporting12.size() != 0 && supporting11.size() == 0 && supporting22.size() == 0) {
            flag12 = true;
        }
        if (supporting21.size() != 0 && supporting22.size() == 0 && supporting11.size() == 0) {
            flag12 = true;
        }
        assert(!(flag11 && flag12));

        int id1, id2, id3, id4;
        if (flag11) {
            id1 = 0, id2 = 0, id3 = 1, id4 = 1;
        }
        else {
            id1 = 0, id2 = 1, id3 = 1, id4 = 0;
        }

        if (!flag11 && !flag12)
            return;
        if (std::abs(graph[incoming_nodes[id1]].outgoing_edges[node1].at(0).multiplicity - graph[node2].outgoing_edges[outgoing_nodes[id2]].at(0).multiplicity) / std::min(graph[incoming_nodes[id1]].outgoing_edges[node1].at(0).multiplicity, graph[node2].outgoing_edges[outgoing_nodes[id2]].at(0).multiplicity) > 0.5)
            return;
        if (std::abs(graph[incoming_nodes[id3]].outgoing_edges[node1].at(0).multiplicity - graph[node2].outgoing_edges[outgoing_nodes[id4]].at(0).multiplicity) / std::min(graph[incoming_nodes[id3]].outgoing_edges[node1].at(0).multiplicity, graph[node2].outgoing_edges[outgoing_nodes[id4]].at(0).multiplicity) > 0.5)
            return;
        if (flag11) {
            std::cout << "Resolving 2-in-2-out edge: path1: " << incoming_nodes[0] << " -> " << node2 << " -> " << outgoing_nodes[0] << ", #supporting ids " << supporting11.size() << "; path2: " << incoming_nodes[1] << "->" << node2 << " -> " << outgoing_nodes[1] << ", #supporting ids " << supporting22.size() << std::endl;
        }
        else {
            std::cout << "Resolving 2-in-2-out edge: path1: " << incoming_nodes[0] << " -> " << node2 << " -> " << outgoing_nodes[1] << ", #supporting ids " << supporting12.size() << "; path2: " << incoming_nodes[1] << "->" << node2 << " -> " << outgoing_nodes[0] << ", #supporting ids " << supporting21.size() << std::endl;
        }
        Path path1, path2;
        this->add_node_to_path(path1, incoming_nodes[id1]);
        this->add_node_to_path(path1, node1);
        this->add_node_to_path(path1, outgoing_nodes[id2]);
        this->add_node_to_path(path2, incoming_nodes[id3]);
        this->add_node_to_path(path2, node1);
        this->add_node_to_path(path2, outgoing_nodes[id4]);

        // get fraction of coverage for the middile edge
        auto edge1 = graph[incoming_nodes[id1]].outgoing_edges[node1].at(0), edge2 = graph[node2].outgoing_edges[outgoing_nodes[id2]].at(0);
        auto edge3 = graph[incoming_nodes[id3]].outgoing_edges[node1].at(0), edge4 = graph[node2].outgoing_edges[outgoing_nodes[id4]].at(0);
        double coverage1 = 1.0 * (edge1.multiplicity * edge1.length + edge2.multiplicity * edge2.length) / path1.length;
        double coverage2 = 1.0 * (edge3.multiplicity * edge3.length + edge4.multiplicity * edge4.length) / path2.length;

        path1.multiplicity = coverage1;
        path2.multiplicity = coverage2;

        Edge edge_final1(path1.sequence.at(k), path1.length, path1.sequence, path1.multiplicity);
        edge_final1.path_edges_in_original_graph = path1.path_edges_in_original_graph;
        edge_final1.path_nodes_in_original_graph = path1.path_nodes_in_original_graph;
        this->graph[incoming_nodes.at(id1)].outgoing_edges[outgoing_nodes.at(id2)].push_back(edge_final1);
        this->graph[outgoing_nodes.at(id2)].incoming_edges[incoming_nodes.at(id1)].push_back(edge_final1);

        Edge edge_final2(path2.sequence.at(k), path2.length, path2.sequence, path2.multiplicity);
        edge_final2.path_edges_in_original_graph = path2.path_edges_in_original_graph;
        edge_final2.path_nodes_in_original_graph = path2.path_nodes_in_original_graph;
        this->graph[incoming_nodes.at(id3)].outgoing_edges[outgoing_nodes.at(id4)].push_back(edge_final2);
        this->graph[outgoing_nodes.at(id4)].incoming_edges[incoming_nodes.at(id3)].push_back(edge_final2);

        for (int i = 0; i < path1.nodes.size() - 1; ++i) {
            this->graph[path1.nodes[i]].outgoing_edges.erase(path1.nodes[i + 1]);
            this->graph[path1.nodes[i + 1]].incoming_edges.erase(path1.nodes[i]);
        }

        for (int i = 0; i < path2.nodes.size() - 1; ++i) {
            this->graph[path2.nodes[i]].outgoing_edges.erase(path2.nodes[i + 1]);
            this->graph[path2.nodes[i + 1]].incoming_edges.erase(path2.nodes[i]);
        }

        nodes_to_remove.push_back(node1);
    }
    else {
        if (incoming_nodes.size() == 2 && outgoing_nodes.size() == 2) {
            std::unordered_set<std::string> supporting_ids, tmp_ids, supporting11, supporting12, supporting21, supporting22;
            // check incoming 1
            supporting_ids.clear();
            tmp_ids.clear();
            for (auto&& i : graph[incoming_nodes[0]].outgoing_edges[node1].at(0).ref_ids) {
                supporting_ids.insert(i.substr(0, i.find(' ')));
            }
            for (auto&& i : graph[node1].outgoing_edges[node2].at(0).ref_ids) {
                std::string edge_id = i.substr(0, i.find(' '));
                if (supporting_ids.find(edge_id) != supporting_ids.end())
                    tmp_ids.insert(edge_id);
            }
            // out 1
            for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[0]].at(0).ref_ids) {
                std::string edge_id = i.substr(0, i.find(' '));
                if (tmp_ids.find(edge_id) != tmp_ids.end())
                    supporting11.insert(edge_id);
            }
            // out 2
            for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[1]].at(0).ref_ids) {
                std::string edge_id = i.substr(0, i.find(' '));
                if (tmp_ids.find(edge_id) != tmp_ids.end())
                    supporting12.insert(edge_id);
            }

            // check incoming 2
            supporting_ids.clear();
            tmp_ids.clear();
            for (auto&& i : graph[incoming_nodes[1]].outgoing_edges[node1].at(0).ref_ids) {
                supporting_ids.insert(i.substr(0, i.find(' ')));
            }
            for (auto&& i : graph[node1].outgoing_edges[node2].at(0).ref_ids) {
                std::string edge_id = i.substr(0, i.find(' '));
                if (supporting_ids.find(edge_id) != supporting_ids.end())
                    tmp_ids.insert(edge_id);
            }
            // out 1
            for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[0]].at(0).ref_ids) {
                std::string edge_id = i.substr(0, i.find(' '));
                if (tmp_ids.find(edge_id) != tmp_ids.end())
                    supporting21.insert(edge_id);
            }
            // out 2
            for (auto&& i : graph[node2].outgoing_edges[outgoing_nodes[1]].at(0).ref_ids) {
                std::string edge_id = i.substr(0, i.find(' '));
                if (tmp_ids.find(edge_id) != tmp_ids.end())
                    supporting22.insert(edge_id);
            }

            // check confict
            bool flag11 = false, flag12 = false;
            if (supporting11.size() != 0 && supporting12.size() == 0 && supporting21.size() == 0) {
                flag11 = true;
            }
            if (supporting22.size() != 0 && supporting21.size() == 0 && supporting12.size() == 0) {
                flag11 = true;
            }
            if (supporting12.size() != 0 && supporting11.size() == 0 && supporting22.size() == 0) {
                flag12 = true;
            }
            if (supporting21.size() != 0 && supporting22.size() == 0 && supporting11.size() == 0) {
                flag12 = true;
            }
            assert(!(flag11 && flag12));

            int id1, id2, id3, id4;
            if (flag11) {
                id1 = 0, id2 = 0, id3 = 1, id4 = 1;
            }
            else {
                id1 = 0, id2 = 1, id3 = 1, id4 = 0;
            }

            if (!flag11 && !flag12)
                return;
            if (std::abs(graph[incoming_nodes[id1]].outgoing_edges[node1].at(0).multiplicity - graph[node2].outgoing_edges[outgoing_nodes[id2]].at(0).multiplicity) / std::min(graph[incoming_nodes[id1]].outgoing_edges[node1].at(0).multiplicity, graph[node2].outgoing_edges[outgoing_nodes[id2]].at(0).multiplicity) > 0.5)
                return;
            if (std::abs(graph[incoming_nodes[id3]].outgoing_edges[node1].at(0).multiplicity - graph[node2].outgoing_edges[outgoing_nodes[id4]].at(0).multiplicity) / std::min(graph[incoming_nodes[id3]].outgoing_edges[node1].at(0).multiplicity, graph[node2].outgoing_edges[outgoing_nodes[id4]].at(0).multiplicity) > 0.5)
                return;
            if (flag11) {
                std::cout << "Resolving 2-in-2-out edge: path1: " << incoming_nodes[0] << " -> " << node1 << " -> " << node2 << " -> " << outgoing_nodes[0] << ", #supporting ids " << supporting11.size() << "; path2: " << incoming_nodes[1] << " -> " << node1 << "->" << node2 << " -> " << outgoing_nodes[1] << ", #supporting ids " << supporting22.size() << std::endl;
            }
            else {
                std::cout << "Resolving 2-in-2-out edge: path1: " << incoming_nodes[0] << " -> " << node1 << " -> " << node2 << " -> " << outgoing_nodes[1] << ", #supporting ids " << supporting12.size() << "; path2: " << incoming_nodes[1] << " -> " << node1 << "->" << node2 << " -> " << outgoing_nodes[0] << ", #supporting ids " << supporting21.size() << std::endl;
            }
            Path path1, path2;
            this->add_node_to_path(path1, incoming_nodes[id1]);
            this->add_node_to_path(path1, node1);
            this->add_node_to_path(path1, node2);
            this->add_node_to_path(path1, outgoing_nodes[id2]);
            this->add_node_to_path(path2, incoming_nodes[id3]);
            this->add_node_to_path(path2, node1);
            this->add_node_to_path(path2, node2);
            this->add_node_to_path(path2, outgoing_nodes[id4]);

            // get fraction of coverage for the middile edge
            auto edge1 = graph[incoming_nodes[id1]].outgoing_edges[node1].at(0), edge2 = graph[node2].outgoing_edges[outgoing_nodes[id2]].at(0);
            auto edge3 = graph[incoming_nodes[id3]].outgoing_edges[node1].at(0), edge4 = graph[node2].outgoing_edges[outgoing_nodes[id4]].at(0);
            double coverage1 = 1.0 * (edge1.multiplicity * edge1.length + edge2.multiplicity * edge2.length) / path1.length;
            double coverage2 = 1.0 * (edge3.multiplicity * edge3.length + edge4.multiplicity * edge4.length) / path2.length;

            path1.multiplicity = coverage1 + graph[node1].outgoing_edges[node2].at(0).multiplicity * graph[node1].outgoing_edges[node2].at(0).length * coverage1 / (coverage1 + coverage2) / path1.length;
            path2.multiplicity = coverage2 + graph[node1].outgoing_edges[node2].at(0).multiplicity * graph[node1].outgoing_edges[node2].at(0).length * coverage2 / (coverage1 + coverage2) / path2.length;

            Edge edge_final1(path1.sequence.at(k), path1.length, path1.sequence, path1.multiplicity);
            edge_final1.path_edges_in_original_graph = path1.path_edges_in_original_graph;
            edge_final1.path_nodes_in_original_graph = path1.path_nodes_in_original_graph;
            this->graph[incoming_nodes.at(id1)].outgoing_edges[outgoing_nodes.at(id2)].push_back(edge_final1);
            this->graph[outgoing_nodes.at(id2)].incoming_edges[incoming_nodes.at(id1)].push_back(edge_final1);

            Edge edge_final2(path2.sequence.at(k), path2.length, path2.sequence, path2.multiplicity);
            edge_final2.path_edges_in_original_graph = path2.path_edges_in_original_graph;
            edge_final2.path_nodes_in_original_graph = path2.path_nodes_in_original_graph;
            this->graph[incoming_nodes.at(id3)].outgoing_edges[outgoing_nodes.at(id4)].push_back(edge_final2);
            this->graph[outgoing_nodes.at(id4)].incoming_edges[incoming_nodes.at(id3)].push_back(edge_final2);

            for (int i = 0; i < path1.nodes.size() - 1; ++i) {
                this->graph[path1.nodes[i]].outgoing_edges.erase(path1.nodes[i + 1]);
                this->graph[path1.nodes[i + 1]].incoming_edges.erase(path1.nodes[i]);
            }

            for (int i = 0; i < path2.nodes.size() - 1; ++i) {
                this->graph[path2.nodes[i]].outgoing_edges.erase(path2.nodes[i + 1]);
                this->graph[path2.nodes[i + 1]].incoming_edges.erase(path2.nodes[i]);
            }

            nodes_to_remove.push_back(node1);
            nodes_to_remove.push_back(node2);
        }
        else {
            // from incoming to outgoing
            std::unordered_map<int, std::unordered_map<int, std::unordered_set< std::string>>> supporting_ids_in_to_out;
            std::vector<double> in_multis, out_multis;
            // deduplicate
            std::vector<std::vector<std::string>> incoming_refs, outgoing_refs;
            std::unordered_set<std::string>all_refs, all_duplicate;
            incoming_refs.resize(incoming_nodes.size());
            outgoing_refs.resize(outgoing_nodes.size());

            for (int i = 0; i < incoming_nodes.size(); ++i) {
                for (auto&& r : graph[incoming_nodes[i]].outgoing_edges[node1].at(0).ref_ids) {
                    if (all_refs.find(r.substr(0, r.find(' '))) != all_refs.end())
                        all_duplicate.insert(r.substr(0, r.find(' ')));
                    all_refs.insert(r.substr(0, r.find(' ')));
                }
            }
            for (int i = 0; i < incoming_nodes.size(); ++i) {
                for (auto&& r : graph[incoming_nodes[i]].outgoing_edges[node1].at(0).ref_ids) {
                    if (all_duplicate.find(r.substr(0, r.find(' '))) == all_duplicate.end())
                        incoming_refs[i].push_back(r.substr(0, r.find(' ')));
                }
            }

            all_refs.clear(), all_duplicate.clear();
            for (int i = 0; i < outgoing_nodes.size(); ++i) {
                for (auto&& r : graph[node2].outgoing_edges[outgoing_nodes[i]].at(0).ref_ids) {
                    if (all_refs.find(r.substr(0, r.find(' '))) != all_refs.end())
                        all_duplicate.insert(r.substr(0, r.find(' ')));
                    all_refs.insert(r.substr(0, r.find(' ')));
                }
            }
            for (int i = 0; i < outgoing_nodes.size(); ++i) {
                for (auto&& r : graph[node2].outgoing_edges[outgoing_nodes[i]].at(0).ref_ids) {
                    if (all_duplicate.find(r.substr(0, r.find(' '))) == all_duplicate.end())
                        outgoing_refs[i].push_back(r.substr(0, r.find(' ')));
                }
            }

            bool flag = false;
            for (int i = 0; i < incoming_nodes.size(); ++i) {
                in_multis.push_back(graph[incoming_nodes[i]].outgoing_edges[node1].at(0).multiplicity);
                std::unordered_set<std::string> supporting_ids, tmp_ids;
                for (auto&& r : incoming_refs[i]) {
                    supporting_ids.insert(r);
                }
                for (auto&& r : graph[node1].outgoing_edges[node2].at(0).ref_ids) {
                    std::string edge_id = r.substr(0, r.find(' '));
                    if (supporting_ids.find(edge_id) != supporting_ids.end())
                        tmp_ids.insert(edge_id);
                }
                for (int j = 0; j < outgoing_nodes.size(); ++j) {
                    for (auto&& r : outgoing_refs[j]) {
                        if (tmp_ids.find(r) != tmp_ids.end()) {
                            supporting_ids_in_to_out[i][j].insert(r);
                            flag = true;
                        }
                    }
                }
            }
            if (!flag)
                return;
            for (int j = 0; j < outgoing_nodes.size(); ++j) {
                out_multis.push_back(graph[node2].outgoing_edges[outgoing_nodes[j]].at(0).multiplicity);
            }
            for (int i = 0; i < incoming_nodes.size(); ++i) {
                for (int j = 0; j < outgoing_nodes.size(); ++j) {
                    if (supporting_ids_in_to_out[i][j].size() != 0) {
                        if (in_multis[i] == 0 || out_multis[j] == 0)
                            return;
                        if (in_multis[i] >= out_multis[j]) {
                            if (in_multis[i] / out_multis[j] <= 1.5) {
                                in_multis[i] = 0;
                                out_multis[j] = 0;
                            }
                            else {
                                in_multis[i] -= out_multis[j];
                                out_multis[j] = 0;
                            }
                        }
                        else if (out_multis[j] / in_multis[i] <= 1.5) {
                            in_multis[i] = 0;
                            out_multis[j] = 0;
                        }
                        else {
                            out_multis[j] -= in_multis[i];
                            in_multis[i] = 0;
                        }
                    }
                }
            }

            double in_multi = 0, out_multi = 0;
            for (auto&& m : in_multis)
                in_multi += m;
            for (auto&& m : out_multis)
                out_multi += m;

            if (in_multi + out_multi != 0) {
                if (in_multi >= out_multi && out_multi / in_multi < 0.5)
                    return;
                if (in_multi <= out_multi && in_multi / out_multi < 0.5)
                    return;
            }

            std::cout << "Resolving m-in-n-out edge: " << "... -> " << node1 << " -> " << node2 << " -> ..." << std::endl;

            for (int i = 0; i < incoming_nodes.size(); ++i) {
                for (int j = 0; j < outgoing_nodes.size(); ++j) {
                    if (supporting_ids_in_to_out[i][j].size() != 0) {
                        Path path;
                        this->add_node_to_path(path, incoming_nodes[i]);
                        this->add_node_to_path(path, node1);
                        this->add_node_to_path(path, node2);
                        this->add_node_to_path(path, outgoing_nodes[j]);

                        auto edge1 = graph[incoming_nodes[i]].outgoing_edges[node1].at(0), edge2 = graph[node2].outgoing_edges[outgoing_nodes[j]].at(0);
                        path.multiplicity = 1.0 * (edge1.multiplicity * edge1.length + edge2.multiplicity * edge2.length) / (edge1.length + edge2.length);
                        std::cout << "    " << incoming_nodes[i] << " -> " << node1 << " -> " << node2 << " -> " << outgoing_nodes[j] << ", #supporting ids " << supporting_ids_in_to_out[i][j].size() << " multi " << path.multiplicity << std::endl;

                        Edge edge_final(path.sequence.at(k), path.length, path.sequence, path.multiplicity);
                        edge_final.path_edges_in_original_graph = path.path_edges_in_original_graph;
                        edge_final.path_nodes_in_original_graph = path.path_nodes_in_original_graph;
                        this->graph[incoming_nodes.at(i)].outgoing_edges[outgoing_nodes.at(j)].push_back(edge_final);
                        this->graph[outgoing_nodes.at(j)].incoming_edges[incoming_nodes.at(i)].push_back(edge_final);
                    }
                }
            }

            for (int i = 0; i < in_multis.size(); ++i) {
                if (in_multis[i] == 0) {
                    graph[incoming_nodes.at(i)].outgoing_edges.erase(node1);
                    graph[node1].incoming_edges.erase(incoming_nodes.at(i));
                }
                else {
                    graph[incoming_nodes.at(i)].outgoing_edges[node1].at(0).multiplicity = in_multis[i];
                    graph[node1].incoming_edges[incoming_nodes.at(i)].at(0).multiplicity = in_multis[i];
                }
            }
            for (int i = 0; i < out_multis.size(); ++i) {
                if (out_multis[i] == 0) {
                    graph[node2].outgoing_edges.erase(outgoing_nodes[i]);
                    graph[outgoing_nodes[i]].incoming_edges.erase(node2);
                }
                else {
                    graph[node2].outgoing_edges[outgoing_nodes[i]].at(0).multiplicity = out_multis[i];
                    graph[outgoing_nodes[i]].incoming_edges[node2].at(0).multiplicity = out_multis[i];
                }
            }
            if (in_multi == 0) {
                assert(out_multi == 0);
                nodes_to_remove.push_back(node1);
                nodes_to_remove.push_back(node2);
            }
        }
    }
}

void Graph::remove_whirl(Path& unambiguous_path, std::vector<std::string>& nodes_to_remove) {
    std::cout << "General whirl (min multi " << unambiguous_path.min_multi << ", len " << unambiguous_path.length << "): " << unambiguous_path.nodes.at(0);

    if (unambiguous_path.safe_to_extract) {
        // change multiplicity of cyclic unambiguous path
        for (int i = 1; i < unambiguous_path.nodes.size(); ++i) {
            int bulge_leg = unambiguous_path.bulge_legs[i - 1];
            graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].at(bulge_leg).remove_multi_from_path(unambiguous_path);
            graph[unambiguous_path.nodes[i]].incoming_edges[unambiguous_path.nodes[i - 1]].at(bulge_leg).remove_multi_from_path(unambiguous_path);
            // if multiplicity is reduced to 0, remove edge
            if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].at(bulge_leg).multiplicity <= MIN_MULTI) {
                graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].erase(graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].begin() + bulge_leg);
                graph[unambiguous_path.nodes[i]].incoming_edges[unambiguous_path.nodes[i - 1]].erase(graph[unambiguous_path.nodes[i]].incoming_edges[unambiguous_path.nodes[i - 1]].begin() + bulge_leg);
            }
            if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges[unambiguous_path.nodes[i]].empty()) {
                graph[unambiguous_path.nodes[i - 1]].outgoing_edges.erase(unambiguous_path.nodes[i]);
                graph[unambiguous_path.nodes[i]].incoming_edges.erase(unambiguous_path.nodes[i - 1]);
                // add nodes without edges to remove list; we should never remove node.first because there is a self-loop
                if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges.size() == 0 && graph[unambiguous_path.nodes[i - 1]].incoming_edges.size() == 0)
                    nodes_to_remove.push_back(unambiguous_path.nodes[i - 1]);
            }
            // if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges.size() == 0 && graph[unambiguous_path.nodes[i - 1]].incoming_edges.size() != 0)
            //     std::cout << "(Tip)";
            // if (graph[unambiguous_path.nodes[i - 1]].outgoing_edges.size() != 0 && graph[unambiguous_path.nodes[i - 1]].incoming_edges.size() == 0)
            //     std::cout << "(Tip)";
            std::cout << "->" << unambiguous_path.nodes.at(i);
        }
        std::cout << std::endl;
        // add a self-loop to node.first
        Edge edge = Edge(unambiguous_path.sequence.at(k), unambiguous_path.length, unambiguous_path.sequence, unambiguous_path.min_multi);
        edge.path_edges_in_original_graph = unambiguous_path.path_edges_in_original_graph;
        edge.path_nodes_in_original_graph = unambiguous_path.path_nodes_in_original_graph;
        graph[unambiguous_path.nodes[0]].outgoing_edges[unambiguous_path.nodes[0]].push_back(edge);
        graph[unambiguous_path.nodes[0]].incoming_edges[unambiguous_path.nodes[0]].push_back(edge);
    }
    else {
        assert(unambiguous_path.nodes.size() == 3);
        std::cout << "->" << unambiguous_path.nodes[1] << "->" << unambiguous_path.nodes[2] << std::endl;

        graph[unambiguous_path.nodes[0]].outgoing_edges[unambiguous_path.nodes[1]].erase(graph[unambiguous_path.nodes[0]].outgoing_edges[unambiguous_path.nodes[1]].begin() + unambiguous_path.bulge_legs[0]);
        graph[unambiguous_path.nodes[1]].incoming_edges[unambiguous_path.nodes[0]].erase(graph[unambiguous_path.nodes[1]].incoming_edges[unambiguous_path.nodes[0]].begin() + unambiguous_path.bulge_legs[0]);
        if (graph[unambiguous_path.nodes[0]].outgoing_edges[unambiguous_path.nodes[1]].empty()) {
            graph[unambiguous_path.nodes[0]].outgoing_edges.erase(unambiguous_path.nodes[1]);
            graph[unambiguous_path.nodes[1]].incoming_edges.erase(unambiguous_path.nodes[0]);
        }

        int max_len = 0;
        double multi = 0;
        for (auto&& n : graph[unambiguous_path.nodes[1]].incoming_edges) {
            for (auto&& e : n.second) {
                if (n.first != unambiguous_path.nodes[0] && n.first != unambiguous_path.nodes[1] && e.length > max_len) {
                    multi = e.multiplicity;
                    max_len = e.length;
                }
            }
        }
        for (auto&& n : graph[unambiguous_path.nodes[0]].outgoing_edges) {
            for (auto&& e : n.second) {
                if (n.first != unambiguous_path.nodes[1] && n.first != unambiguous_path.nodes[0] && e.length > max_len) {
                    multi = e.multiplicity;
                    max_len = e.length;
                }
            }
        }
        graph[unambiguous_path.nodes[1]].outgoing_edges[unambiguous_path.nodes[2]].at(unambiguous_path.bulge_legs[1]).multiplicity = multi;

        Edge edge = Edge(unambiguous_path.sequence.at(k), unambiguous_path.length, unambiguous_path.sequence, unambiguous_path.min_multi);
        edge.path_edges_in_original_graph = unambiguous_path.path_edges_in_original_graph;
        edge.path_nodes_in_original_graph = unambiguous_path.path_nodes_in_original_graph;
        graph[unambiguous_path.nodes[0]].outgoing_edges[unambiguous_path.nodes[0]].push_back(edge);
        graph[unambiguous_path.nodes[0]].incoming_edges[unambiguous_path.nodes[0]].push_back(edge);

        Path unambiguous_path_r;
        assert(get_reverse_path(unambiguous_path, unambiguous_path_r));
        std::cout << "General whirl (min multi " << unambiguous_path_r.min_multi << ", len " << unambiguous_path_r.length << "): " << unambiguous_path_r.nodes.at(0) << "->" << unambiguous_path_r.nodes[1] << "->" << unambiguous_path_r.nodes[2] << std::endl;
        graph[unambiguous_path_r.nodes[1]].outgoing_edges[unambiguous_path_r.nodes[2]].erase(graph[unambiguous_path_r.nodes[1]].outgoing_edges[unambiguous_path_r.nodes[2]].begin() + unambiguous_path_r.bulge_legs[1]);
        graph[unambiguous_path_r.nodes[2]].incoming_edges[unambiguous_path_r.nodes[1]].erase(graph[unambiguous_path_r.nodes[2]].incoming_edges[unambiguous_path_r.nodes[1]].begin() + unambiguous_path_r.bulge_legs[1]);
        if (graph[unambiguous_path_r.nodes[1]].outgoing_edges[unambiguous_path_r.nodes[2]].empty()) {
            graph[unambiguous_path_r.nodes[1]].outgoing_edges.erase(unambiguous_path_r.nodes[2]);
            graph[unambiguous_path_r.nodes[2]].incoming_edges.erase(unambiguous_path_r.nodes[1]);
        }

        max_len = 0;
        double multi_r = 0;
        for (auto&& n : graph[unambiguous_path_r.nodes[1]].outgoing_edges) {
            for (auto&& e : n.second) {
                if (n.first != unambiguous_path_r.nodes[1] && n.first != unambiguous_path_r.nodes[0] && e.length > max_len) {
                    multi_r = e.multiplicity;
                    max_len = e.length;
                }
            }
        }
        for (auto&& n : graph[unambiguous_path_r.nodes[0]].incoming_edges) {
            for (auto&& e : n.second) {
                if (n.first != unambiguous_path_r.nodes[0] && n.first != unambiguous_path_r.nodes[1] && e.length > max_len) {
                    multi_r = e.multiplicity;
                    max_len = e.length;
                }
            }
        }
        if (multi_r != multi) {
            write_graph("debug");
            std::cout << multi_r << " " << multi << std::endl;
        }
        assert(multi_r == multi);
        graph[unambiguous_path_r.nodes[0]].outgoing_edges[unambiguous_path_r.nodes[1]].at(unambiguous_path_r.bulge_legs[0]).multiplicity = multi_r;

        Edge edge_r = Edge(unambiguous_path_r.sequence.at(k), unambiguous_path_r.length, unambiguous_path_r.sequence, unambiguous_path_r.min_multi);
        edge_r.path_edges_in_original_graph = unambiguous_path_r.path_edges_in_original_graph;
        edge_r.path_nodes_in_original_graph = unambiguous_path_r.path_nodes_in_original_graph;
        graph[unambiguous_path_r.nodes[2]].outgoing_edges[unambiguous_path_r.nodes[2]].push_back(edge_r);
        graph[unambiguous_path_r.nodes[2]].incoming_edges[unambiguous_path_r.nodes[2]].push_back(edge_r);
        assert(edge.multiplicity == edge_r.multiplicity && edge.sequence == reverse_complementary(edge_r.sequence));
    }
}

// remove general whirls
void Graph::general_whirl_removal(unsigned& removed_whirls, bool simple_whirl, bool force) {
    unsigned removed_bulges = 1;
    while (removed_bulges) {
        this->multi_bulge_removal(removed_bulges);
    }
    std::vector<std::string> nodes_to_remove;

    removed_whirls = 0;
    for (auto&& node : graph) {
        std::vector<std::string> outgoing_edges;
        for (auto&& successor : node.second.outgoing_edges) {
            // skip bulges
            if (successor.second.size() == 1)
                outgoing_edges.push_back(successor.first);
        }

        for (auto&& successor : outgoing_edges) {
            Path unambiguous_path;
            std::string successor_node = successor;
            if (successor_node == node.first || (force == false && successor_node == reverse_complementary_node(node.first)))
                continue;
            if (check_non_branching(successor_node))
                continue;
            if (!this->add_node_to_path(unambiguous_path, node.first))
                continue;

            std::string prev_node = node.first;
            // successor_node has unambiguous outgoing edge
            while (graph[successor_node].outgoing_edges.size() == 1) {
                // add successor_node to the unambiguous path
                if (!this->add_node_to_path(unambiguous_path, successor_node))
                    break;

                std::string node_after_successor;
                for (auto&& n : graph[successor_node].outgoing_edges) {
                    node_after_successor = n.first;
                }
                // skip bulges
                if (this->graph[successor_node].outgoing_edges[node_after_successor].size() > 1)
                    break;

                // no back edges
                bool flag = false;
                for (int i = 1; i < unambiguous_path.nodes.size(); ++i) {
                    if (node_after_successor == unambiguous_path.nodes.at(i) || (force == false && node_after_successor == reverse_complementary_node(unambiguous_path.nodes.at(i)))) {
                        flag = true;
                        break;
                    }
                }
                if (flag || (force == false && node_after_successor == reverse_complementary_node(unambiguous_path.nodes.at(0))))
                    break;

                if (check_non_branching(node_after_successor))
                    continue;

                // find a cyclic unambiguous path
                if (node_after_successor == node.first) {
                    if (!this->add_node_to_path(unambiguous_path, node_after_successor))
                        break;
                    // process the reverse complementary whirl
                    Path unambiguous_path_reverse;
                    assert(get_reverse_path(unambiguous_path, unambiguous_path_reverse));

                    // assert(unambiguous_path.color == unambiguous_path_reverse.color && unambiguous_path.multiplicity == unambiguous_path_reverse.multiplicity);
                    if (simple_whirl && unambiguous_path.nodes.size() > 3)
                        break;
                    if (unambiguous_path.nodes.size() > 3 && !unambiguous_path.safe_to_extract)
                        break;
                    if (force == false && !unambiguous_path.safe_to_extract)
                        break;
                    // std::cout << unambiguous_path.safe_to_extract << std::endl;
                    this->remove_whirl(unambiguous_path, nodes_to_remove);
                    if (unambiguous_path.safe_to_extract)
                        this->remove_whirl(unambiguous_path_reverse, nodes_to_remove);
                    // this->write_graph("jumbodbg.chromosomes.k8001.consensus_new/whirl." + unambiguous_path.nodes.at(0));
                    removed_whirls += 2;
                    break;
                }

                // update previous node and successor node
                prev_node = successor_node;
                successor_node = node_after_successor;
            }
        }
    }

    for (auto&& n : nodes_to_remove) {
        this->graph.erase(n);
    }
    this->merge_non_branching_paths();
}

// deleting elements from the start and end of second sequence is "free"
int Graph::matches_by_edlib(std::string sequence1, std::string sequence2) {
    EdlibAlignResult result = edlibAlign(sequence1.c_str(), sequence1.size(), sequence2.c_str(), sequence2.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        std::string cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        edlibFreeAlignResult(result);
        int lcs_len = count_matches(cigar);
        return lcs_len;
    }
    else {
        std::cout << "edlib failed" << std::endl;
        return -1;
    }
}

std::string Graph::collapse_complex_bulge_two_multi_edge_paths(Path p1, Path p2, bool p1_2_in_2_out, bool p2_2_in_2_out, double max_identity, std::vector<std::string>& nodes_to_remove) {
    std::cout << "Complex bulge: path1 (length " << p1.length << ", min multi " << p1.min_multi << "): " << p1.nodes[0];
    for (int i = 1; i < p1.nodes.size();++i)
        std::cout << "->" << p1.nodes[i];
    std::cout << "; path2 (length " << p2.length << ", min multi " << p2.min_multi << "): " << p2.nodes[0];
    for (int i = 1; i < p2.nodes.size();++i)
        std::cout << "->" << p2.nodes[i];
    std::cout << "; " << "identity " << max_identity << std::endl;

    int extract_index = 0;
    if (p1_2_in_2_out == false && p2_2_in_2_out == false) {
        if (p1.safe_to_extract && p2.safe_to_extract)
            extract_index = 3;
        else if (p1.safe_to_extract)
            extract_index = 1;
        else
            extract_index = 2;
    }
    else if (p2_2_in_2_out == false) {
        extract_index = 2;
        if (!p2.safe_to_extract)
            return "";
    }
    else if (p1_2_in_2_out == false) {
        extract_index = 1;
        if (!p1.safe_to_extract)
            return "";
    }
    else {
        return "";
    }

    if (extract_index == 3) {
        // decrease multiplicity for p1
        if (!p1.safe_to_extract) {
            std::cout << "Potential tip" << std::endl;
        }
        for (int i = 1; i < p1.nodes.size();++i) {
            int bulge_leg = p1.bulge_legs.at(i - 1);
            graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p1);
            graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p1);
            if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).multiplicity <= MIN_MULTI) {
                graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].erase(graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].begin() + bulge_leg);
                graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].erase(graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].begin() + bulge_leg);
            }
            if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].empty()) {
                graph[p1.nodes.at(i - 1)].outgoing_edges.erase(p1.nodes.at(i));
                graph[p1.nodes.at(i)].incoming_edges.erase(p1.nodes.at(i - 1));
                if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() == 0)
                    nodes_to_remove.emplace_back(p1.nodes.at(i - 1));
            }
            if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() != 0)
                std::cout << "Tip induced: " << p1.nodes.at(i - 1) << std::endl;
            if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() != 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() == 0)
                std::cout << "Tip induced: " << p1.nodes.at(i - 1) << std::endl;
        }

        // decrease multiplicity for p2
        if (!p2.safe_to_extract) {
            std::cout << "Potential tip" << std::endl;
        }
        for (int i = 1; i < p2.nodes.size();++i) {
            int bulge_leg = p2.bulge_legs.at(i - 1);
            graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p2);
            graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p2);
            if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).multiplicity <= MIN_MULTI) {
                graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].erase(graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].begin() + bulge_leg);
                graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].erase(graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].begin() + bulge_leg);
            }
            if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].empty()) {
                graph[p2.nodes.at(i - 1)].outgoing_edges.erase(p2.nodes.at(i));
                graph[p2.nodes.at(i)].incoming_edges.erase(p2.nodes.at(i - 1));
                if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() == 0)
                    nodes_to_remove.emplace_back(p2.nodes.at(i - 1));
            }
            if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() != 0)
                std::cout << "Tip induced: " << p2.nodes.at(i - 1) << std::endl;
            if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() != 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() == 0)
                std::cout << "Tip induced: " << p2.nodes.at(i - 1) << std::endl;
        }

        // add an edge from node.first to node.last with multiplicity p1 + p2
        std::string seq;
        if (p1.nodes.size() == 2 || (p2.nodes.size() > 2 && p1.length > p2.length)) {
            Edge edge(p1.sequence.at(k), p1.length, p1.sequence, p1.min_multi);
            edge.path_nodes_in_original_graph = p1.path_nodes_in_original_graph;
            edge.path_edges_in_original_graph = p1.path_edges_in_original_graph;
            // p2.min_multi = 1.0 * p2.min_multi * p2.length / p1.length;
            edge.add_multi_from_edge_or_path(p2);
            graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].emplace_back(edge);
            graph[p1.nodes.at(p1.nodes.size() - 1)].incoming_edges[p1.nodes.at(0)].emplace_back(edge);
            seq = edge.sequence;
            // std::cout << "Resulting multi " << edge.multiplicity << std::endl;
        }
        else {
            Edge edge(p2.sequence.at(k), p2.length, p2.sequence, p2.min_multi);
            edge.path_edges_in_original_graph = p2.path_edges_in_original_graph;
            edge.path_nodes_in_original_graph = p2.path_nodes_in_original_graph;
            // p1.min_multi = 1.0 * p1.min_multi * p1.length / p2.length;
            edge.add_multi_from_edge_or_path(p1);
            graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].emplace_back(edge);
            graph[p1.nodes.at(p1.nodes.size() - 1)].incoming_edges[p1.nodes.at(0)].emplace_back(edge);
            seq = edge.sequence;
            // std::cout << "Resulting multi " << edge.multiplicity << std::endl;
        }

        // if it forms a simple bulge, remove it immediately
        // if (graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].size() >= 2 && p1.nodes.at(0) != reverse_complementary_node(p1.nodes.at(p1.nodes.size() - 1))) {
        //     unsigned removed_bulges = 0;
        //     seq = this->collapse_bulge(p1.nodes.at(0), p1.nodes.at(p1.nodes.size() - 1), removed_bulges);
        // }

        return seq;
    }
    else if (extract_index == 1) {
        // decrease multiplicity for p1
        if (!p1.safe_to_extract) {
            std::cout << "Potential tip" << std::endl;
        }
        for (int i = 1; i < p1.nodes.size();++i) {
            int bulge_leg = p1.bulge_legs.at(i - 1);
            graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p1);
            graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p1);
            if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).multiplicity <= MIN_MULTI) {
                graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].erase(graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].begin() + bulge_leg);
                graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].erase(graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].begin() + bulge_leg);
            }
            if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].empty()) {
                graph[p1.nodes.at(i - 1)].outgoing_edges.erase(p1.nodes.at(i));
                graph[p1.nodes.at(i)].incoming_edges.erase(p1.nodes.at(i - 1));
                if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() == 0)
                    nodes_to_remove.emplace_back(p1.nodes.at(i - 1));
            }
            if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() != 0)
                std::cout << "Tip induced: " << p1.nodes.at(i - 1) << std::endl;
            if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() != 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() == 0)
                std::cout << "Tip induced: " << p1.nodes.at(i - 1) << std::endl;
        }

        // add multiplicity for p2
        // p1.min_multi = 1.0 * p1.min_multi * p1.length / p2.length;
        for (int i = 1; i < p2.nodes.size();++i) {
            int bulge_leg = p2.bulge_legs.at(i - 1);
            graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).add_multi_from_edge_or_path(p1);
            graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].at(bulge_leg).add_multi_from_edge_or_path(p1);
        }
        std::string seq = p2.sequence;
        // std::cout << "Multi change " << p1.min_multi << std::endl;
        return seq;
    }
    else if (extract_index == 2) {
        // decrease multiplicity for p2
        if (!p2.safe_to_extract) {
            std::cout << "Potential tip" << std::endl;
        }
        for (int i = 1; i < p2.nodes.size();++i) {
            int bulge_leg = p2.bulge_legs.at(i - 1);
            graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p2);
            graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p2);
            if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).multiplicity <= MIN_MULTI) {
                graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].erase(graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].begin() + bulge_leg);
                graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].erase(graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].begin() + bulge_leg);
            }
            if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].empty()) {
                graph[p2.nodes.at(i - 1)].outgoing_edges.erase(p2.nodes.at(i));
                graph[p2.nodes.at(i)].incoming_edges.erase(p2.nodes.at(i - 1));
                if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() == 0)
                    nodes_to_remove.emplace_back(p2.nodes.at(i - 1));
            }
            if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() != 0)
                std::cout << "Tip induced: " << p2.nodes.at(i - 1) << std::endl;
            if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() != 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() == 0)
                std::cout << "Tip induced: " << p2.nodes.at(i - 1) << std::endl;
        }
        // add multiplicity for p1
        // p2.min_multi = 1.0 * p2.min_multi * p2.length / p1.length;
        for (int i = 1; i < p1.nodes.size();++i) {
            int bulge_leg = p1.bulge_legs.at(i - 1);
            graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).add_multi_from_edge_or_path(p2);
            graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].at(bulge_leg).add_multi_from_edge_or_path(p2);
        }
        std::string seq = p1.sequence;
        // std::cout << "Multi change " << p2.min_multi << std::endl;
        return seq;
    }
    else
        return "";
}

bool Graph::process_palindromic_bulges(Path& p1, Path& p2, std::vector<std::string>& nodes_to_remove, bool verbose, std::string output) {
    std::cout << "Palindromic bulge: path1 (length " << p1.length << ", min multi " << p1.multiplicity << "): " << p1.nodes[0];
    for (int i = 1; i < p1.nodes.size();++i)
        std::cout << "->" << p1.nodes[i];
    std::cout << "; path2 (length " << p2.length << ", min multi " << p2.multiplicity << ", color " << "): " << p2.nodes[0];
    for (int i = 1; i < p2.nodes.size();++i)
        std::cout << "->" << p2.nodes[i];
    std::cout << std::endl;

    // decrease multiplicity for p1
    for (int i = 1; i < p1.nodes.size();++i) {
        int bulge_leg = p1.bulge_legs.at(i - 1);
        if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].size() == 1)
            bulge_leg = 0;
        graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p1);
        graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p1);
        if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].at(bulge_leg).multiplicity == 0) {
            graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].erase(graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].begin() + bulge_leg);
            graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].erase(graph[p1.nodes.at(i)].incoming_edges[p1.nodes.at(i - 1)].begin() + bulge_leg);
        }
        if (graph[p1.nodes.at(i - 1)].outgoing_edges[p1.nodes.at(i)].empty()) {
            graph[p1.nodes.at(i - 1)].outgoing_edges.erase(p1.nodes.at(i));
            graph[p1.nodes.at(i)].incoming_edges.erase(p1.nodes.at(i - 1));
            if (p1.nodes.at(i - 1) != p1.nodes.at(0) && graph[p1.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p1.nodes.at(i - 1)].outgoing_edges.size() == 0)
                nodes_to_remove.emplace_back(p1.nodes.at(i - 1));
        }
    }

    // decrease multiplicity for p2
    for (int i = 1; i < p2.nodes.size();++i) {
        int bulge_leg = p2.bulge_legs.at(i - 1);
        if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].size() == 1)
            bulge_leg = 0;
        graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).remove_multi_from_path(p2);
        graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].at(bulge_leg).remove_multi_from_path(p2);
        if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].at(bulge_leg).multiplicity == 0) {
            graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].erase(graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].begin() + bulge_leg);
            graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].erase(graph[p2.nodes.at(i)].incoming_edges[p2.nodes.at(i - 1)].begin() + bulge_leg);
        }
        if (graph[p2.nodes.at(i - 1)].outgoing_edges[p2.nodes.at(i)].empty()) {
            graph[p2.nodes.at(i - 1)].outgoing_edges.erase(p2.nodes.at(i));
            graph[p2.nodes.at(i)].incoming_edges.erase(p2.nodes.at(i - 1));
            if (p2.nodes.at(i - 1) != p2.nodes.at(0) && graph[p2.nodes.at(i - 1)].incoming_edges.size() == 0 && graph[p2.nodes.at(i - 1)].outgoing_edges.size() == 0)
                nodes_to_remove.emplace_back(p2.nodes.at(i - 1));
        }
    }

    Edge edge1(p1.sequence.at(k), p1.length, p1.sequence, p1.multiplicity);
    graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].emplace_back(edge1);
    graph[p1.nodes.at(p1.nodes.size() - 1)].incoming_edges[p1.nodes.at(0)].emplace_back(edge1);

    Edge edge2(p2.sequence.at(k), p2.length, p2.sequence, p2.multiplicity);
    graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].emplace_back(edge2);
    graph[p1.nodes.at(p1.nodes.size() - 1)].incoming_edges[p1.nodes.at(0)].emplace_back(edge2);

    if (verbose)
        write_graph(output + "/graph.palindromic_bulge." + p1.nodes.at(0) + "_" + p1.nodes.at(p1.nodes.size() - 1));

    return true;
}

bool Graph::get_reverse_path(Path& path, Path& path_reverse) {
    bool flag = true;
    this->add_node_to_path(path_reverse, reverse_complementary_node(path.nodes[path.nodes.size() - 1]));

    for (int i = int(path.nodes.size()) - 2; i >= 0; --i) {
        // to process palindromic simple bulge
        if (graph[path_reverse.nodes[path_reverse.nodes.size() - 1]].outgoing_edges[reverse_complementary_node(path.nodes[i])].size() > 1 && path_reverse.nodes[path_reverse.nodes.size() - 1] == path.nodes[i]) {
            std::vector<Edge>& edges = graph[path_reverse.nodes[path_reverse.nodes.size() - 1]].outgoing_edges[reverse_complementary_node(path.nodes[i])];

            std::string seq_forward = path.index2palindromic_seq.at(i);
            std::string seq_reverse = reverse_complementary(seq_forward);

            int bulge_leg = 0;
            bool find_reverse = false;
            for (int e = 0; e < edges.size(); ++e) {
                if (edges.at(e).sequence == seq_reverse) {
                    bulge_leg = e;
                    find_reverse = true;
                    break;
                }
            }
            if (!find_reverse) {
                flag = false;
                std::cout << "Error: reverse edge does not exist" << std::endl;
                break;
            }

            if (!this->add_node_to_path(path_reverse, reverse_complementary_node(path.nodes[i]), bulge_leg)) {
                flag = false;
                break;
            }

            path_reverse.index2palindromic_seq[path_reverse.nodes.size() - 2] = seq_reverse;
        }
        else if (!this->add_node_to_path(path_reverse, reverse_complementary_node(path.nodes[i]), path.bulge_legs[i])) {
            flag = false;
            break;
        }
    }
    if (!flag)
        return false;

    if (path.min_multi - path_reverse.min_multi > MIN_MULTI) {
        write_graph("debug");
        std::cout << "p_f" << std::endl;
        for (auto&& n : path.nodes) {
            std::cout << n << std::endl;
        }
        std::cout << "p_r" << std::endl;
        for (auto&& n : path_reverse.nodes) {
            std::cout << n << std::endl;
        }
        std::cout << path.min_multi << ' ' << path_reverse.min_multi << std::endl;
    }
    assert(path.min_multi - path_reverse.min_multi <= MIN_MULTI);
    path_reverse.min_multi = path.min_multi;
    // assert(path.safe_to_extract == path_reverse.safe_to_extract);
    // assert(path.color == path_reverse.color);
    return true;
}

// collapse bulges with two paths of multiple edges (at most x)
void Graph::resolving_bulge_with_two_multi_edge_paths(unsigned& removed_paths, int x, double identity, bool use_length, int security_level, bool allow_reverse_comp, bool verbose) {
    // this->write_graph("debug");
    unsigned removed_whirls = 1;
    while (removed_whirls != 0) {
        general_whirl_removal(removed_whirls);
    }

    if (allow_reverse_comp)
        multi_bulge_removal(removed_whirls, false);

    removed_paths = 0;
    std::vector<std::string> nodes_to_remove;

    // store all multi-edge bulges
    std::vector<Bulge> bulges;
    for (auto&& node : graph) {
        if (node.second.outgoing_edges.size() <= 1)
            continue;
        // store all paths starts from node and and at key
        std::unordered_map<std::string, std::vector<Path>> all_paths_ending_at_key;
        std::queue<Path> paths_bfs;
        Path path;
        this->add_node_to_path(path, node.first);
        paths_bfs.push(path);

        // maximum length is limited at x
        while (paths_bfs.front().nodes.size() <= x) {
            if (paths_bfs.empty()) break;
            Path prev_path = paths_bfs.front();
            paths_bfs.pop();

            std::string prev_node = prev_path.nodes.at(prev_path.nodes.size() - 1);
            for (auto&& successor : this->graph[prev_node].outgoing_edges) {
                // ignore back edges, this will ignore self-loops as well
                bool flag = false;
                for (auto&& i : prev_path.nodes) {
                    if (successor.first == i)
                        flag = true;
                    if (!allow_reverse_comp && successor.first == reverse_complementary_node(i))
                        flag = true;
                }
                if (flag) continue;

                for (int i = 0; i < successor.second.size(); ++i) {
                    Path current_path = prev_path;
                    if (!this->add_node_to_path(current_path, successor.first, i))
                        continue;

                    paths_bfs.push(current_path);
                    all_paths_ending_at_key[successor.first].emplace_back(current_path);
                }
            }
        }

        // add bulges to container
        for (auto&& ending_node : all_paths_ending_at_key) {
            // check whether there are multiple paths
            if (ending_node.second.size() < 2)
                continue;

            // store all the remaining paths
            std::vector<Path> paths = ending_node.second;
            // check whether there are multiple remaining paths
            if (paths.size() < 2)
                continue;

            // select two paths with the highest identity
            double max_identity = 0;
            int max_lcs_len = 0;
            Path p1, p2;
            for (int i = 0;i < paths.size();++i) {
                std::set<std::string> s1;
                for (int k = 1; k < paths[i].nodes.size() - 1; ++k)
                    s1.insert(paths[i].nodes[k]);
                for (int j = i + 1; j < paths.size(); ++j) {
                    std::string sequence1 = paths[i].sequence, sequence2 = paths[j].sequence;
                    // the two paths should have no shared nodes
                    if (allow_reverse_comp) {
                        bool flag = false;
                        std::set<std::string> s2;
                        for (int k = 1; k < paths[j].nodes.size() - 1; ++k)
                            s2.insert(paths[j].nodes[k]);
                        for (auto&& n : s1) {
                            if (s2.find(n) != s2.end() || s2.find(reverse_complementary_node(n)) != s2.end())
                                flag = true;
                        }
                        std::set<std::pair<std::string, std::string>> p1_pairs;

                        for (size_t x = 0; x < paths[i].nodes.size() - 1; ++x) {
                            p1_pairs.insert({ paths[i].nodes[x], paths[i].nodes[x + 1] });
                            p1_pairs.insert({ reverse_complementary_node(paths[i].nodes[x + 1]), reverse_complementary_node(paths[i].nodes[x]) });
                        }

                        for (size_t x = 0; x < paths[j].nodes.size() - 1; ++x) {
                            if (p1_pairs.find({ paths[j].nodes[x], paths[j].nodes[x + 1] }) != p1_pairs.end()) {
                                flag = true; // Pair found
                            }
                        }
                        if (flag) continue;
                    }
                    else {
                        bool flag = false;
                        for (auto&& n : paths[j].nodes) {
                            if (s1.find(n) != s1.end())
                                flag = true;
                        }
                        if (flag) continue;
                        std::set<std::string> s2;
                        for (int k = 1; k < paths[j].nodes.size() - 1; ++k)
                            s2.insert(paths[j].nodes[k]);
                        for (auto&& n : paths[i].nodes) {
                            if (s2.find(n) != s2.end())
                                flag = true;
                        }
                        if (flag) continue;
                    }

                    // check reverse paths
                    Path p1_reverse, p2_reverse;
                    assert(get_reverse_path(paths[i], p1_reverse));
                    assert(get_reverse_path(paths[j], p2_reverse));

                    // std::cout << "Run edlib for path1 (length " << paths[i].length << "): " << paths[i].nodes[0] << std::flush;
                    // for (int k = 1; k < paths[i].nodes.size();++k)
                    //     std::cout << "->" << paths[i].nodes[k];
                    // std::cout << "; path2 (length " << paths[j].length << "): " << paths[j].nodes[0];
                    // for (int k = 1; k < paths[j].nodes.size();++k)
                    //     std::cout << "->" << paths[j].nodes[k];

                    // calculate identity
                    double alignment_identity = 0;
                    int lcs_len = 0;
                    if (sequence1.size() > sequence2.size() && 1.0 * (sequence1.size() - sequence2.size()) / sequence1.size() > 1 - identity) {
                        // std::cout << std::endl;
                        continue;
                    }
                    if (sequence2.size() > sequence1.size() && 1.0 * (sequence2.size() - sequence1.size()) / sequence2.size() > 1 - identity) {
                        // std::cout << std::endl;
                        continue;
                    }
                    if (use_length)
                        lcs_len = std::min(sequence1.size(), sequence2.size());
                    else
                        lcs_len = matches_by_edlib(sequence1, sequence2);
                    alignment_identity = 1.0 * lcs_len / std::max(sequence1.size(), sequence2.size());

                    // record max identity paths to p1 and p2
                    // std::cout << "; identity " << alignment_identity << std::endl;
                    if (alignment_identity > max_identity) {
                        p1 = paths[i];
                        p2 = paths[j];
                        max_identity = alignment_identity;
                        max_lcs_len = lcs_len;
                    }
                }
            }
            // the two paths cannot pass similarity check
            if (p1.nodes.empty() || p2.nodes.empty())
                continue;
            if (p1.sequence != reverse_complementary(p2.sequence) && max_identity < identity)
                continue;

            // check confict
            Bulge b(p1, p2, max_identity);
            this->find_2_in_2_out(b);
            for (int i = 0; i < bulges.size(); ++i) {
                b.check_conflict(bulges[i]);
            }
            // allow bulges without 2-in-2-out edges
            if (security_level == 3) {
                if (b.path_resolutions1.size() + b.path_resolutions2.size() == 0)
                    bulges.emplace_back(b);
            }
            else if (security_level == 2) {
                if (b.path_resolutions1.size() == 0 || b.path_resolutions2.size() == 0)
                    bulges.emplace_back(b);
            }
            else {
                bulges.emplace_back(b);
            }
        }
    }

    std::sort(bulges.begin(), bulges.end(),
        [](const Bulge& lhs, const Bulge& rhs) {
            if (lhs.path_resolutions1.size() + lhs.path_resolutions2.size() == 0 && rhs.path_resolutions1.size() + rhs.path_resolutions2.size() > 0)
                return true;
            if (lhs.path_resolutions1.size() + lhs.path_resolutions2.size() > 0 && rhs.path_resolutions1.size() + rhs.path_resolutions2.size() == 0)
                return false;
            if (lhs.is_confict == false && rhs.is_confict == true)
                return true;
            if (lhs.is_confict == true && rhs.is_confict == false)
                return false;
            if (lhs.num_edges < rhs.num_edges)
                return true;
            if (lhs.num_edges > rhs.num_edges)
                return false;
            return lhs.max_identity > rhs.max_identity;
        }
    );

    // process the paths
    for (auto&& bulge : bulges) {
        Path p1, p2;
        bool flag = true;
        // check forward paths
        for (int i = 0; i < bulge.leg1.nodes.size(); ++i) {
            if (i == 0) {
                if (!this->add_node_to_path(p1, bulge.leg1.nodes[i])) {
                    flag = false;
                    break;
                }
            }
            else if (!this->add_node_to_path(p1, bulge.leg1.nodes[i], bulge.leg1.bulge_legs[i - 1])) {
                flag = false;
                break;
            }
        }
        for (int i = 0; i < bulge.leg2.nodes.size(); ++i) {
            if (i == 0) {
                if (!this->add_node_to_path(p2, bulge.leg2.nodes[i])) {
                    flag = false;
                    break;
                }
            }
            else if (!this->add_node_to_path(p2, bulge.leg2.nodes[i], bulge.leg2.bulge_legs[i - 1])) {
                flag = false;
                break;
            }
        }

        if (!flag)
            continue;

        // check reverse paths
        Path p1_reverse, p2_reverse;
        assert(get_reverse_path(p1, p1_reverse));
        assert(get_reverse_path(p2, p2_reverse));

        // modify multiplicity in case of reverse complementary -1990600->-5726575 in 5726575->1990600->741685->-1990600->-5726575->6492462
        std::set<std::pair<std::string, std::string>> p_pairs;
        for (size_t x = 0; x < p1.nodes.size() - 1; ++x) {
            p_pairs.insert({ reverse_complementary_node(p1.nodes[x + 1]), reverse_complementary_node(p1.nodes[x]) });
        }
        for (size_t x = 0; x < p1_reverse.nodes.size() - 1; ++x) {
            if (p_pairs.find({ p1_reverse.nodes[x], p1_reverse.nodes[x + 1] }) != p_pairs.end()) {
                p1.min_multi = std::min(p1.min_multi, graph[p1_reverse.nodes[x]].outgoing_edges[p1_reverse.nodes[x + 1]].at(p1_reverse.bulge_legs[x]).multiplicity / 2);
                p1_reverse.min_multi = std::min(p1_reverse.min_multi, graph[p1_reverse.nodes[x]].outgoing_edges[p1_reverse.nodes[x + 1]].at(p1_reverse.bulge_legs[x]).multiplicity / 2);
            }
        }
        assert(std::abs(p1.min_multi - p1_reverse.min_multi) < MIN_MULTI);
        p_pairs.clear();
        for (size_t x = 0; x < p2.nodes.size() - 1; ++x) {
            p_pairs.insert({ reverse_complementary_node(p2.nodes[x + 1]), reverse_complementary_node(p2.nodes[x]) });
        }
        for (size_t x = 0; x < p2_reverse.nodes.size() - 1; ++x) {
            if (p_pairs.find({ p2_reverse.nodes[x], p2_reverse.nodes[x + 1] }) != p_pairs.end()) {
                p2.min_multi = std::min(p2.min_multi, graph[p2_reverse.nodes[x]].outgoing_edges[p2_reverse.nodes[x + 1]].at(p2_reverse.bulge_legs[x]).multiplicity / 2);
                p2_reverse.min_multi = std::min(p2_reverse.min_multi, graph[p2_reverse.nodes[x]].outgoing_edges[p2_reverse.nodes[x + 1]].at(p2_reverse.bulge_legs[x]).multiplicity / 2);
            }
        }
        assert(std::abs(p2.min_multi - p2_reverse.min_multi) < MIN_MULTI);

        if (!p1.safe_to_extract && !p2.safe_to_extract)
            continue;


        if (verbose)
            write_graph("debug_" + p1.nodes.at(0) + "_" + p1.nodes.at(p1.nodes.size() - 1) + "_before");
        bool p1_2_in_2_out = check_2_in_2_out(p1);
        bool p2_2_in_2_out = check_2_in_2_out(p2);
        std::string seq1 = this->collapse_complex_bulge_two_multi_edge_paths(p1, p2, p1_2_in_2_out, p2_2_in_2_out, bulge.max_identity, nodes_to_remove);

        // Path p1_reverse, p2_reverse;
        // assert(get_reverse_path(p1, p1_reverse));
        // assert(get_reverse_path(p2, p2_reverse));
        if (p2_reverse.nodes != p1.nodes) {
            std::string seq2 = this->collapse_complex_bulge_two_multi_edge_paths(p1_reverse, p2_reverse, p1_2_in_2_out, p2_2_in_2_out, bulge.max_identity, nodes_to_remove);

            // if (seq1 != reverse_complementary(seq2))
            //     std::cout << "The two bulges did not result in reverse complementary sequence" << std::endl;
            if (seq1 == "" || seq2 == "")
                continue;
            removed_paths += 1;
        }
        if (seq1 == "")
            continue;
        if (graph[p1.nodes.at(0)].outgoing_edges.find(p1.nodes.at(p1.nodes.size() - 1)) != graph[p1.nodes.at(0)].outgoing_edges.end()) {
            if (graph[p1.nodes.at(0)].outgoing_edges[p1.nodes.at(p1.nodes.size() - 1)].size() >= 2) {
                unsigned removed_bulges = 0;
                this->collapse_bulge(p1.nodes.at(0), p1.nodes.at(p1.nodes.size() - 1), removed_bulges);
            }
            if (graph[p1_reverse.nodes.at(0)].outgoing_edges[p1_reverse.nodes.at(p1_reverse.nodes.size() - 1)].size() >= 2) {
                unsigned removed_bulges = 0;
                this->collapse_bulge(p1_reverse.nodes.at(0), p1_reverse.nodes.at(p1_reverse.nodes.size() - 1), removed_bulges);
            }
        }
        removed_paths += 1;

        for (auto&& node : nodes_to_remove) {
            this->graph.erase(node);
        }
        nodes_to_remove.clear();
        merge_non_branching_paths();
        if (verbose)
            write_graph("debug_" + p1.nodes.at(0) + "_" + p1.nodes.at(p1.nodes.size() - 1) + "_after");
    }
}

void Graph::decouple_strands(std::string node1, std::string node2, int& decoupled_stands, std::vector<std::string>& nodes_to_remove, bool strict) {
    std::vector<std::string> incoming_nodes, outgoing_nodes_tmp, outgoing_nodes;

    for (auto&& e : this->graph[node1].incoming_edges) {
        incoming_nodes.push_back(e.first);
    }
    for (auto&& e : this->graph[node2].outgoing_edges) {
        outgoing_nodes_tmp.push_back(e.first);
    }
    if (incoming_nodes.size() != 2 || outgoing_nodes_tmp.size() != 2)
        return;

    if (incoming_nodes[0] == reverse_complementary_node(outgoing_nodes_tmp[0]) || incoming_nodes[1] == reverse_complementary_node(outgoing_nodes_tmp[1])) {
        outgoing_nodes = outgoing_nodes_tmp;
    }
    else if (incoming_nodes[0] == reverse_complementary_node(outgoing_nodes_tmp[1]) || incoming_nodes[1] == reverse_complementary_node(outgoing_nodes_tmp[0])) {
        outgoing_nodes.push_back(outgoing_nodes_tmp[1]);
        outgoing_nodes.push_back(outgoing_nodes_tmp[0]);
    }
    else
        return;

    if (node1 == reverse_complementary_node(incoming_nodes[0])
        || node1 == reverse_complementary_node(incoming_nodes[1])
        || node2 == reverse_complementary_node(outgoing_nodes[0])
        || node2 == reverse_complementary_node(outgoing_nodes[1])
        )
        return;
    if (strict) {
        if (incoming_nodes[0] != reverse_complementary_node(outgoing_nodes[0]))
            return;
        if (incoming_nodes[1] != reverse_complementary_node(outgoing_nodes[1]))
            return;
    }
    std::cout << "Decouple strands: " << "path1: " << incoming_nodes[0] << "->" << node1 << "->...->" << node2 << "->" << outgoing_nodes[1] << "; path2: " << incoming_nodes[1] << "->" << node1 << "->...->" << node2 << "->" << outgoing_nodes[0] << std::endl;

    // process for forward strand
    // remove a paired path from the graph
    Path path1;
    this->add_node_to_path(path1, incoming_nodes.at(0));
    std::string node = node1;
    while (node != node2) {
        this->add_node_to_path(path1, node);
        for (auto&& n : graph[node].outgoing_edges) {
            node = n.first;
        }
    }
    this->add_node_to_path(path1, node2);
    this->add_node_to_path(path1, outgoing_nodes.at(1));
    path1.multiplicity = 1.0 * (graph[incoming_nodes.at(0)].outgoing_edges[node1].at(0).multiplicity * graph[incoming_nodes.at(0)].outgoing_edges[node1].at(0).length
        + graph[node2].outgoing_edges[outgoing_nodes.at(1)].at(0).multiplicity * graph[node2].outgoing_edges[outgoing_nodes.at(1)].at(0).length)
        / (graph[incoming_nodes.at(0)].outgoing_edges[node1].at(0).length + graph[node2].outgoing_edges[outgoing_nodes.at(1)].at(0).length);

    Path path2;
    this->add_node_to_path(path2, incoming_nodes.at(1));
    node = node1;
    while (node != node2) {
        this->add_node_to_path(path2, node);
        for (auto&& n : graph[node].outgoing_edges) {
            node = n.first;
        }
    }
    this->add_node_to_path(path2, node2);
    this->add_node_to_path(path2, outgoing_nodes.at(0));
    path2.multiplicity = 1.0 * (graph[incoming_nodes.at(1)].outgoing_edges[node1].at(0).multiplicity * graph[incoming_nodes.at(1)].outgoing_edges[node1].at(0).length
        + graph[node2].outgoing_edges[outgoing_nodes.at(0)].at(0).multiplicity * graph[node2].outgoing_edges[outgoing_nodes.at(0)].at(0).length)
        / (graph[incoming_nodes.at(1)].outgoing_edges[node1].at(0).length + graph[node2].outgoing_edges[outgoing_nodes.at(0)].at(0).length);

    Path path_1_r, path_2_r;
    assert(get_reverse_path(path1, path_1_r));
    assert(get_reverse_path(path2, path_2_r));
    path_1_r.multiplicity = path1.multiplicity;
    path_2_r.multiplicity = path2.multiplicity;

    Edge edge1(path1.sequence.at(k), path1.length, path1.sequence, path1.multiplicity);
    edge1.path_edges_in_original_graph = path1.path_edges_in_original_graph;
    edge1.path_nodes_in_original_graph = path1.path_nodes_in_original_graph;
    this->graph[incoming_nodes.at(0)].outgoing_edges[outgoing_nodes.at(1)].push_back(edge1);
    this->graph[outgoing_nodes.at(1)].incoming_edges[incoming_nodes.at(0)].push_back(edge1);

    Edge edge2(path2.sequence.at(k), path2.length, path2.sequence, path2.multiplicity);
    edge2.path_edges_in_original_graph = path2.path_edges_in_original_graph;
    edge2.path_nodes_in_original_graph = path2.path_nodes_in_original_graph;
    this->graph[incoming_nodes.at(1)].outgoing_edges[outgoing_nodes.at(0)].push_back(edge2);
    this->graph[outgoing_nodes.at(0)].incoming_edges[incoming_nodes.at(1)].push_back(edge2);

    for (int i = 0; i < path1.nodes.size() - 1; ++i) {
        this->graph[path1.nodes[i]].outgoing_edges.erase(path1.nodes[i + 1]);
        this->graph[path1.nodes[i + 1]].incoming_edges.erase(path1.nodes[i]);
    }

    for (int i = 0; i < path2.nodes.size() - 1; ++i) {
        this->graph[path2.nodes[i]].outgoing_edges.erase(path2.nodes[i + 1]);
        this->graph[path2.nodes[i + 1]].incoming_edges.erase(path2.nodes[i]);
    }

    for (int i = 0; i < path1.nodes.size(); ++i) {
        if (graph[path1.nodes[i]].incoming_edges.empty() && graph[path1.nodes[i]].outgoing_edges.empty() && std::find(nodes_to_remove.begin(), nodes_to_remove.end(), path1.nodes[i]) == nodes_to_remove.end())
            nodes_to_remove.push_back(path1.nodes[i]);
    }
    for (int i = 0; i < path2.nodes.size(); ++i) {
        if (graph[path2.nodes[i]].incoming_edges.empty() && graph[path2.nodes[i]].outgoing_edges.empty() && std::find(nodes_to_remove.begin(), nodes_to_remove.end(), path2.nodes[i]) == nodes_to_remove.end())
            nodes_to_remove.push_back(path2.nodes[i]);
    }


    // write_graph("debugging_" + path1.nodes.at(0) + "_" + path1.nodes.at(path1.nodes.size() - 1) + "_middle");

    if (path_1_r.nodes != path2.nodes) {
        std::vector<std::string> in_rev, out_rev;
        in_rev.push_back(reverse_complementary_node(outgoing_nodes.at(1)));
        in_rev.push_back(reverse_complementary_node(outgoing_nodes.at(0)));

        out_rev.push_back(reverse_complementary_node(incoming_nodes.at(0)));
        out_rev.push_back(reverse_complementary_node(incoming_nodes.at(1)));
        std::cout << "Decouple strands: " << "path1: " << in_rev[0] << "->" << reverse_complementary_node(node2) << "->...->" << reverse_complementary_node(node1) << "->" << out_rev[0] << "; path2: " << in_rev[1] << "->" << reverse_complementary_node(node2) << "->...->" << reverse_complementary_node(node1) << "->" << out_rev[1] << std::endl;


        Edge edge1(path_1_r.sequence.at(k), path_1_r.length, path_1_r.sequence, path_1_r.multiplicity);
        edge1.path_edges_in_original_graph = path_1_r.path_edges_in_original_graph;
        edge1.path_nodes_in_original_graph = path_1_r.path_nodes_in_original_graph;

        this->graph[in_rev[0]].outgoing_edges[out_rev[0]].push_back(edge1);
        this->graph[out_rev[0]].incoming_edges[in_rev[0]].push_back(edge1);

        Edge edge2(path_2_r.sequence.at(k), path_2_r.length, path_2_r.sequence, path_2_r.multiplicity);
        edge2.path_edges_in_original_graph = path_2_r.path_edges_in_original_graph;
        edge2.path_nodes_in_original_graph = path_2_r.path_nodes_in_original_graph;
        this->graph[in_rev[1]].outgoing_edges[out_rev[1]].push_back(edge2);
        this->graph[out_rev[1]].incoming_edges[in_rev[1]].push_back(edge2);

        for (int i = 0; i < path_1_r.nodes.size() - 1; ++i) {
            this->graph[path_1_r.nodes[i]].outgoing_edges.erase(path_1_r.nodes[i + 1]);
            this->graph[path_1_r.nodes[i + 1]].incoming_edges.erase(path_1_r.nodes[i]);
        }

        for (int i = 0; i < path_2_r.nodes.size() - 1; ++i) {
            this->graph[path_2_r.nodes[i]].outgoing_edges.erase(path_2_r.nodes[i + 1]);
            this->graph[path_2_r.nodes[i + 1]].incoming_edges.erase(path_2_r.nodes[i]);
        }

        for (int i = 0; i < path_1_r.nodes.size(); ++i) {
            if (graph[path_1_r.nodes[i]].incoming_edges.empty() && graph[path_1_r.nodes[i]].outgoing_edges.empty() && std::find(nodes_to_remove.begin(), nodes_to_remove.end(), path_1_r.nodes[i]) == nodes_to_remove.end())
                nodes_to_remove.push_back(path_1_r.nodes[i]);
        }
        for (int i = 0; i < path_2_r.nodes.size(); ++i) {
            if (graph[path_2_r.nodes[i]].incoming_edges.empty() && graph[path_2_r.nodes[i]].outgoing_edges.empty() && std::find(nodes_to_remove.begin(), nodes_to_remove.end(), path_2_r.nodes[i]) == nodes_to_remove.end())
                nodes_to_remove.push_back(path_2_r.nodes[i]);
        }

    }

    unsigned bulges = 1;
    multi_bulge_removal(bulges, false);

    decoupled_stands += 2;
    // write_graph("debugging_" + path1.nodes.at(0) + "_" + path1.nodes.at(path1.nodes.size() - 1) + "_after");
}

void Graph::decoupling(int& decoupled_stands, bool strict) {
    unsigned removed_paths = 1;
    unsigned bulges = 1;
    if (!strict) {
        bulges = 1;
        while (bulges)
            multi_bulge_removal(bulges, false);
    }
    while (removed_paths) {
        resolving_bulge_with_two_multi_edge_paths(removed_paths, 5, 0, true, 2, !strict);
    }
    removed_paths = 1;
    while (removed_paths)
        multi_bulge_removal(removed_paths, false);

    decoupled_stands = 0;
    std::vector<std::string> nodes_to_remove;
    std::vector<std::string> source_nodes, sink_nodes;
    //search for 2-in-2-out edge
    for (auto&& node : graph) {
        if (node.second.incoming_edges.size() == 2 && node.second.outgoing_edges.size() == 1) {
            std::string node_sink = node.first;
            while (graph[node_sink].outgoing_edges.size() == 1) {
                bool flag = false;
                for (auto&& n : graph[node_sink].outgoing_edges) {
                    if (n.first != node_sink) {
                        node_sink = n.first;
                        flag = true;
                    }
                }
                if (flag == false)
                    break;
                if (graph[node_sink].incoming_edges.size() != 1)
                    break;
            }
            // find a 2-in-2-out component
            if (this->graph[node_sink].incoming_edges.size() == 1 && this->graph[node_sink].outgoing_edges.size() == 2 && graph[node_sink].outgoing_edges.find(node_sink) == graph[node_sink].outgoing_edges.end()) {
                source_nodes.push_back(node.first);
                sink_nodes.push_back(node_sink);
            }
        }
    }

    for (int i = 0; i < source_nodes.size();++i)
        this->decouple_strands(source_nodes[i], sink_nodes[i], decoupled_stands, nodes_to_remove, strict);

    for (auto&& node : nodes_to_remove) {
        this->graph.erase(node);
    }
    nodes_to_remove.clear();
    // write_graph("debug" + std::to_string(cnt_rounds) + "before_non");
    this->merge_non_branching_paths(true);
    // write_graph("debug" + std::to_string(cnt_rounds) + "end");
}


void remove_items_from_vector(std::vector<int>& vec, const std::vector<int>& discontinued_indices) {
    // Make a copy of the indices and sort them in descending order
    std::vector<int> sorted_indices = discontinued_indices;
    std::sort(sorted_indices.rbegin(), sorted_indices.rend());

    // Remove elements at the specified indices from the vector
    for (int index : sorted_indices) {
        if (index >= 0 && index < vec.size()) {
            vec.erase(vec.begin() + index);
        }
    }
}
template<typename T>
void Graph::remove_items_from_vector(std::vector<T>& vec, const std::vector<int>& discontinued_indices) {
    // Make a copy of the indices and sort them in descending order
    std::vector<int> sorted_indices = discontinued_indices;
    std::sort(sorted_indices.rbegin(), sorted_indices.rend());

    // Remove elements at the specified indices from the vector
    for (int index : sorted_indices) {
        if (index >= 0 && index < vec.size()) {
            vec.erase(vec.begin() + index);
        }
    }
}

void Graph::remove_low_cov_on_node(std::string node, unsigned removed_edges, double coverage, std::vector<std::string>& nodes_to_remove) {
    std::vector<std::string> out_to_remove;
    for (auto&& n2 : graph[node].outgoing_edges) {
        std::vector<int> indices;
        for (int i = 0; i < n2.second.size(); ++i) {
            if (n2.second.at(i).multiplicity <= 10 || (n2.second.at(i).multiplicity < coverage && n2.second.at(i).length < 10000 && n2.first != reverse_complementary_node(node))) {
                indices.push_back(i);
                removed_edges += 1;
            }
        }
        remove_items_from_vector(n2.second, indices);
        if (n2.second.empty())
            out_to_remove.push_back(n2.first);

        indices.clear();
        for (int i = 0; i < graph[n2.first].incoming_edges[node].size(); ++i) {
            if (graph[n2.first].incoming_edges[node].at(i).multiplicity <= 10 || (graph[n2.first].incoming_edges[node].at(i).multiplicity < coverage && graph[n2.first].incoming_edges[node].at(i).length < 10000 && n2.first != reverse_complementary_node(node))) {
                indices.push_back(i);
            }
        }
        remove_items_from_vector(graph[n2.first].incoming_edges[node], indices);
        if (graph[n2.first].incoming_edges[node].empty())
            graph[n2.first].incoming_edges.erase(node);
    }

    std::vector<std::string> in_to_remove;
    for (auto&& n2 : graph[node].incoming_edges) {
        std::vector<int> indices;
        for (int i = 0; i < n2.second.size(); ++i) {
            if (n2.second.at(i).multiplicity <= 10 || (n2.second.at(i).multiplicity < coverage && n2.second.at(i).length < 10000 && n2.first != reverse_complementary_node(node))) {
                indices.push_back(i);
                removed_edges += 1;
            }
        }
        remove_items_from_vector(n2.second, indices);
        if (n2.second.empty())
            in_to_remove.push_back(n2.first);

        indices.clear();
        for (int i = 0; i < graph[n2.first].outgoing_edges[node].size(); ++i) {
            if (graph[n2.first].outgoing_edges[node].at(i).multiplicity <= 10 || (graph[n2.first].outgoing_edges[node].at(i).multiplicity < coverage && graph[n2.first].outgoing_edges[node].at(i).length < 10000 && n2.first != reverse_complementary_node(node))) {
                indices.push_back(i);
            }
        }
        remove_items_from_vector(graph[n2.first].outgoing_edges[node], indices);
        if (graph[n2.first].outgoing_edges[node].empty())
            graph[n2.first].outgoing_edges.erase(node);
    }

    for (auto&& e : out_to_remove) {
        graph[node].outgoing_edges.erase(e);
    }
    for (auto&& e : in_to_remove) {
        graph[node].incoming_edges.erase(e);
    }
    if (graph[node].incoming_edges.empty() && graph[node].outgoing_edges.empty())
        nodes_to_remove.push_back(node);

}

void Graph::remove_low_coverage_edges(unsigned& removed_edges, double coverage, bool tips) {
    removed_edges = 0;
    std::vector<std::string> nodes_to_remove;
    std::unordered_set<std::string> scanned_nodes;
    for (auto&& node : graph) {
        if (scanned_nodes.find(node.first) != scanned_nodes.end())
            continue;
        if (tips == true && node.second.incoming_edges.size() >= 1 && node.second.outgoing_edges.size() >= 1)
            continue;
        scanned_nodes.insert(node.first);
        scanned_nodes.insert(reverse_complementary_node(node.first));
        remove_low_cov_on_node(node.first, removed_edges, coverage, nodes_to_remove);
        remove_low_cov_on_node(reverse_complementary_node(node.first), removed_edges, coverage, nodes_to_remove);
        if (node.second.outgoing_edges.size() == 1 && node.second.incoming_edges.size() == 1 && node.second.outgoing_edges.find(node.first) != node.second.outgoing_edges.end()) {
            nodes_to_remove.push_back(node.first);
            nodes_to_remove.push_back(reverse_complementary_node(node.first));
        }
    }
    for (auto&& n : nodes_to_remove)
        graph.erase(n);

    merge_non_branching_paths(true);
}
