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