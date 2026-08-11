#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "dot_graph.hpp"
#include <cassert>
#include <set>
#include <queue>
#include <filesystem>
namespace fs = std::filesystem;

using namespace dbg;

int Graph::write_reads(const std::string& prefix) {
    int max_length = 0;
    std::ofstream fout(prefix + ".fasta");
    for (auto&& r : read2aln) {
        if (r.second.nodes_path.size() < 2)
            continue;
        std::vector<std::string> nodes;
        std::string sequence = find_path_from_start_bases(r.second.nodes_path[0], r.second.start_base_path, r.first, nodes, r.second.prefix, r.second.suffix);
        // std::cout << "Length for " << r.first << ": " << sequence.size() << std::endl;
        fout << ">" << r.first << "\n" << sequence << "\n";
        if (sequence.size() > max_length)
            max_length = sequence.size();
    }
    fout.close();
    return max_length;
}

void Graph::load_read_path(const std::string& graph_aln) {
    if (!fs::is_regular_file(graph_aln))
        return;
    std::ifstream aln_file(graph_aln);
    if (graph_aln.find("dbg") != std::string::npos) {
        std::cout << "[ReadGraph] Load corrected and pseudo reads from " << graph_aln << std::endl;
        std::string line;
        getline(aln_file, line);
        int cnt_frames = std::atoi(line.c_str());
        assert(cnt_frames == 2);
        getline(aln_file, line);
        int cnt_corrected = std::atoi(line.c_str());;
        for (int i = 0; i < cnt_corrected; ++i) {
            getline(aln_file, line);
            std::size_t pos = 0;
            std::vector<std::string> items;
            while (pos != std::string::npos) {
                std::size_t pos_next = line.find(' ', pos);
                items.push_back(line.substr(pos, pos_next - pos));
                if (pos_next != std::string::npos)
                    pos = pos_next + 1;
                else
                    pos = pos_next;
            }
            read2aln[items[0]].prefix = std::atoi(items[3].c_str());
            read2aln[items[0]].suffix = std::atoi(items[4].c_str());
            read2aln[items[0]].start_base_path = items[2].substr(2);
            find_path_from_start_bases(items[1], read2aln[items[0]].start_base_path, items[0], read2aln[items[0]].nodes_path, read2aln[items[0]].prefix, read2aln[items[0]].suffix);
            if (items[1] != "0")
                assert(read2aln[items[0]].start_base_path.size() + 1 == read2aln[items[0]].nodes_path.size());
        }
        getline(aln_file, line);
        int cnt_pseudo = std::atoi(line.c_str());;
        for (int i = 0; i < cnt_pseudo; ++i) {
            getline(aln_file, line);
            std::size_t pos = 0;
            std::vector<std::string> items;
            while (pos != std::string::npos) {
                std::size_t pos_next = line.find(' ', pos);
                items.push_back(line.substr(pos, pos_next - pos));
                if (pos_next != std::string::npos)
                    pos = pos_next + 1;
                else
                    pos = pos_next;
            }
            pseudo2aln[items[0]].prefix = std::atoi(items[3].c_str());
            pseudo2aln[items[0]].suffix = std::atoi(items[4].c_str());
            pseudo2aln[items[0]].start_base_path = items[2].substr(2);
            find_path_from_start_bases(items[1], pseudo2aln[items[0]].start_base_path, items[0], pseudo2aln[items[0]].nodes_path, pseudo2aln[items[0]].prefix, pseudo2aln[items[0]].suffix);
            if (items[1] != "0")
                assert(pseudo2aln[items[0]].start_base_path.size() + 1 == pseudo2aln[items[0]].nodes_path.size());
        }
        std::cout << "[ReadGraph] Load " << cnt_corrected << " corrected reads, " << cnt_pseudo << " pseudo reads." << std::endl;
    }
    else {
        std::cout << "[ReadGraph] Load corrected reads from " << graph_aln << std::endl;
        std::string line;
        getline(aln_file, line);
        int cnt_corrected = std::atoi(line.c_str());;
        for (int i = 0; i < cnt_corrected; ++i) {
            getline(aln_file, line);
            std::size_t pos = 0;
            std::vector<std::string> items;
            while (pos != std::string::npos) {
                std::size_t pos_next = line.find(' ', pos);
                items.push_back(line.substr(pos, pos_next - pos));
                if (pos_next != std::string::npos)
                    pos = pos_next + 1;
                else
                    pos = pos_next;
            }
            read2aln[items[0]].prefix = std::atoi(items[3].c_str());
            read2aln[items[0]].suffix = std::atoi(items[4].c_str());
            read2aln[items[0]].start_base_path = items[2].substr(2);
            find_path_from_start_bases(items[1], read2aln[items[0]].start_base_path, items[0], read2aln[items[0]].nodes_path, read2aln[items[0]].prefix, read2aln[items[0]].suffix);
            if (items[1] != "0")
                assert(read2aln[items[0]].start_base_path.size() + 1 == read2aln[items[0]].nodes_path.size());
        }
        std::cout << "[ReadGraph] Load " << cnt_corrected << " corrected reads" << std::endl;
    }
}

std::string Graph::find_path_from_start_bases(std::string start_node, std::string start_bases, std::string read_name, std::vector<std::string>& nodes_path, int prefix, int suffix) {
    if (start_node == "0")
        return "";

    std::string seq = "";
    if (graph.find(start_node) == graph.end()) {
        std::cout << "[ReadGraph] No start node found on graph for " << start_node << ", P:" << start_bases << std::endl;
    }
    assert(graph.find(start_node) != graph.end());

    int index = 0;
    nodes_path.push_back(start_node);
    while (index != start_bases.size()) {
        bool flag = false;
        std::string selected_next;
        for (auto&& node_next : graph[start_node].outgoing_edges) {
            for (auto&& e : node_next.second) {
                if (e.start_base == start_bases[index]) {
                    if (flag == true) {
                        std::cout << "[ReadGraph] Multiple paths found for " << start_node << ", P:" << start_bases << " at index " << index << std::endl;
                        continue;
                    }
                    e.reads.insert(read_name);
                    nodes_path.push_back(node_next.first);
                    selected_next = node_next.first;
                    flag = true;

                    if (index == 0 && index != start_bases.size() - 1) {
                        assert(prefix < e.sequence.size());
                        seq += e.sequence.substr(prefix);
                    }
                    else if (index == 0 && index == start_bases.size() - 1) {
                        assert(e.sequence.size() > prefix + suffix);
                        seq += e.sequence.substr(prefix, e.sequence.size() - prefix - suffix);
                    }
                    else if (index != 0 && index != start_bases.size() - 1) {
                        seq += e.sequence.substr(k);
                    }
                    else {
                        assert(e.sequence.size() - suffix - k > 0);
                        seq += e.sequence.substr(k, e.sequence.size() - suffix - k);
                    }
                }
            }
        }
        if (flag == false) {
            throw std::runtime_error("No path found for " + start_node + ", P:" + start_bases + " at index " + std::to_string(index));
        }
        index += 1;
        start_node = selected_next;
    }
    return seq;
}

void Graph::reroute_reads_from_edge_to_edge(Edge& edge_des, Edge& edge_ori, std::string node_s, std::string node_e) {
    edge_des.reads.insert(edge_ori.reads.begin(), edge_ori.reads.end());
    for (auto&& r : edge_ori.reads) {
        ReadAln& aln = read2aln.find(r) != read2aln.end() ? read2aln[r] : pseudo2aln[r];
        for (int i = 0; i < aln.start_base_path.size();++i) {
            if (aln.nodes_path[i] == node_s && aln.nodes_path[i + 1] == node_e && aln.start_base_path[i] == edge_ori.start_base) {
                // update start base to new edge
                aln.start_base_path[i] = edge_des.start_base;
                // if it's in the start/end, need update the prefix/suffix
                if (i == 0) {
                    aln.prefix = int(1.0 * aln.prefix / edge_ori.length * edge_des.length);
                }
                if (i == aln.start_base_path.size() - 1) {
                    aln.suffix = int(1.0 * aln.suffix / edge_ori.length * edge_des.length);
                }
                if (i == 0 && i == aln.start_base_path.size() - 1)
                    assert(aln.prefix + aln.suffix <= edge_des.length);
            }
        }
    }
}

void Graph::reroute_reads_from_path_to_edge(Edge& edge_des, std::string node_s, std::string node_m, std::string node_e) {
    Edge& edge_ori_1 = graph[node_s].outgoing_edges[node_m].at(0);
    Edge& edge_ori_2 = graph[node_m].outgoing_edges[node_e].at(0);

    edge_des.reads.insert(edge_ori_1.reads.begin(), edge_ori_1.reads.end());
    edge_des.reads.insert(edge_ori_2.reads.begin(), edge_ori_2.reads.end());

    for (auto&& r : edge_des.reads) {
        ReadAln& aln = read2aln.find(r) != read2aln.end() ? read2aln[r] : pseudo2aln[r];
        int i = 1;
        while (i + 1 < aln.nodes_path.size()) {
            if (aln.nodes_path[i] == node_m) {
                assert(aln.nodes_path[i - 1] == node_s);
                assert(aln.nodes_path[i + 1] == node_e);
                aln.nodes_path.erase(aln.nodes_path.begin() + i);
                aln.start_base_path.erase(aln.start_base_path.begin() + i);
            }
            else {
                i += 1;
            }
        }
        if (aln.nodes_path.size() > 0 && aln.nodes_path[0] == node_m) {
            if (aln.nodes_path.size() > 1)
                assert(aln.nodes_path[1] == node_e);
            aln.nodes_path[0] = node_s;
            aln.start_base_path[0] = edge_ori_1.start_base;
            aln.prefix += edge_ori_1.length;
        }

        if (aln.nodes_path.size() > 0 && aln.nodes_path[aln.nodes_path.size() - 1] == node_m) {
            if (aln.nodes_path.size() > 1)
                assert(aln.nodes_path[aln.nodes_path.size() - 2] == node_s);
            aln.nodes_path[aln.nodes_path.size() - 1] = node_e;
            aln.suffix += edge_ori_2.length;
        }
    }
}

void Graph::reroute_reads_from_outtip_to_edge(Edge& edge_des, std::string node_s, std::string node_e, std::string node_t) {
    // tip from node_s to node_t
    Edge& tip = graph[node_s].outgoing_edges[node_t][0];
    edge_des.reads.insert(tip.reads.begin(), tip.reads.end());

    for (auto&& r : tip.reads) {
        // std::cout << "Process for outtip " << node_s << " to " << node_t << " read " << r << std::endl;
        ReadAln& aln = read2aln.find(r) != read2aln.end() ? read2aln[r] : pseudo2aln[r];
        if (aln.nodes_path.size() >= 2) {
            assert(aln.nodes_path[aln.nodes_path.size() - 1] == node_t &&
                aln.nodes_path[aln.nodes_path.size() - 2] == node_s &&
                aln.start_base_path[aln.start_base_path.size() - 1] == tip.start_base
            );
            aln.nodes_path[aln.nodes_path.size() - 1] = node_e;
            aln.start_base_path[aln.start_base_path.size() - 1] = graph[node_s].outgoing_edges[node_e].at(0).start_base;
            if (aln.nodes_path.size() == 2 && aln.prefix >= edge_des.length) {
                aln.nodes_path.clear();
                aln.prefix = 0;
                aln.suffix = 0;
                aln.start_base_path = "";
                // std::cout << "Prefix too long, the alignment is removed." << std::endl;
            }
            else if (edge_des.length > (tip.length - aln.suffix)) {
                // std::cout << "Suffix updated from " << aln.suffix;
                aln.suffix = edge_des.length - (tip.length - aln.suffix);
                // std::cout << " to " << aln.suffix << ", tip " << tip.length << " edge " << edge_des.length << std::endl;
            }
            else {
                aln.suffix = 0;
                // std::cout << "Suffix updated to 0, tip " << tip.length << " edge " << edge_des.length << " previous suffix " << aln.suffix << std::endl;
            }
        }
    }
}

void Graph::reroute_reads_from_intip_to_edge(Edge& edge_des, std::string node_s, std::string node_e, std::string node_t) {
    Edge& tip = graph[node_t].outgoing_edges[node_e][0];
    edge_des.reads.insert(tip.reads.begin(), tip.reads.end());

    for (auto&& r : tip.reads) {
        // std::cout << "Process for intip " << node_t << " to " << node_e << " read " << r << std::endl;
        ReadAln& aln = read2aln.find(r) != read2aln.end() ? read2aln[r] : pseudo2aln[r];
        if (aln.nodes_path.size() >= 2) {
            assert(aln.nodes_path[0] == node_t &&
                aln.nodes_path[1] == node_e &&
                aln.start_base_path[0] == tip.start_base);
            aln.nodes_path[0] = node_s;
            aln.start_base_path[0] = graph[node_s].outgoing_edges[node_e].at(0).start_base;
            if (aln.nodes_path.size() == 2 && aln.suffix >= edge_des.length) {
                aln.nodes_path.clear();
                aln.prefix = 0;
                aln.suffix = 0;
                aln.start_base_path = "";
                // std::cout << "Suffix too long, the alignment is removed." << std::endl;
            }
            else if (edge_des.length > (tip.length - aln.prefix)) {
                // std::cout << "Prefix updated from " << aln.prefix;
                aln.prefix = edge_des.length - (tip.length - aln.prefix);
                // std::cout << " to " << aln.prefix << ", tip " << tip.length << " edge " << edge_des.length << std::endl;
            }
            else {
                aln.prefix = 0;
                // std::cout << "Prefix updated to 0, tip " << tip.length << " edge " << edge_des.length << " previous prefix " << aln.prefix << std::endl;
            }
        }
    }
}

void Graph::add_virtual_reads(unsigned& removed_paths, int x, double identity, bool use_length) {
    removed_paths = 0;
    std::vector<std::string> nodes_to_remove;
    std::unordered_set<std::string> traversed_bulges;
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

                    // check reverse paths
                    Path p1_reverse, p2_reverse;
                    get_reverse_path(paths[i], p1_reverse);
                    get_reverse_path(paths[j], p2_reverse);

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

            Bulge b(p1, p2, max_identity);
            std::string reads = p1.nodes.at(0) + "_" + p1.nodes.at(p1.nodes.size() - 1);
            std::string reads_rev = reverse_complementary_node(p1.nodes.at(p1.nodes.size() - 1)) + "_" + reverse_complementary_node(p1.nodes.at(0));
            if (traversed_bulges.find(reads) != traversed_bulges.end())
                continue;
            traversed_bulges.insert(reads);
            traversed_bulges.insert(reads_rev);
            if (p1.nodes.size() > 2) {
                pseudo2aln[reads + "_p1"].nodes_path = p1.nodes;
                pseudo2aln[reads + "_p1"].prefix = 0;
                pseudo2aln[reads + "_p1"].suffix = 0;
                for (int i = 0; i < p1.nodes.size() - 1; ++i) {
                    pseudo2aln[reads + "_p1"].start_base_path += graph[p1.nodes.at(i)].outgoing_edges[p1.nodes.at(i + 1)].at(p1.bulge_legs.at(i)).start_base;
                    graph[p1.nodes.at(i)].outgoing_edges[p1.nodes.at(i + 1)][p1.bulge_legs.at(i)].reads.insert(reads + "_p1");
                }
                removed_paths += 1;
            }
            if (p2.nodes.size() > 2) {
                pseudo2aln[reads + "_p2"].nodes_path = p2.nodes;
                pseudo2aln[reads + "_p2"].prefix = 0;
                pseudo2aln[reads + "_p2"].suffix = 0;
                for (int i = 0; i < p2.nodes.size() - 1; ++i) {
                    pseudo2aln[reads + "_p2"].start_base_path += graph[p2.nodes.at(i)].outgoing_edges[p2.nodes.at(i + 1)].at(p2.bulge_legs.at(i)).start_base;
                    graph[p2.nodes.at(i)].outgoing_edges[p2.nodes.at(i + 1)][p2.bulge_legs.at(i)].reads.insert(reads + "_p2");
                }
                removed_paths += 1;
            }
        }
    }
}

void Graph::add_complementary_virtual_reads(unsigned& added_reads) {
    added_reads = 0;
    std::vector<std::string> source_nodes, sink_nodes;
    std::unordered_set<std::string> traversed_nodes;
    //search for 2-in-2-out edge
    for (auto&& node : graph) {
        if (traversed_nodes.find(node.first) != traversed_nodes.end())
            continue;
        if (node.second.incoming_edges.size() == 2 && node.second.outgoing_edges.size() == 1) {
            std::string node_sink = node.first;
            if (graph[node_sink].outgoing_edges.size() == 1) {
                bool flag = false;
                for (auto&& n : graph[node_sink].outgoing_edges) {
                    if (n.first != node_sink) {
                        node_sink = n.first;
                        flag = true;
                    }
                }
                if (flag == false)
                    continue;
                if (graph[node_sink].incoming_edges.size() != 1)
                    continue;
            }
            // find a 2-in-2-out component
            if (this->graph[node_sink].incoming_edges.size() == 1 && this->graph[node_sink].outgoing_edges.size() == 2 && graph[node_sink].outgoing_edges.find(node_sink) == graph[node_sink].outgoing_edges.end()) {
                source_nodes.push_back(node.first);
                sink_nodes.push_back(node_sink);
                traversed_nodes.insert(node.first);
                traversed_nodes.insert(node_sink);
                traversed_nodes.insert(reverse_complementary_node(node_sink));
                traversed_nodes.insert(reverse_complementary_node(node.first));
            }
        }
    }

    for (int i = 0; i < source_nodes.size(); ++i) {
        std::string node1 = source_nodes[i], node2 = sink_nodes[i];
        std::vector<std::string> incoming_nodes, outgoing_nodes;
        if (graph.find(node1) == graph.end() || graph.find(node2) == graph.end())
            continue;
        for (auto&& e : this->graph[node1].incoming_edges) {
            incoming_nodes.push_back(e.first);
        }
        for (auto&& e : this->graph[node2].outgoing_edges) {
            outgoing_nodes.push_back(e.first);
        }
        if (incoming_nodes.size() != 2 || outgoing_nodes.size() != 2)
            continue;

        // avoid any simple bulges
        if (graph[incoming_nodes[0]].outgoing_edges[node1].size() > 1)
            continue;
        if (graph[incoming_nodes[1]].outgoing_edges[node1].size() > 1)
            continue;
        if (graph[node1].outgoing_edges[node2].size() > 1)
            continue;
        if (graph[node2].outgoing_edges[outgoing_nodes[0]].size() > 1)
            continue;
        if (graph[node2].outgoing_edges[outgoing_nodes[1]].size() > 1)
            continue;

        std::vector<std::string> starts, ends;
        for (auto&& r : graph[node1].outgoing_edges[node2][0].reads) {
            if (pseudo2aln.find(r) == pseudo2aln.end())
                continue;
            ReadAln& aln = pseudo2aln[r];
            for (int i = 1; i < aln.nodes_path.size() - 2; ++i) {
                if (aln.nodes_path[i] == node1) {
                    assert(aln.nodes_path[i + 1] == node2);
                    starts.push_back(aln.nodes_path[i - 1]);
                    ends.push_back(aln.nodes_path[i + 2]);
                    assert(incoming_nodes[0] == aln.nodes_path[i - 1] || incoming_nodes[1] == aln.nodes_path[i - 1]);
                }
            }
        }
        if (node1 != reverse_complementary_node(node2)) {
            for (auto&& r : graph[reverse_complementary_node(node2)].outgoing_edges[reverse_complementary_node(node1)][0].reads) {
                if (pseudo2aln.find(r) == pseudo2aln.end())
                    continue;
                ReadAln& aln = pseudo2aln[r];
                for (int i = 1; i < aln.nodes_path.size() - 2; ++i) {
                    if (aln.nodes_path[i] == reverse_complementary_node(node2)) {
                        assert(aln.nodes_path[i + 1] == reverse_complementary_node(node1));
                        starts.push_back(reverse_complementary_node(aln.nodes_path[i + 2]));
                        ends.push_back(reverse_complementary_node(aln.nodes_path[i - 1]));
                        assert(incoming_nodes[0] == reverse_complementary_node(aln.nodes_path[i + 2]) || incoming_nodes[1] == reverse_complementary_node(aln.nodes_path[i + 2]));
                    }
                }
            }
        }

        std::cout << "[VirtualRead] Virtual read traversing 2-in-2-out: " << starts.size() << "(" << node1 << " -> " << node2 << ")" << std::endl;

        if (starts.size() == 1) {
            std::string virtual_start, virtual_end;
            bool flag = false;
            for (auto&& s : incoming_nodes) {
                if (s != starts[0]) {
                    assert(!flag);
                    virtual_start = s;
                    flag = true;
                }
            }
            flag = false;
            for (auto&& e : outgoing_nodes) {
                if (e != ends[0]) {
                    assert(!flag);
                    virtual_end = e;
                    flag = true;
                }
            }
            assert(flag);

            // add virtual read
            std::string read_name = node1 + "_" + node2 + "_complement";
            assert(pseudo2aln.find(read_name) == pseudo2aln.end());
            pseudo2aln[read_name].nodes_path.push_back(virtual_start);
            pseudo2aln[read_name].nodes_path.push_back(node1);
            pseudo2aln[read_name].nodes_path.push_back(node2);
            pseudo2aln[read_name].nodes_path.push_back(virtual_end);

            pseudo2aln[read_name].prefix = 0;
            pseudo2aln[read_name].suffix = 0;
            pseudo2aln[read_name].start_base_path = pseudo2aln[read_name].start_base_path + graph[virtual_start].outgoing_edges[node1].at(0).start_base +
                graph[node1].outgoing_edges[node2].at(0).start_base +
                graph[node2].outgoing_edges[virtual_end].at(0).start_base;

            graph[virtual_start].outgoing_edges[node1][0].reads.insert(read_name);
            graph[node1].outgoing_edges[node2][0].reads.insert(read_name);
            graph[node2].outgoing_edges[virtual_end][0].reads.insert(read_name);

            std::cout << "[VirtualRead] Virtual read add: " << virtual_start << " -> " << node1 << " -> " << node2 << " -> " << virtual_end << " " << pseudo2aln[read_name].start_base_path << std::endl;
            added_reads += 1;
        }
    }
}

void Graph::detect_chimeric_reads() {
    std::unordered_set<std::string> chimeric_reads;

    struct ErrorEdge
    {
        std::string node1, node2;
        std::unordered_set<std::string> chimeric_reads;
        int leg_forward, leg_reverse;

        std::string reverse_complementary(std::string& seq) {
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

        std::string get_new_node(std::unordered_map<std::string, Node>& graph, std::string seq) {
            std::string n;
            if (node2.at(0) == '-')
                n = node2.substr(1);
            else
                n = node2;
            n = n + "1";
            while (graph.find(n) != graph.end()) {
                n += "1";
            }
            if (seq <= reverse_complementary(seq))
                return n;
            else
                return "-" + n;
        }

        ErrorEdge(std::string n1, std::string n2, std::unordered_set<std::string>& chim_r, int edge, int edge_r) {
            node1 = n1;
            node2 = n2;
            chimeric_reads = chim_r;
            leg_forward = edge;
            leg_reverse = edge_r;
        }
    };

    std::vector<ErrorEdge> erroneous_edges;

    for (auto&& node1 : graph) {
        for (auto&& node2 : node1.second.outgoing_edges) {
            for (int i = 0; i < node2.second.size(); ++i) {
                auto& e1 = node2.second[i];
                if (graph[node2.first].outgoing_edges.size() != 1)
                    continue;
                for (auto&& node3 : graph[node2.first].outgoing_edges) {
                    if (node3.second.size() > 1)
                        continue;
                    for (int j = 0;j < node3.second.size(); ++j) {
                        if (node1.first == node2.first && node2.first == node3.first)
                            continue;
                        auto& e2 = node3.second[j];
                        Path p;
                        add_node_to_path(p, node1.first, 0);
                        add_node_to_path(p, node2.first, i);
                        add_node_to_path(p, node3.first, j);

                        Path p_reverse;
                        get_reverse_path(p, p_reverse);

                        assert(p_reverse.nodes.size() == 3);
                        auto& e1_r = graph[p_reverse.nodes[1]].outgoing_edges[p_reverse.nodes[2]][p_reverse.bulge_legs[1]];
                        auto& e2_r = graph[p_reverse.nodes[0]].outgoing_edges[p_reverse.nodes[1]][p_reverse.bulge_legs[0]];

                        std::unordered_set<std::string> e1_all = e1.reads;
                        e1_all.insert(e1_r.reads.begin(), e1_r.reads.end());

                        std::unordered_set<std::string> e2_all = e2.reads;
                        e2_all.insert(e2_r.reads.begin(), e2_r.reads.end());

                        auto inters = get_intersection(e1_all, e2_all);

                        // if (node1.first == "-33853")
                        //     std::cout << node1.first << " -> (" << e1.length << " " << e1.multiplicity << ") -> " << node2.first << " -> (" << e2.length << " " << e2.multiplicity << ") -> " << node3.first << " reads1 " << e1_all.size() << " reads2 " << e2_all.size() << " shared " << inters.size() << ":" << std::endl;

                        // if (inters.size() <= 1 && (e1_all.size() >= 10 || e2_all.size() >= 10)) {
                        if (inters.size() <= 1) {
                            erroneous_edges.emplace_back(ErrorEdge(node1.first, node2.first, inters, i, p_reverse.bulge_legs[1]));

                            std::cout << "[RemoveChimeric] Find erroneuous edge " << node1.first << " -> (" << e1.length << " " << e1.multiplicity << ") -> " << node2.first << " -> (" << e2.length << " " << e2.multiplicity << ") -> " << node3.first << " reads1 " << e1_all.size() << " reads2 " << e2_all.size() << " shared " << inters.size() << ":" << std::endl;
                            for (auto r : inters) {
                                if (chimeric_reads.find(r) == chimeric_reads.end())
                                    std::cout << "    Find chimeric read " << r << ", " << node1.first << " (" << e1.length << " " << e1.multiplicity << ") " << node2.first << " (" << e2.length << " " << e2.multiplicity << ") " << node3.first << std::endl;
                                chimeric_reads.insert(r);
                            }
                        }
                    }
                }
            }

        }
    }

    for (auto&& error_edge : erroneous_edges) {

        auto& out_forward = graph[error_edge.node1].outgoing_edges[error_edge.node2][error_edge.leg_forward];
        auto& in_forward = graph[error_edge.node2].incoming_edges[error_edge.node1][error_edge.leg_forward];
        auto& out_reverse = graph[reverse_complementary_node(error_edge.node2)].outgoing_edges[reverse_complementary_node(error_edge.node1)][error_edge.leg_reverse];
        auto& in_reverse = graph[reverse_complementary_node(error_edge.node1)].incoming_edges[reverse_complementary_node(error_edge.node2)][error_edge.leg_reverse];

        if (out_forward.multiplicity == 0)
            continue;

        for (auto&& r : error_edge.chimeric_reads) {
            ReadAln& aln = read2aln.find(r) != read2aln.end() ? read2aln[r] : pseudo2aln[r];
            aln.prefix = 0;
            aln.suffix = 0;
            aln.start_base_path = "";
            aln.nodes_path.clear();

            out_forward.reads.erase(r);
            in_forward.reads.erase(r);
            out_reverse.reads.erase(r);
            in_reverse.reads.erase(r);
        }

        out_forward.sequence = out_forward.sequence.substr(0, out_forward.sequence.size() - k);
        in_forward.sequence = in_forward.sequence.substr(0, in_forward.sequence.size() - k);
        out_reverse.sequence = out_reverse.sequence.substr(k);
        in_reverse.sequence = in_reverse.sequence.substr(k);

        assert(out_forward.sequence.size() == out_reverse.sequence.size());
        assert(out_forward.sequence.size() == in_forward.sequence.size());
        assert(in_reverse.sequence.size() == out_reverse.sequence.size());

        if (out_forward.sequence.size() > k) {
            out_forward.length -= k;
            in_forward.length -= k;
            out_reverse.length -= k;
            in_reverse.length -= k;

            out_reverse.start_base = out_reverse.sequence.at(k);
            in_reverse.start_base = in_reverse.sequence.at(k);

            assert(out_forward.length == out_forward.sequence.size() - k);
            assert(in_forward.length == in_forward.sequence.size() - k);
            assert(out_reverse.length == out_reverse.sequence.size() - k);
            assert(in_reverse.length == in_reverse.sequence.size() - k);

            out_forward.label = "";
            in_forward.label = "";
            out_reverse.label = "";
            in_reverse.label = "";
            assert(out_forward.sequence == reverse_complementary(out_reverse.sequence));

            std::string new_node = error_edge.get_new_node(graph, out_forward.sequence.substr(out_forward.sequence.size() - k));

            for (auto&& r : out_forward.reads) {
                ReadAln& aln = read2aln.find(r) != read2aln.end() ? read2aln[r] : pseudo2aln[r];
                if (aln.nodes_path.size() != 0) {
                    assert(aln.nodes_path.at(aln.nodes_path.size() - 1) == error_edge.node2);
                    aln.nodes_path[aln.nodes_path.size() - 1] = new_node;
                    assert(aln.start_base_path[aln.start_base_path.size() - 1] == out_forward.start_base);
                    aln.suffix = aln.suffix < k ? 0 : aln.suffix - k;
                }
            }

            for (auto&& r : out_reverse.reads) {
                ReadAln& aln = read2aln.find(r) != read2aln.end() ? read2aln[r] : pseudo2aln[r];
                if (aln.nodes_path.size() != 0) {
                    assert(aln.nodes_path.at(0) == reverse_complementary_node(error_edge.node2));
                    aln.nodes_path[0] = reverse_complementary_node(new_node);
                    aln.start_base_path[0] = out_reverse.start_base;
                    aln.prefix = aln.prefix < k ? 0 : aln.prefix - k;
                }
            }

            graph[error_edge.node1].outgoing_edges[new_node].push_back(out_forward);
            graph[new_node].incoming_edges[error_edge.node1].push_back(in_forward);
            // std::cout << "[RemoveChimeric] Add new edge: " << error_edge.node1 << " -> " << new_node << std::endl;
            graph[reverse_complementary_node(new_node)].outgoing_edges[reverse_complementary_node(error_edge.node1)].push_back(out_reverse);
            graph[reverse_complementary_node(error_edge.node1)].incoming_edges[reverse_complementary_node(new_node)].push_back(in_reverse);

            // std::cout << "[RemoveChimeric] Add new edge: " << reverse_complementary_node(new_node) << " -> " << reverse_complementary_node(error_edge.node1) << std::endl;
        }

        out_forward.multiplicity = 0;
        in_forward.multiplicity = 0;
        out_reverse.multiplicity = 0;
        in_reverse.multiplicity = 0;
    }
}

std::unordered_set<std::string> Graph::get_intersection(std::unordered_set<std::string>& set1, std::unordered_set<std::string>& set2) {
    std::unordered_set<std::string> result;

    // Always iterate through the smaller set for better performance
    const auto& smaller = set1.size() <= set2.size() ? set1 : set2;
    const auto& larger = set1.size() <= set2.size() ? set2 : set1;

    for (const auto& elem : smaller) {
        if (larger.find(elem) != larger.end()) {
            result.insert(elem);
        }
    }

    return result;
}