#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#define MIN_MULTI 0.1

class Path;
class Edge;

class Edge {
public:
    char start_base;
    unsigned length;
    std::string label;
    std::string rc_label;
    std::string sequence;
    std::vector<std::string> ref_ids;
    std::vector<std::string> path_nodes_in_original_graph;
    std::vector<std::string> path_edges_in_original_graph;
    double multiplicity = 0;
    void add_multi_from_edge_or_path(Edge& edge_or_path);
    void add_multi_from_edge_or_path(Path& edge_or_path);
    void remove_multi_from_path(Path& edge_or_path);
    Edge(char& base, const unsigned& length, const std::string& sequence, const double& multi);
    Edge();
};

class Node {
public:
    std::unordered_map<std::string, std::vector<Edge>> incoming_edges;
    std::unordered_map<std::string, std::vector<Edge>> outgoing_edges;
    int number_of_contracted_edge = 0;
    int length_of_contracted_edge = 0;
    double max_in_multi = 0;
    std::string max_in_node;
    template<typename T>
    void mergeMaps(std::unordered_map<std::string, std::vector<T>>& map1, const std::unordered_map<std::string, std::vector<T>>& map2);
    Node();
};

class Genome {
public:
    std::vector<std::string> genome_path;
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<int>>> node2order;
    void add_node_to_end(std::string node);
    void erase_node(std::string node);
    void construct_edge_orders();
    bool find_path(Path& path, std::vector<int>& pos);
    Genome();
};

class Path {
public:
    std::vector<std::string> nodes;
    std::vector<int> bulge_legs;
    std::vector<std::string> path_nodes_in_original_graph;
    std::vector<std::string> path_edges_in_original_graph;
    std::string sequence;
    std::unordered_map<int, std::string> index2palindromic_seq;
    unsigned length = 0;
    double multiplicity = 0;
    double min_multi = 0;
    bool is_unambiguous = true;
    bool safe_to_extract = true;

    void update_min_multi(Edge edge);
    Path();
};

class Bulge {
public:
    Path leg1, leg2;
    std::vector<Path> path_resolutions1;
    std::vector<Path> path_resolutions2;
    double max_identity = 0;
    int num_edges = 0;
    bool is_confict = false;
    Bulge();
    Bulge(Path p1, Path p2, double identity);
    void check_conflict(Bulge& bulge);
};

class Graph {
public:
    int k = 0;
    std::unordered_map<std::string, Node> graph;
    std::unordered_map<std::string, std::string> label2rc;
    void read_graph(std::string& output, std::string& restart_from, std::string& graph_dot, const std::string& graph_fasta);
    void read_from_dot(const std::string& graph_dot, const std::string& graph_fasta);
    void write_graph(const std::string& prefix, int thick = 1000000, bool contracted = false, bool colored = false, std::unordered_set<std::string> nodes = std::unordered_set<std::string>());
    void write_graph_contracted(const std::string& prefix, int min_length = 10000);
    void write_graph_colored(const std::string& prefix, const std::string& genomes);
    void write_graph_colored_from_bam(const std::string& prefix, const std::string& bam_processed);
    void decoupling(std::string multidbg, std::string output);
    int get_num_nodes();

    void multi_bulge_removal(unsigned& removed_bulges, bool skip_rc_bulges = true);
    void merge_non_branching_paths(bool merge_self_loop = false);
    void gluing_broken_bulges(unsigned& removed_bulges);
    void merge_tips_into_edges(unsigned& num_tips);
    void merge_tips(unsigned& num_tips, bool restrict_length = true, bool only_edges = false);

    void general_whirl_removal(unsigned& removed_whirls, bool simple_whirl = false, bool force = false);
    void resolving_bulge_with_two_multi_edge_paths(unsigned& removed_paths, int x, double identity, bool use_length = false, int security_level = 3, bool allow_reverse_comp = false, bool verbose = false);
    void remove_low_coverage_edges(unsigned& removed_edges, double coverage = 10, bool tips = false);
    void remove_low_cov_on_node(std::string node, unsigned removed_edges, double coverage, std::vector<std::string>& nodes_to_remove);

    void resolve_edges_in_reverse_complement(int& resolved_edges, bool strict = false);
    template<typename T>
    void merge_vecs(std::vector<T>& e1, std::vector<T>& e2);


    Graph();
private:
    std::string get_unique_label(std::unordered_set<std::string>& labels);
    bool get_reverse_path(Path& path, Path& path_reverse);
    std::string collapse_bulge(std::string node1, std::string node2, unsigned& removed_bulges, double similarity = 0);

    bool check_non_branching(std::string node, bool merge_self_loop = false);
    std::string merge_edges(std::string node, bool merge_self_loop);
    bool add_node_to_path(Path& path, std::string node, int bulge_leg = 0, bool reduce_reverse = true);

    void remove_whirl(Path& unambiguous_path, std::vector<std::string>& nodes_to_remove);

    std::string collapse_complex_bulge_two_multi_edge_paths(Path p1, Path p2, bool p1_2_in_2_out, bool p2_2_in_2_out, double max_identity, std::vector<std::string>& nodes_to_remove);
    bool process_palindromic_bulges(Path& p1, Path& p2, std::vector<std::string>& nodes_to_remove, bool verbose = false, std::string output = "");
    void find_2_in_2_out(Bulge& bulge);
    bool check_2_in_2_out(Path& leg1);

    void resolve_edges_rc(std::string node1, std::string node2, int& resolved_edges, std::vector<std::string>& nodes_to_remove, bool strict = false);
    void resolve_2_in_2_out(std::string node1, std::string node2, std::vector<std::string>& nodes_to_remove);

    int matches_by_edlib(std::string sequence1, std::string sequence2);
    int count_matches(std::string cigar);
    std::string reverse_complementary_node(std::string node);
    std::string reverse_complementary(std::string& s);
    std::string doubleToString(double value);
    template<typename T>
    void remove_items_from_vector(std::vector<T>& vec, const std::vector<int>& discontinued_indices);
};