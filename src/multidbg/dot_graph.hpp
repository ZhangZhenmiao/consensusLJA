#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>

namespace multidbg {
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
        std::string sequence;
        int number_of_contracted_edge = 0;
        int length_of_contracted_edge = 0;
        double max_in_multi = 0;
        std::string max_in_node;
        // for contracted nodes
        int circles = 0;
        long total_length = 0;
        double median_length = 0;

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
        struct Read {
            std::string name;
            int aln_start = 0;
            int eln_end = 0;
        };
        struct Path_for_MEdge {
            std::vector<std::string> dbg_nodes;
            std::vector<std::string> dbg_edges;
            std::vector<int> dbg_edge_lengths;
            std::string dbg_start_base_string;
            int dbg_path_length;
            std::string dbg_path_sequence;
            int mdbg_s_start = 0;
            int mdbg_s_end = 0;
            int mdbg_e_start = 0;
            int mdbg_e_end = 0;
            std::string mdbg_s;
            std::string mdbg_e;
            std::string mdbg_seq;
            double multi = 0;
            std::vector<Read> aln_reads;
        };
        std::unordered_map<std::string, Node> graph;
        std::unordered_map<std::string, std::string> label2rc;
        std::unordered_map<std::string, std::string> nodeid2Rev;
        std::unordered_map<std::string, std::vector<std::string>> edge_initial_to_path_in_dbg;
        std::unordered_map<std::string, Path_for_MEdge> edge_initial_to_DBG_reads;
        int error_peak = -1, first_minima = -1;
        int mean_cov = -1;
        std::map<int, int> create_histogram(const std::vector<double>& edges);
        void analyze_histogram(const std::map<int, int>& hist);
        void read_graph(std::string& output, std::string& graph_dot, const std::string& graph_fasta, const std::string& nodes_fasta, const std::string& graph_dbg, const std::string& paths_dbg);
        void read_from_dot(const std::string& graph_dot, const std::string& graph_fasta, const std::string& nodes_fasta, const std::string& graph_dbg, const std::string& paths_dbg);
        void write_graph(const std::string& prefix, int thick = 1000000, bool contracted = false, bool colored = false, std::unordered_set<std::string> nodes = std::unordered_set<std::string>());
        void write_graph_contracted(const std::string& prefix, int min_length = 10000, bool simplify = false);
        void write_graph_colored(const std::string& prefix, const std::string& genomes);
        void write_graph_colored_from_bam(const std::string& prefix, const std::string& bam_processed);
        void write_graph_gfa(const std::string& prefix);
        void decoupling(std::string multidbg, std::string output);
        void remove_contained_contigs(const double sim = 0.8);
        int get_num_nodes();
        std::string get_contracted_name(std::string node);
        std::string get_contracted_label(std::string node);

        void multi_bulge_removal(unsigned& removed_bulges, bool skip_rc_bulges = true);
        void merge_non_branching_paths(bool merge_self_loop = false);
        void gluing_broken_bulges(unsigned& removed_bulges);
        void merge_tips_into_edges(unsigned& num_tips, double ratio = 0.8);
        void merge_tips(unsigned& num_tips);

        void general_whirl_removal(unsigned& removed_whirls, bool simple_whirl = false, bool force = false);
        void resolving_bulge_with_two_multi_edge_paths(unsigned& removed_paths, int x, double identity, bool use_length = false, int security_level = 3, bool allow_reverse_comp = false, bool allow_tip = false, bool verbose = false);
        void remove_low_coverage_edges(unsigned& removed_edges, double coverage = 10, bool tips = false, bool force = false, bool avoid_chromosome = false);
        void remove_low_cov_on_node(std::string node, unsigned& removed_edges, double coverage, std::vector<std::string>& nodes_to_remove, bool force = false);
        void remove_chimeric_edge(std::string chimeric_path);

        void resolve_edges_in_reverse_complement(int& resolved_edges, bool strict = false);
        template<typename T>
        void merge_vecs(std::vector<T>& e1, std::vector<T>& e2);
        std::string getExecutablePath();
        void get_annotation(std::string prefix);

        Graph();
    private:
        std::string get_unique_label(std::unordered_set<std::string>& labels);
        bool get_reverse_path(Path& path, Path& path_reverse);
        std::string collapse_bulge(std::string node1, std::string node2, unsigned& removed_bulges, double similarity = 0);

        bool check_non_branching(std::string node, bool merge_self_loop = false);
        std::string merge_edges(std::string node, bool merge_self_loop);
        bool add_node_to_path(Path& path, std::string node, int bulge_leg = 0, bool reduce_reverse = true);

        void remove_whirl(Path& unambiguous_path, std::vector<std::string>& nodes_to_remove);

        std::string collapse_complex_bulge_two_multi_edge_paths(Path p1, Path p2, bool p1_2_in_2_out, bool p2_2_in_2_out, double max_identity, std::vector<std::string>& nodes_to_remove, bool allow_tip);
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
        void remove_items_from_vector(std::vector<int>& vec, const std::vector<int>& discontinued_indices);
        void remove_items_from_vector(std::vector<Edge>& vec, const std::vector<int>& discontinued_indices);
    };
}