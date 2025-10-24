import re
import argparse

import os
import subprocess
from pathlib import Path

EDGE_RE = re.compile(r'label="[^"]*\s+(\d+)\(')


def read_fasta(fasta_path):
    """Yield (contig_name, sequence) from a FASTA file."""
    name, seq = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            yield name, "".join(seq)


def write_single_fasta(contig_name, seq, out_dir):
    """Write one contig to a FASTA file and return the path."""
    fasta_path = Path(out_dir) / f"{contig_name}.fa"
    with open(fasta_path, "w") as f:
        f.write(f">{contig_name}\n{seq}\n")
    return fasta_path


def run_jumbodbg(contig_name, fasta_path, out_dir, threads=10):
    """Run jumboDBG on a single FASTA file and return path to graph.dot."""
    contig_out = Path(out_dir) / contig_name
    if not os.path.exists(contig_out):
        contig_out.mkdir(parents=True, exist_ok=True)
        cmd = [
            "/Poppy/zmzhang/Consensus_Assembly/cLJA/lib/LJA/bin/jumboDBG",
            "-k", "101",
            "--reads", str(fasta_path),
            "-t", str(threads),
            "--coverage",
            "-o", str(contig_out)
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return contig_out / "graph.dot"


def parse_graph_dot(dot_path):
    """Return the sum of edge lengths from a graph.dot file."""
    total = 0
    with open(dot_path) as f:
        for line in f:
            m = EDGE_RE.search(line)
            if m:
                total += int(m.group(1))
    return total/2


def compute_dbg_ratios(fasta_path, out_dir="jumbodbg_tmp", threads=10):
    """
    Run jumboDBG for each contig in a FASTA file and compute
    (sum of edge lengths / contig length) ratio.

    Returns:
        dict mapping name parts -> ratio, e.g.:
        {
            'contigA': 1.02,
            '001': 1.02,
            ...
        }
    """
    os.makedirs(out_dir, exist_ok=True)
    results = {}

    for name, seq in read_fasta(fasta_path):
        print(f"Processing {name} ...", flush=True)
        contig_len = len(seq)
        contig_fa = write_single_fasta(name, seq, out_dir)
        dot_file = run_jumbodbg(name, contig_fa, out_dir, threads)
        edge_sum = parse_graph_dot(dot_file)
        ratio = edge_sum / contig_len if contig_len > 0 else 0

        parts = name.split("_")
        if len(parts) > 0:
            results[parts[0]] = ratio
        if len(parts) > 1:
            results[parts[1]] = ratio
    
    print(results)

    return results

def parse_bam_stats(stats_file):
    """Parse BAM stats file into a dict: edge_id -> (chrom, annotation)"""
    bam_stats = {}
    with open(stats_file) as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) < 4:
                continue
            edge_id = items[0]
            chrom = items[1]
            annotation = " ".join(items[1:])

            species = "bonobo" if chrom.strip('-').startswith("bonobo") else "human"
            if edge_id not in bam_stats:
                bam_stats[edge_id] = {"chrom": species, "annotations": []}
            if species != bam_stats[edge_id]["chrom"]:
                bam_stats[edge_id]["chrom"] = "mixed"
            bam_stats[edge_id]["annotations"].append(annotation)
    return bam_stats

chrom_color = {
    "human": "red",
    "bonobo": "blue",
    "mixed": "purple"
}

def extract_length(s):
    match = re.search(r'([\-\d,\.]+)\s*\(', s)
    if match:
        num_str = match.group(1).replace(',', '')
        try:
            return float(num_str)
        except ValueError:
            return 0
    return 0

def process_dot_file(dot_path, output_path, bam_stats, DBG_ratio):
    with open(dot_path) as f:
        lines = f.readlines()

    header, nodes_dict, edges, footer = [], {}, [], []
    in_graph = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("digraph"):
            header.append(line)
            in_graph = True
        elif in_graph and stripped == "}":
            footer.append(line)
            in_graph = False
        elif in_graph:
            if '->' in line:
                edges.append(line)
            elif '[' in line:
                node_id = line.strip().split()[0]
                if "\"" in node_id:
                    nodes_dict[node_id] = line
                else:
                    nodes_dict["\"" + node_id + "\""] = line
            else:
                header.append(line)
        else:
            header.append(line)

    # Sort edges by decreasing length
    edges.sort(key=lambda l: extract_length(l), reverse=True)

    # Replace annotations and color for edges using BAM stats
    new_edges = []
    used_nodes = set()
    for line in edges:
        if 'label="' in line:
            label_start = line.find('label="') + len('label="')
            label_end = line.find(' ', label_start)
            if label_end == -1:
                new_edges.append(line)
                continue
            edge_id = line[label_start:label_end]
            if edge_id in bam_stats:
                chrom = bam_stats[edge_id]["chrom"]
                annotations = bam_stats[edge_id]["annotations"]
                color = chrom_color.get(chrom.lstrip('-'), "#000000")
                annotation_str = "".join(["\\n" + a for a in annotations])
                label_quote = line.find('"', label_start)
                label_prefix = line[:label_start]
                label_suffix = line[label_quote:]
                new_line = f'{label_prefix}{line[label_start: line.find(")")+1]}{annotation_str}'f'{label_suffix}'
                color_start = new_line.find('color="')
                color_end = new_line.find('"', color_start + 7)
                if color_start != -1 and color_end != -1:
                    color_prefix = new_line[:color_start]
                    color_suffix = new_line[color_end+1:]
                    new_line = f'{color_prefix}color="{color}"{color_suffix}'
                new_edges.append(new_line)
            else:
                color = "#000000"
                label_quote = line.find('"', label_start)
                label_prefix = line[:label_start]
                label_suffix = line[label_quote:]
                new_line = f'{label_prefix}{line[label_start: line.find(")")+1]}'f'{label_suffix}'
                color_start = new_line.find('color="')
                color_end = new_line.find('"', color_start + 7)
                if color_start != -1 and color_end != -1:
                    color_prefix = new_line[:color_start]
                    color_suffix = new_line[color_end+1:]
                    new_line = f'{color_prefix}color="{color}"{color_suffix}'
                new_edges.append(new_line)
        else:
            new_edges.append(line)

        # Track used node IDs from the edge line
        parts = line.split('->')
        if len(parts) == 2:
            src = parts[0].strip()
            dst = parts[1].split()[0].strip()
            used_nodes.add(src)
            used_nodes.add(dst)

    with open(output_path, "w") as f:
        for line in header:
            f.write(line)
        for line in new_edges:
            f.write(line)
        for node_id in used_nodes:
            if node_id in nodes_dict:
                f.write(nodes_dict[node_id])
        for line in footer:
            f.write(line)

def main():
    parser = argparse.ArgumentParser(description="Replace DOT edge annotations and colors using BAM stats, and reorder nodes/edges by length.")
    parser.add_argument('-d', '--dot', required=True, help='Input DOT file')
    parser.add_argument('-s', '--stats', required=True, help='Input BAM stats file')
    parser.add_argument('-o', '--output', required=True, help='Output DOT file')
    parser.add_argument('-f', '--fasta', required=False, help='Input FASTA file for computing DBG ratios (optional)')
    args = parser.parse_args()

    bam_stats = parse_bam_stats(args.stats)
    DBG_ratio = compute_dbg_ratios(args.fasta, args.output + ".jumboDBG") if args.fasta else {}
    process_dot_file(args.dot, args.output, bam_stats, DBG_ratio)

if __name__ == "__main__":
    main()
