import re
import argparse

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
            if edge_id not in bam_stats:
                bam_stats[edge_id] = {"chrom": chrom, "annotations": []}
            bam_stats[edge_id]["annotations"].append(annotation)
    return bam_stats

chrom_color = {
    "1A": "#325527", "1B": "#325527", "2A": "#628DCF", "2B": "#628DCF", "3A": "#41496B", "3B": "#41496B",
    "4A": "#12CCD6", "4B": "#12CCD6", "5A": "#3E16F3", "5B": "#3E16F3", "6A": "#E46C0A", "6B": "#E46C0A",
    "7A": "#446768", "7B": "#446768", "8A": "#FF0000", "8B": "#FF0000", "9A": "#3C06A6", "9B": "#3C06A6",
    "10A": "#6CB9AB", "10B": "#6CB9AB", "11A": "#988430", "11B": "#988430", "12A": "#4BAA54", "12B": "#4BAA54",
    "13A": "#154E54", "13B": "#154E54", "14A": "#A74C5D", "14B": "#A74C5D", "15A": "#528444", "15B": "#528444",
    "16A": "#B61664", "16B": "#B61664", "17A": "#8F3296", "17B": "#8F3296", "18A": "#E1A9E7", "18B": "#E1A9E7",
    "19A": "#54340D", "19B": "#54340D", "20A": "#316260", "20B": "#316260", "21A": "#8041AF", "21B": "#8041AF",
    "22A": "#5AB499", "22B": "#5AB499", "23A": "#952395", "23B": "#952395", "24A": "#70229F", "24B": "#70229F",
    "25A": "#4D4050", "25B": "#4D4050", "26A": "#969696", "26B": "#969696", "X": "#969696", "Y": "#969696",
    "1M": "#325527", "1P": "#325527", "2M": "#628DCF", "2P": "#628DCF", "3M": "#41496B", "3P": "#41496B",
    "4M": "#12CCD6", "4P": "#12CCD6", "5M": "#3E16F3", "5P": "#3E16F3", "6M": "#E46C0A", "6P": "#E46C0A",
    "7M": "#446768", "7P": "#446768", "8M": "#FF0000", "8P": "#FF0000", "9M": "#3C06A6", "9P": "#3C06A6",
    "10M": "#6CB9AB", "10P": "#6CB9AB", "11M": "#988430", "11P": "#988430", "12M": "#4BAA54", "12P": "#4BAA54",
    "13M": "#154E54", "13P": "#154E54", "14M": "#A74C5D", "14P": "#A74C5D", "15M": "#528444", "15P": "#528444",
    "16M": "#B61664", "16P": "#B61664", "17M": "#8F3296", "17P": "#8F3296", "18M": "#E1A9E7", "18P": "#E1A9E7",
    "19M": "#54340D", "19P": "#54340D", "20M": "#316260", "20P": "#316260", "21M": "#8041AF", "21P": "#8041AF",
    "22M": "#5AB499", "22P": "#5AB499", "23M": "#952395", "23P": "#952395", "Chr1": "#325527", "Chr2": "#628DCF",
    "Chr3": "#41496B", "Chr4": "#12CCD6", "Chr5": "#3E16F3", "Chr6": "#E46C0A", "Chr7": "#446768",
    "Chr8": "#FF0000", "Chr9": "#3C06A6", "Chr10": "#6CB9AB", "Chr11": "#988430", "Chr12": "#4BAA54",
    "Chr13": "#154E54", "Chr14": "#A74C5D"
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

def process_dot_file(dot_path, output_path, bam_stats):
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
                new_line = f'{label_prefix}{line[label_start: line.find(")")+1]}{annotation_str}{label_suffix}'
                color_start = new_line.find('color="')
                color_end = new_line.find('"', color_start + 7)
                if color_start != -1 and color_end != -1:
                    color_prefix = new_line[:color_start]
                    color_suffix = new_line[color_end+1:]
                    new_line = f'{color_prefix}color="{color}"{color_suffix}'
                new_edges.append(new_line)
            else:
                new_edges.append(line)
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
    args = parser.parse_args()

    bam_stats = parse_bam_stats(args.stats)
    process_dot_file(args.dot, args.output, bam_stats)

if __name__ == "__main__":
    main()
