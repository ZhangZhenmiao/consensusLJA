import re
import sys
from collections import defaultdict

def classify_dot_edges(dot_file):
    with open(dot_file, 'r') as f:
        dot_text = f.read()
    # Pattern to extract edge info
    edge_pattern = re.compile(r'"?([^"]+)"?\s*->\s*"?(.*?)"?\s*\[label="([^"]+)"')

    edges = []
    node_counts = defaultdict(int)

    # First pass: parse all edges and count node usage
    for line in dot_text.strip().splitlines():
        match = edge_pattern.search(line)
        if match:
            src, tgt, label = match.groups()
            label = label[:label.find(" ")]
            edges.append((src, tgt, label))
            node_counts[src] += 1
            if tgt != src:
                node_counts[tgt] += 1

    # Second pass: classify edges
    results = defaultdict(str)
    for src, tgt, label in edges:
        if node_counts[src] > 1 or node_counts[tgt] > 1:
            classification = "multi-edge"
        else:
            classification = "isolated"
        results[label] = classification

    return results

def main():
    if len(sys.argv) != 2:
        print("Usage: python classify_dot_edges.py input.dot")
        sys.exit(1)

    filename = sys.argv[1]

    

    results = classify_dot_edges(filename)

    for label, kind in results.items():
        print(f"{label}: {kind}")

if __name__ == "__main__":
    main()
