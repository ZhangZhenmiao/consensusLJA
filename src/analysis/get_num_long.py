import re
import sys

def count_edges(dot_file, length_threshold=20000):
    # Regex to extract edge length: number before '('
    pattern = re.compile(r'label="[^"]*?(\d+)\([^)]*\)"')

    total_edges = 0
    long_edges = 0

    with open(dot_file, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                total_edges += 1
                length = int(match.group(1))
                if length > length_threshold:
                    long_edges += 1

    return total_edges, long_edges

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python {sys.argv[0]} <graph.dot>")
        sys.exit(1)

    dot_file = sys.argv[1]
    total, long = count_edges(dot_file)
    print(f"Total edges: {total}")
    print(f"Long edges (>20k): {long}")
