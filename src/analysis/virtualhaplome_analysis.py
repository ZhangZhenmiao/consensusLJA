#!/usr/bin/env python
import argparse
import re
import sys

def merge_intervals(intervals):
    """Calculates the union length of overlapping intervals."""
    if not intervals: return 0
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current
        if curr_start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, curr_end))
        else:
            merged.append(current)
    return sum(end - start for start, end in merged)

def get_rev_node(node_id, full_node_map):
    """Swaps side A and B for a node ID based on the A_B format."""
    parts = full_node_map.get(node_id, "").split('_')
    if len(parts) < 2: return node_id
    return parts[1] if node_id == parts[0] else parts[0]

def filter_contained_paths(paths):
    """Removes paths that are fully contained within another path's coordinates."""
    if not paths: return []
    # Sort by start (ascending) then by length (descending) to keep the largest wrapper first
    paths.sort(key=lambda x: (x['start'], -(x['end'] - x['start'])))
    
    keep = []
    for i, current in enumerate(paths):
        is_contained = False
        for j, other in enumerate(paths):
            if i == j: continue
            # If current is inside other, mark for removal
            if current['start'] >= other['start'] and current['end'] <= other['end']:
                is_contained = True
                break
        if not is_contained:
            keep.append(current)
    return keep

def process_gaf(input_gaf, identity_threshold=0.9):
    contig_segments = {} 
    contig_lengths = {}
    contig_intervals = {} 
    node_to_full = {}

    with open(input_gaf, 'r') as g:
        for line in g:
            items = line.strip().split('\t')
            if len(items) < 12: continue
            
            contig = items[0]
            q_len, q_start, q_end = int(items[1]), int(items[2]), int(items[3])
            path_str = items[5]
            
            try:
                idt_tag = next(x for x in items if x.startswith("id:f:"))
                idt = float(idt_tag[5:])
            except (StopIteration, ValueError):
                idt = 0.0

            if idt >= identity_threshold:
                contig_lengths[contig] = q_len
                contig_intervals.setdefault(contig, []).append((q_start, q_end))
                
                raw_segments = re.findall(r'([><])([^><]+)', path_str)
                processed_nodes = []
                for orientation, node_id in raw_segments:
                    parts = node_id.split('_')
                    side = parts[0] if orientation == '>' else (parts[1] if len(parts) > 1 else parts[0])
                    processed_nodes.append(side)
                    node_to_full[side] = node_id
                
                contig_segments.setdefault(contig, []).append({
                    'nodes': processed_nodes,
                    'start': q_start,
                    'end': q_end
                })

    def glue_paths(segments):
        changed = True
        while changed:
            changed = False
            i = 0
            while i < len(segments):
                p1 = segments[i]
                merged_flag = False
                for j in range(len(segments)):
                    if i == j: continue
                    p2 = segments[j]
                    p2_rev = [get_rev_node(n, node_to_full) for n in reversed(p2['nodes'])]
                    
                    new_nodes = None
                    if p1['nodes'][-1] == p2['nodes'][0]:
                        new_nodes = p1['nodes'] + p2['nodes'][1:]
                    # elif p1['nodes'][-1] == p2_rev[0]:
                    #     new_nodes = p1['nodes'] + p2_rev[1:]
                    elif p2['nodes'][-1] == p1['nodes'][0]:
                        new_nodes = p2['nodes'] + p1['nodes'][1:]
                    # elif p2_rev[-1] == p1['nodes'][0]:
                    #     new_nodes = p2_rev + p1['nodes'][1:]

                    if new_nodes:
                        segments[i] = {
                            'nodes': new_nodes,
                            'start': min(p1['start'], p2['start']),
                            'end': max(p1['end'], p2['end'])
                        }
                        segments.pop(j)
                        changed = merged_flag = True
                        break
                if not merged_flag: i += 1
        return segments

    file_rows = []
    print(f"{'Contig':<25} | {'Overall Cov':<12} | {'Paths Kept'}")
    print("-" * 52)

    for contig in sorted(contig_lengths.keys()):
        c_len = contig_lengths[contig]
        overall_cov_bases = merge_intervals(contig_intervals.get(contig, []))
        overall_cov_ratio = overall_cov_bases / c_len
        
        # 1. Glue based on graph nodes
        glued = glue_paths(contig_segments.get(contig, []))
        
        # 2. Filter out contained paths based on coordinates
        filtered = filter_contained_paths(glued)
        
        # 3. Final sort by Path_Start
        filtered.sort(key=lambda x: x['start'])
        
        print(f"{contig:<25} | {overall_cov_ratio:<12.4%} | {len(filtered)}")
        
        for p in filtered:
            p_start, p_end = p['start'], p['end']
            p_cov = (p_end - p_start) / c_len
            file_rows.append([
                contig, c_len, p_start, p_end, f"{p_cov:.4f}", ",".join(p['nodes'])
            ])
            
    return file_rows

def main():
    parser = argparse.ArgumentParser(description="Process GAF paths, glue edges, and remove contained paths.")
    parser.add_argument("input", help="Input GAF file")
    parser.add_argument("-o", "--output", help="Output TSV file", required=True)
    parser.add_argument("-t", "--threshold", type=float, default=0.9, help="Identity threshold")
    args = parser.parse_args()

    results = process_gaf(args.input, args.threshold)
    
    headers = ["Contig", "Contig_Len", "Path_Start_on_Contig", "Path_End_on_Contig", "Path_Cov_on_Contig", "Edges"]
    
    with open(args.output, 'w') as f:
        f.write("\t".join(headers) + "\n")
        for row in results:
            f.write("\t".join(map(str, row)) + "\n")
    
    print(f"\nSaved {len(results)} non-contained paths to {args.output}")

if __name__ == "__main__":
    main()