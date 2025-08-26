#!/usr/bin/env python
import argparse
import re

from collections import Counter

def Most_Common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]

def get_paths(graphaligner):
    contig2paths = {}
    with open(graphaligner) as g:
            line = g.readline()
            while(line):
                items = line.strip().split()
                contig = items[0]
                genome_path = items[5]
                genome_path = re.split("<|>", genome_path[1:])
                idt = float(items[-2][5:])
                if idt >= 0.9:
                    if contig in contig2paths:
                        contig2paths[contig].append(genome_path)
                    else:
                        contig2paths[contig] = [genome_path]
                line = g.readline()

    def get_edge_ids(path):
        # Get the first ID of the first edge and the second ID of the last edge
        start_edge = path[0]
        end_edge = path[-1]
        
        start_edge = start_edge[:start_edge.find('.')] + '_' + start_edge[start_edge.find('_') + 1 : start_edge.rfind('.')]
        end_edge = end_edge[:end_edge.find('.')] + '_' + end_edge[end_edge.find('_') + 1 : end_edge.rfind('.')]
        return start_edge, end_edge

    def glue_paths(paths):
        """
        Function to glue as many paths together as possible if the end ID of one path matches
        the start ID of another path. It handles unordered paths by attempting to glue until no more
        changes occur.
        """
        changed = True
        while changed:  # Keep iterating until no paths can be glued further
            changed = False
            i = 0
            while i < len(paths):
                current_path = paths[i]
                current_start, current_end = get_edge_ids(current_path)
                merged = False
                
                for j in range(len(paths)):
                    if i != j:  # Avoid self-comparison
                        next_path = paths[j]
                        next_start, next_end = get_edge_ids(next_path)

                        # Check if current path's end matches next path's start
                        if current_end == next_start or current_end.split('_')[0] in next_start or current_end.split('_')[1] in next_start:
                            paths[i] = [*current_path,*next_path[1:]]  # Glue the paths
                            paths.pop(j)  # Remove the glued path
                            changed = True
                            merged = True
                            break
                        
                        # Check if current path's start matches next path's end (reverse order)
                        elif current_start == next_end or current_start.split('_')[0] in next_end or current_start.split('_')[1] in next_end:
                            paths[i] = [*next_path,*current_path[1:]]  # Glue paths in reverse order
                            paths.pop(j)  # Remove the glued path
                            changed = True
                            merged = True
                            break

                if not merged:
                    i += 1  # Move to the next path if no merging occurred
        
        return paths

    for c in contig2paths:
        contig2paths[c] = glue_paths(contig2paths[c])
    
    return contig2paths

def get_switches(graphaligner, paths, dot, out):
    edge2len = {}
    edge2cov = {}
    with open(dot) as d:
        line = d.readline()
        while line:
            if "->" in line:
                label = line[line.find("label="):]
                label = label[7:label.find("\" ")]
                la, _, le = label.split()
                edge2len[la] = le[:le.find('(')]
                edge2cov[la] = le[le.find('(')+1:le.find(')')]
            line = d.readline()

    edge2ref = {}
    with open(graphaligner) as g:
        line = g.readline()
        sum_idt = 0
        smallest_idt = 100
        cnt = 0
        cnt_pass = 0
        cnt_edges = 0
        while(line):
            items = line.strip().split()
            chr = items[0]
            if "tig" not in chr:
                chr = chr[chr.find('_') + 1:]
                genome_path = items[5]
                genome_path = re.split("<|>", genome_path[1:])
                idt = float(items[-2][5:])
                sum_idt += idt
                cnt += 1
                if idt < smallest_idt:
                    smallest_idt = idt
                if idt >= 0.99:
                    cnt_pass += 1
                    for edge in genome_path:
                        cnt_edges += 1
                        if edge not in edge2ref:
                            edge2ref[edge] = set()
                        edge2ref[edge].add(chr)
                    # print(chr, cnt_edges)
            line = g.readline()
    print(f"Ref stats: average idt {sum_idt/cnt}, smallest idt {smallest_idt}, paths with valid idt {cnt_pass}, all paths {cnt}, edges in ref {cnt_edges}")
    for e in edge2ref:
        edge2ref[e] = list(edge2ref[e])
        # print(e, edge2ref[e])

    is_switch = 0
    not_switch = 0
    is_special = 0
    contig2paths = get_paths(paths)
    with open(out, 'w') as w:
        w.write('\t'.join(["Path_ID_in_GA","Contig_name", "Edge_in_DBG", "Edge_len", "Edge_cov", "Ref_list", "Ref_of_the_path", "Status"]) + '\n')
        pathid = 0
        for c in contig2paths:
            for genome_path in contig2paths[c]:
                if len(genome_path) < 2:
                    continue
                refs = []
                for edge in genome_path:
                    if edge in edge2ref:
                        for x in edge2ref[edge]:
                            refs.append(x[:-1])
                if len(refs) > 0:
                    pathid += 1
                    majority = Most_Common(refs)
                    if majority == "17": continue
                    previous = ""
                    for e in genome_path:
                        if e in edge2ref:
                            ecov = edge2cov[e.split('_')[0]]
                            elen = edge2len[e.split('_')[0]]
                            w.write('\t'.join([str(pathid),c, e, elen, ecov, str(edge2ref[e]), majority]) + '\t')
                            if previous == "":
                                previous = edge2ref[e]
                                w.write("NOT SWITCH\n")
                                not_switch += 1
                            elif len(list(set(previous) & set(edge2ref[e]))) != 0:
                                previous = list(set(previous) & set(edge2ref[e]))
                                w.write("NOT SWITCH\n")
                                not_switch += 1
                            else:
                                prev_chr = [x[:-1] for x in previous]
                                curr_chr = [x[:-1] for x in edge2ref[e]]
                                if len(list(set(prev_chr) & set(curr_chr))) != 0:
                                    w.write("IS SWITCH\n")
                                    is_switch += 1
                                else:
                                    w.write("IS SWITCH (SPECIAL)\n")
                                    is_special += 1
                                previous = edge2ref[e]
                        else:
                            ecov = edge2cov[e.split('_')[0]]
                            elen = edge2len[e.split('_')[0]]
                            w.write('\t'.join([str(pathid), c, e, elen, ecov,  "NO LABEL", majority, "NOT SWITCH"]) + '\n')
    
    print(f"Evaluated switches\n\
            IS SWITCH {is_switch}\n\
            NOT SWITCH {not_switch}\n\
            SPECIAL SWITCH {is_special}\n\
            SWITCH RATE {0 if is_switch + is_special + not_switch ==0 else (is_switch + is_special)/(is_switch + is_special + not_switch)} ({is_switch + is_special}/{is_switch + is_special + not_switch})\
           ")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Phase switch of consensus assembly.")
    parser.add_argument("ref_gaf", help="Path to ref paths by graphaligner")
    parser.add_argument("asm_gaf", help="Path to the asm paths by graphaligner")
    parser.add_argument("ref_dot", help="Path to ref jumbodbg")
    parser.add_argument("out", help="Path to output")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function with the BAM file path provided by the user
    # read_bam(args.bam_file, args.reads_file)
    get_switches(args.ref_gaf, args.asm_gaf, args.ref_dot, args.out)
