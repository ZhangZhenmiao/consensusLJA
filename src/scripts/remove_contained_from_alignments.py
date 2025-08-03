#!/usr/bin/env python
import pysam
import argparse
import sys
import re
from collections import defaultdict
import copy

def true_query_start_end(cigar):
    # Parse the CIGAR string into list of (length, op)
    tokens = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    tokens = [(int(length), op) for length, op in tokens]

    # Compute the total length of the read (excluding hard clips)
    read_length = sum(length for length, op in tokens if op in ('M', 'I', 'S', '=', 'X', 'H'))

    # Count left and right clip (soft + hard)
    left_clip = 0
    right_clip = 0

    if tokens:
        if tokens[0][1] in ('S', 'H'):
            left_clip = tokens[0][0]
        if tokens[-1][1] in ('S', 'H'):
            right_clip = tokens[-1][0]

    # true query start is where alignment begins in full read
    true_query_start = left_clip
    # true query end is where alignment ends in full read
    true_query_end = read_length - right_clip

    return true_query_start, true_query_end, read_length

def calculate_identities(alignment, gap_threshold=10):
    cigar = alignment.cigartuples
    if not cigar:
        return 0

    matches = sum(length for op, length in cigar if op == 7)  # '='
    mismatches = sum(length for op, length in cigar if op == 8)  # 'X'
    insertions = sum(length for op, length in cigar if op == 1)
    long_gaps1 = sum(length for op, length in cigar if op == 1 and length >= gap_threshold)
    deletions = sum(length for op, length in cigar if op == 2)
    long_gaps2 = sum(length for op, length in cigar if op == 2 and length >= gap_threshold)

    aligned_len = matches + mismatches # query_alignment_length is aligned_len + insertions
    total_bases = aligned_len + insertions + deletions
    pid = matches / total_bases if total_bases > 0 else 0.0

    pid_nogap = matches / (aligned_len + insertions + deletions - long_gaps1 - long_gaps2)
    
    return pid, pid_nogap

def filter_alignments_with_identity(bam_file_path, threshold=0):
    high_identity_alignments = {}
    total_alignments = 0
    processed_alignments = 0

    query_to_ref = {}

    with pysam.AlignmentFile(bam_file_path, "r") as bamfile:
        for alignment in bamfile:
            total_alignments += 1

            if alignment.is_unmapped:
                continue

            identity, identity_nogap = calculate_identities(alignment, gap_threshold=10)
            
            if identity_nogap >= threshold:
                processed_alignments += 1
                
                query_header = alignment.query_name.split('_')
                linear_edge_name = query_header[0] + "_" + query_header[1]
                query_name = ""
                if query_header[2] == "0":
                    query_name = linear_edge_name
                    query_to_ref[query_name] = (linear_edge_name, 0)
                elif query_header[2] == "1":
                    query_name = query_header[0]
                    query_to_ref[query_name] = (linear_edge_name, 1)
                else:
                    query_name = query_header[1]
                    query_to_ref[query_name] = (linear_edge_name, 2)

                ref_id = alignment.reference_name

                # skip records aligned to the same linear edge
                if ref_id == linear_edge_name:
                    continue

                # print(ref_id, linear_edge_name)
                
                ref_len = bamfile.get_reference_length(alignment.reference_name)

                query_alignment_start, query_alignment_end, query_len = true_query_start_end(alignment.cigarstring)

                aln = {
                    'identity': identity,
                    'identity_nogap': identity_nogap,
                    'ref_id': ref_id,
                    'ref_start': alignment.reference_start,
                    'ref_end': alignment.reference_end,
                    'query_start': query_alignment_start,
                    'query_end': query_alignment_end,
                    'reverse': alignment.is_reverse,
                    'aln_length_query': query_alignment_end - query_alignment_start,
                    'aln_length_ref': alignment.reference_end - alignment.reference_start,
                    'length_query': query_len,
                    'length_ref': ref_len,
                    'length_entire_query': bamfile.get_reference_length(linear_edge_name),
                    'alignment': alignment
                }

                # if linear_edge_name == "-5274_4784":
                #     print(query_name, ref_id, aln["aln_length_query"], aln["query_start"], aln["query_end"], aln["ref_start"], aln["ref_end"], aln["reverse"])

                if ref_id not in high_identity_alignments:
                    high_identity_alignments[ref_id] = {}
                
                if linear_edge_name not in high_identity_alignments[ref_id]:
                    high_identity_alignments[ref_id][linear_edge_name] = {
                        "start": [],
                        "end": [],
                        "entire": []
                    }

                if query_header[2] == "0":
                    high_identity_alignments[ref_id][linear_edge_name]["entire"].append(aln)
                elif query_header[2] == "1":
                    high_identity_alignments[ref_id][linear_edge_name]["start"].append(aln)
                else:
                    high_identity_alignments[ref_id][linear_edge_name]["end"].append(aln)
    return high_identity_alignments

def merge_and_filter_alignments_ref(alns, overlap_thresh=0.9):
    if not alns:
        return [], []

    # Work on copies to preserve original input
    forward = [copy.deepcopy(a) for a in alns if not a["reverse"]]
    reverse = [copy.deepcopy(a) for a in alns if a["reverse"]]

    def process_group(group):
        # Sort by ref_start
        group.sort(key=lambda x: (x["ref_start"], -x["ref_end"]))
        merged = []

        for aln in group:
            added = False
            for m in merged:
                # Compute overlap on reference
                start1, end1 = m["ref_start"], m["ref_end"]
                start2, end2 = aln["ref_start"], aln["ref_end"]

                overlap = max(0, min(end1, end2) - max(start1, start2))
                span1 = end1 - start1
                span2 = end2 - start2
                min_span = min(span1, span2)

                if min_span > 0 and overlap / min_span >= overlap_thresh:
                    # Keep the one with larger span
                    if span2 > span1:
                        m.update(copy.deepcopy(aln))
                    added = True
                    break

                # Check if within mergeable distance
                max_gap = max(100000, min(200000, m["length_query"]/5))
                gap = start2 - end1
                if gap <= max_gap:
                    # Merge: extend ref and query coordinates
                    m["ref_end"] = max(m["ref_end"], aln["ref_end"])
                    m["query_end"] = max(m["query_end"], aln["query_end"])
                    m["ref_start"] = min(m["ref_start"], aln["ref_start"])
                    m["query_start"] = min(m["query_start"], aln["query_start"])
                    added = True
                    break

            if not added:
                merged.append(copy.deepcopy(aln))

        if not merged:
            return []

        # Select the record with the largest ref span
        best = max(merged, key=lambda x: x["ref_end"] - x["ref_start"])
        ref_span = best["ref_end"] - best["ref_start"]
        if ref_span < 0.8 * best["length_query"] or best["length_query"] < 0.8 * ref_span:
            return []  # reject if too short

        return [best]

    return process_group(forward), process_group(reverse)

def merge_and_filter_alignments_query(alns, overlap_thresh=0.9):
    if not alns:
        return [], []

    # Work on deep copies to preserve original input
    forward = [copy.deepcopy(a) for a in alns if not a["reverse"]]
    reverse = [copy.deepcopy(a) for a in alns if a["reverse"]]

    def process_group(group):
        # Sort by query_start
        group.sort(key=lambda x: (x["query_start"], -x["query_end"]))
        merged = []

        for aln in group:
            added = False
            for m in merged:
                # Compute overlap on query
                start1, end1 = m["query_start"], m["query_end"]
                start2, end2 = aln["query_start"], aln["query_end"]

                overlap = max(0, min(end1, end2) - max(start1, start2))
                span1 = end1 - start1
                span2 = end2 - start2
                min_span = min(span1, span2)

                if min_span > 0 and overlap / min_span >= overlap_thresh:
                    # Keep the one with larger span
                    if span2 > span1:
                        m.update(copy.deepcopy(aln))
                    added = True
                    break

                # Check if within mergeable distance
                max_gap = max(100000, min(200000, m["length_query"]/5))
                gap = start2 - end1
                if gap <= max_gap:
                    # Merge: extend ref and query coordinates
                    m["ref_end"] = max(m["ref_end"], aln["ref_end"])
                    m["query_end"] = max(m["query_end"], aln["query_end"])
                    m["ref_start"] = min(m["ref_start"], aln["ref_start"])
                    m["query_start"] = min(m["query_start"], aln["query_start"])
                    added = True
                    break

            if not added:
                merged.append(copy.deepcopy(aln))

        if not merged:
            return []

        # Select the record with the largest ref span
        best = max(merged, key=lambda x: x["query_end"] - x["query_start"])
        ref_span = best["ref_end"] - best["ref_start"]
        query_span = best["query_end"] - best["query_start"]
        if ref_span < 0.8 * query_span or query_span < 0.8 * ref_span:
            return []

        return [best]

    return process_group(forward), process_group(reverse)

def find_contained_contigs(high_identity_alignments):
    contained_contigs = set()
    # for ref in high_identity_alignments:
    #     for edge in high_identity_alignments[ref]:
    #         data = high_identity_alignments[ref][edge]
    #         entire_forward, entire_reverse = merge_and_filter_alignments_ref(data["entire"])
    #         start_forward, start_reverse = merge_and_filter_alignments_ref(data["start"])
    #         end_forward, end_reverse = merge_and_filter_alignments_ref(data["end"])
            
    #         flag = False
    #         # process entire alignment
    #         if len(entire_forward) == 1:
    #             aln = entire_forward[0]
    #             span_ref = aln["ref_end"] - aln["ref_start"]
    #             ratio = min(aln["length_entire_query"], span_ref) / max(aln["length_entire_query"], span_ref)
    #             if ratio >= 0.8:
    #                 if aln["length_entire_query"] < aln["length_ref"]:
    #                     print(f'Edge {edge} is contained in {ref}, identity {aln["identity"]}, lengths {aln["length_entire_query"]} and {aln["length_ref"]}')
    #                     contained_contigs.add(edge)
    #                 else:
    #                     print(f'Edge {ref} is contained in {edge}, identity {aln["identity"]}, lengths {aln["length_ref"]} and {aln["length_entire_query"]}')
    #                     contained_contigs.add(ref)
    #                 flag = True
            
    #         if flag: continue
            
    #         if len(entire_reverse) == 1:
    #             aln = entire_reverse[0]
    #             span_ref = max(aln["ref_end"], aln["ref_end"]) - min(aln["ref_start"], aln["ref_start"])
    #             ratio = min(aln["length_entire_query"], span_ref) / max(aln["length_entire_query"], span_ref)
    #             if ratio >= 0.8:
    #                 if aln["length_entire_query"] <= aln["length_ref"]:
    #                     print(f'Edge {edge} is contained in {ref}, identity {aln["identity"]}, lengths {aln["length_entire_query"]} and {aln["length_ref"]}')
    #                     contained_contigs.add(edge)
    #                 else:
    #                     print(f'Edge {ref} is contained in {edge}, identity {aln["identity"]}, lengths {aln["length_ref"]} and {aln["length_entire_query"]}')
    #                     contained_contigs.add(ref)
    #                 flag = True
            
    #         if flag: continue
            
    #         if len(start_reverse) == 1 and len(end_reverse) == 1:
    #             aln_start = start_reverse[0]
    #             aln_end = end_reverse[0]

    #             span_ref = max(aln_start["ref_end"], aln_end["ref_end"]) - min(aln_start["ref_start"], aln_end["ref_start"])
    #             ratio = min(aln_start["length_entire_query"], span_ref) / max(aln_start["length_entire_query"], span_ref)
    #             if ratio >= 0.8:
    #                 if aln_start["length_entire_query"] <= aln_start["length_ref"]:
    #                     print(f'Edge {edge} is contained in {ref}, lengths edge {aln_start["length_entire_query"]} and ref span {span_ref}')
    #                     contained_contigs.add(edge)
    #                     flag = True
            
    #         if flag: continue
            
    #         if len(start_forward) == 1 and len(end_forward) == 1:
    #             aln_start = start_forward[0]
    #             aln_end = end_forward[0]

    #             span_ref = max(aln_start["ref_end"], aln_end["ref_end"]) - min(aln_start["ref_start"], aln_end["ref_start"])
    #             ratio = min(aln_start["length_entire_query"], span_ref) / max(aln_start["length_entire_query"], span_ref)
    #             if ratio >= 0.8:
    #                 if aln_start["length_entire_query"] <= aln_start["length_ref"]:
    #                     print(f'Edge {edge} is contained in {ref}, lengths edge {aln_start["length_entire_query"]} and ref span {span_ref}')
    #                     contained_contigs.add(edge)
    
    for ref in high_identity_alignments:
        for edge in high_identity_alignments[ref]:
            data = high_identity_alignments[ref][edge]
            entire_forward, entire_reverse = merge_and_filter_alignments_query(data["entire"])
            start_forward, start_reverse = merge_and_filter_alignments_query(data["start"])
            end_forward, end_reverse = merge_and_filter_alignments_query(data["end"])
            
            flag = False
            # process entire alignment
            if len(entire_forward) == 1:
                aln = entire_forward[0]
                span_ref = aln["ref_end"] - aln["ref_start"]
                ratio = min(aln["length_entire_query"], span_ref) / max(aln["length_entire_query"], span_ref)
                if ratio >= 0.8:
                    if aln["length_entire_query"] < aln["length_ref"]:
                        print(f'Edge {edge} is contained in {ref}, identity {aln["identity"]}, lengths {aln["length_entire_query"]} and {aln["length_ref"]}')
                        contained_contigs.add(edge)
                    else:
                        print(f'Edge {ref} is contained in {edge}, identity {aln["identity"]}, lengths {aln["length_ref"]} and {aln["length_entire_query"]}')
                        contained_contigs.add(ref)
                    flag = True
            
            if flag: continue
            
            if len(entire_reverse) == 1:
                aln = entire_reverse[0]
                span_ref = max(aln["ref_end"], aln["ref_end"]) - min(aln["ref_start"], aln["ref_start"])
                ratio = min(aln["length_entire_query"], span_ref) / max(aln["length_entire_query"], span_ref)
                if ratio >= 0.8:
                    if aln["length_entire_query"] <= aln["length_ref"]:
                        print(f'Edge {edge} is contained in {ref}, identity {aln["identity"]}, lengths {aln["length_entire_query"]} and {aln["length_ref"]}')
                        contained_contigs.add(edge)
                    else:
                        print(f'Edge {ref} is contained in {edge}, identity {aln["identity"]}, lengths {aln["length_ref"]} and {aln["length_entire_query"]}')
                        contained_contigs.add(ref)
                    flag = True
            
            if flag: continue
            
            if len(start_reverse) == 1 and len(end_reverse) == 1:
                aln_start = start_reverse[0]
                aln_end = end_reverse[0]

                span_ref = max(aln_start["ref_end"], aln_end["ref_end"]) - min(aln_start["ref_start"], aln_end["ref_start"])
                ratio = min(aln_start["length_entire_query"], span_ref) / max(aln_start["length_entire_query"], span_ref)
                if ratio >= 0.8:
                    if aln_start["length_entire_query"] <= aln_start["length_ref"]:
                        print(f'Edge {edge} is contained in {ref}, lengths edge {aln_start["length_entire_query"]} and ref span {span_ref}')
                        contained_contigs.add(edge)
                        flag = True
            
            if flag: continue
            
            if len(start_forward) == 1 and len(end_forward) == 1:
                aln_start = start_forward[0]
                aln_end = end_forward[0]

                span_ref = max(aln_start["ref_end"], aln_end["ref_end"]) - min(aln_start["ref_start"], aln_end["ref_start"])
                ratio = min(aln_start["length_entire_query"], span_ref) / max(aln_start["length_entire_query"], span_ref)
                if ratio >= 0.8:
                    if aln_start["length_entire_query"] <= aln_start["length_ref"]:
                        print(f'Edge {edge} is contained in {ref}, lengths edge {aln_start["length_entire_query"]} and ref span {span_ref}')
                        contained_contigs.add(edge)

    return contained_contigs


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Extract high-identity alignments from a BAM file.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("-t", "--threshold", type=float, default=0.95,
                        help="Identity threshold (default: 0.95).")
    parser.add_argument("-o", "--output", help="Path to the output file. If not specified, prints to stdout.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    bam_file_path = args.bam_file
    threshold = args.threshold
    output_path = args.output

    high_identity_alignments = filter_alignments_with_identity(bam_file_path, threshold=threshold)
    contained_edges = find_contained_contigs(high_identity_alignments)
    if output_path:
        with open(output_path, "w") as out_f:
            for edge in contained_edges:
                out_f.write(edge + "\n")
    else:
        for edge in contained_edges:
            print(edge)
    

if __name__ == "__main__":
    main()
