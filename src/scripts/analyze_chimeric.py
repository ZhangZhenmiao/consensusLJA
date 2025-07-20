#!/usr/bin/env python
from Bio import SeqIO
import pysam
import argparse
import numpy as np
import sys
import os

def calculate_identity(alignment):
    cigar = alignment.cigartuples
    if not cigar:
        return 0

    matches = sum(length for op, length in cigar if op == 7)  # '='
    mismatches = sum(length for op, length in cigar if op == 8)  # 'X'
    insertions = sum(length for op, length in cigar if op == 1)
    deletions = sum(length for op, length in cigar if op == 2)

    aligned_len = matches + mismatches # query_alignment_length is aligned_len + insertions
    total_bases = aligned_len + insertions + deletions
    pid = matches / total_bases if total_bases > 0 else 0.0
    return pid

# def read_bam(bam_file_path, reads_file_path):
def read_bam(bam_file_path, dot_file_path):
    # if (reads != ""):
    #     print(reads, flush=True)
    #     if reads == "-":
    #         handle = sys.stdin
    #     else:
    #         handle = open(reads, "r")

    #     read_lengths = {}
    #     for record in SeqIO.parse(handle, "fasta"):
    #         read_lengths[record.id] = len(record.seq)

    #     print("read", len(read_lengths), "read lengths", flush=True)

    node2len = {}
    with open(dot_file_path) as dot_file:
        line = dot_file.readline()
        while(line):
            if "label" not in line:
                line = dot_file.readline()
                continue
            if "->" in line:
                line = dot_file.readline()
                continue
            node = line.split()[0]
            if "+" not in node:
                node_len = line[line.find('L')+1: line.find('"', line.find('L')+1)]
                node2len[node] = int(node_len)
                # print("Length of " + node, node_len, flush=True)
            line = dot_file.readline()

    # Open the BAM file
    alignments = {}
    chimeric_edges = []
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        # store alignments for each edge in the dict alignments
        for read in bam_file:
            if read.is_unmapped:
                continue
            # if read.reference_name != "84021.66_-85034.68":
            #     continue
            if calculate_identity(read) > 0.99:
                if read.reference_name not in alignments:
                    alignments[read.reference_name] = []
                alignments[read.reference_name].append({
                    "name": read.query_name,
                    "start": read.reference_start,
                    "end": read.reference_end,
                    "ref_length": bam_file.get_reference_length(read.reference_name),
                    "aligned": read.query_alignment_length,
                    "idt": calculate_identity(read),
                    # "aligned_fraction": read.query_alignment_length/read_lengths[read.query_name]
                })

        # detect chimeric at nodes
        tolerant_size = 1000
        # print("Contig", "Contig_len", "Node_start", "Node_end", "Read", "Aln_start", "Aln_end", "Aligned_len", "PI", "AF", sep="\t", flush=True)
        for r in alignments:
            sorted(alignments[r], key=lambda x: x["start"])

            node1, node2 = r.split('_')
            node1 = node1[:node1.find('.')]
            node2 = node2[:node2.find('.')]
            contig_len = bam_file.get_reference_length(r)

            # if r == "84021.66_-85034.68":
            #     for aln in alignments[r]:
            #         print(r, contig_len, node2len[node1], node2len[node2], aln["name"], aln["start"], aln["end"], aln["aligned"], aln["idt"], sep="\t", flush=True)

            split_coordinate1 = min(node2len[node1], contig_len - node2len[node2])
            split_coordinate2 = max(node2len[node1], contig_len - node2len[node2])
            if split_coordinate1 <= 2*tolerant_size or split_coordinate2 >= contig_len-2*tolerant_size:
                continue
            if node2len[node1] + node2len[node2] - contig_len >= 10000 - tolerant_size:
                continue

            cnt_f = 0
            cnt_r = 0
            cnt_all = 0
            for aln in alignments[r]:
                if aln["start"] <= max(split_coordinate1 - tolerant_size, 0) and aln["end"] > node2len[node1]:
                    cnt_f += 1
                if aln["start"] < contig_len - node2len[node2] and aln["end"] >= min(split_coordinate2 + tolerant_size, contig_len):
                    cnt_r += 1
                cnt_all += 1
            # print("contig name", r, "contig len", contig_len, "len 1", node2len[node1], "len 2", node2len[node2], "supporting reads", cnt_f, cnt_r, flush=True)
            if cnt_f <= 0 or cnt_r <= 0:
                print("contig name", r, "contig len", contig_len, "len 1", node2len[node1], "len 2", node2len[node2], "supporting reads", cnt_f, cnt_r, cnt_all, flush=True)
                print(r, "is chimeric", flush=True)
                chimeric_edges.extend(r.split('_'))

        # detect internal chimeric
        tolerant_size = 100
        for r in alignments:
            coverages = np.zeros(bam_file.get_reference_length(r), dtype=int)
            for i in alignments[r]:
                if i["end"]-tolerant_size > i["start"]:
                    coverages[i["start"]:i["end"]-tolerant_size] += 1
            
            potential_pos = []
            if np.average(coverages) >= 5:
                for i in range(bam_file.get_reference_length(r)):
                    if coverages[i] <= 1 and i >= 5001- tolerant_size and i <= bam_file.get_reference_length(r)-5001+tolerant_size:
                        potential_pos.append(i)
                        # print("low coverage:", r, i)
            
            reads = set()
            for i in alignments[r]:
                for j in potential_pos:
                    if i["start"] <= j-5001+tolerant_size and i["end"] >= j+5001-tolerant_size:
                        reads.add(i["name"])
            
            # print IDs of chimeric reads
            if len(reads) != 0:
                print(r, "is chimeric (internal)", flush=True)
                chimeric_edges.extend(r.split('_'))
    return chimeric_edges
            
if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Read a BAM file and print alignment details.")
    parser.add_argument("bam_file", help="Path to the input BAM file")
    parser.add_argument("dot_file", help="Path to the input dot file (for node sizes)")
    parser.add_argument("output", help="Path to the output file")
    # parser.add_argument("reads", help="Path to reads", default="", required=False)
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function with the BAM file path provided by the user
    # read_bam(args.bam_file, args.reads_file)
    # chimeric_edges = read_bam(args.bam_file, args.dot_file, args.reads)
    if not os.path.isfile(args.output):
        chimeric_edges = read_bam(args.bam_file, args.dot_file)
        with open(args.output, "w") as fout:
            for e in chimeric_edges:
                fout.write(e + '\n')

