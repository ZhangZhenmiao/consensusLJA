#!/usr/bin/env python
import os
import argparse
import subprocess
from Bio import SeqIO
from collections import defaultdict
import pysam
import re

def calculate_identity(alignment):
    cigar = alignment.cigartuples
    if not cigar:
        return 0

    matches = sum(length for op, length in cigar if op == 7)  # '='
    mismatches = sum(length for op, length in cigar if op == 8)  # 'X'
    insertions = sum(length for op, length in cigar if op == 1)
    long_gaps_i = sum(length for op, length in cigar if op == 1 and length >= 10)
    deletions = sum(length for op, length in cigar if op == 2)
    long_gaps_d = sum(length for op, length in cigar if op == 2 and length >= 10)

    aligned_len = matches + mismatches # query_alignment_length is aligned_len + insertions
    total_bases = aligned_len + insertions + deletions
    pid = matches / total_bases if total_bases > 0 else 0.0
    total_bases = total_bases - long_gaps_i - long_gaps_d  # Exclude long gaps from total bases
    pid_nogap = matches / total_bases if total_bases > 0 else 0.0
    return pid, pid_nogap

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

def extract_flanks(input_fasta, output_fasta, flank_size=20000):
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_len = len(record.seq)
            left = record.seq[:flank_size]
            right = record.seq[-flank_size:] if seq_len >= flank_size else record.seq
            out_f.write(f">{record.id}_head\n{left}\n")
            out_f.write(f">{record.id}_tail\n{right}\n")
    print(f"[Connect] Flanks written to: {output_fasta}")

def run_minimap2_to_bam(flank_fasta, hifi_reads, output_bam, compress, threads):
    cmd = f"{compress} --dimer-compress 32,32,1 --reads {hifi_reads} | minimap2 -ax map-hifi -Y --eqx --sam-hit-only --secondary=no -t {threads} {flank_fasta} - | samtools sort -@ {threads} -o {output_bam}"
    
    os.system(cmd)
    subprocess.run(f"samtools index {output_bam}", shell=True, check=True)
    print(f"[Connect] Sorted BAM + index written: {output_bam}, {output_bam}.bai")

def filter_alignments_with_identity(bam_file_path, threshold=0.9):
    high_identity_alignments = defaultdict(list)
    total_alignments = 0
    processed_alignments = 0
    # spanning_reads_per_ref = defaultdict(set)

    with pysam.AlignmentFile(bam_file_path, "r") as bamfile:
        for alignment in bamfile:
            total_alignments += 1

            if alignment.is_unmapped:
                continue

            identity, identity_nogap = calculate_identity(alignment)
            
            if identity >= threshold:
                processed_alignments += 1
                ref_id = alignment.reference_name
                query_name = alignment.query_name
                
                ref_len = bamfile.get_reference_length(alignment.reference_name)

                query_alignment_start, query_alignment_end, query_len = true_query_start_end(alignment.cigarstring)

                if query_alignment_end - query_alignment_start >= 3000:
                    ref_start = alignment.reference_start if alignment.is_forward else ref_len - alignment.reference_end
                    ref_end = alignment.reference_end if alignment.is_forward else ref_len - alignment.reference_start
                    query_start = query_alignment_start if alignment.is_forward else query_len - query_alignment_end
                    query_end = query_alignment_end if alignment.is_forward else query_len - query_alignment_start

                    flag = False
                    if ref_start <= 20 and query_start > ref_start:
                        flag = True
                    if ref_len - ref_end <= 20 and query_len - query_end > ref_len - ref_end:
                        flag = True

                    if not flag:
                        continue

                    node1, node2, _ = ref_id.split("_")
                    if node1 == node2:
                        continue

                    high_identity_alignments[query_name].append({
                        'identity': identity,
                        'identity_nogap': identity_nogap,
                        'ref_id': (ref_id, "+") if alignment.is_forward else (ref_id, "-"),
                        'ref_start': ref_start,
                        'ref_end': ref_end,
                        'query_start': query_start,
                        'query_end': query_end,
                        'reverse': alignment.is_reverse,
                        'aln_length_query': query_alignment_end - query_alignment_start,
                        'aln_length_ref': ref_end - ref_start,
                        'length_query': query_len,
                        'length_ref': ref_len,
                        'alignment': alignment
                    })

                    # if "head" in ref_id and alignment.is_forward and ref_start <= 20:
                    #     spanning_reads_per_ref[ref_id].add(query_name)
                    #     print(f"Spanning read {query_name} on ({ref_id}, '+'), ref len {ref_len}, ref span {ref_start}-{ref_end}, query span {query_start}-{query_end}, idt {identity*100:.2f}/{identity_nogap*100:.2f}")
                    # if "head" in ref_id and alignment.is_reverse and ref_len - ref_end <= 20:
                    #     spanning_reads_per_ref[ref_id].add(query_name)
                    #     print(f"Spanning read {query_name} on ({ref_id}, '-'),, ref len {ref_len} ref span {ref_start}-{ref_end}, query span {query_start}-{query_end}, idt {identity*100:.2f}/{identity_nogap*100:.2f}")
                    # if "tail" in ref_id and alignment.is_forward and ref_len - ref_end <= 20:
                    #     spanning_reads_per_ref[ref_id].add(query_name)
                    #     print(f"Spanning read {query_name} on ({ref_id}, '+'), ref len {ref_len}, ref span {ref_start}-{ref_end}, query span {query_start}-{query_end}, idt {identity*100:.2f}/{identity_nogap*100:.2f}")
                    # if "tail" in ref_id and alignment.is_reverse and ref_start <= 20:
                    #     spanning_reads_per_ref[ref_id].add(query_name)
                    #     print(f"Spanning read {query_name} on ({ref_id}, '-'), ref len {ref_len}, ref span {ref_start}-{ref_end}, query span {query_start}-{query_end}, idt {identity*100:.2f}/{identity_nogap*100:.2f}")
    
    # print("Spanning reads per reference (passing all filters):")
    # for ref, reads in sorted(spanning_reads_per_ref.items()):
    #     print(f"{ref}\t{len(reads)} reads")

    return high_identity_alignments

def summarize_full_connected_sequences(high_identity_alignments, ref_fasta, output_result):
    fasta = pysam.FastaFile(ref_fasta)
    ref_pair_connections = defaultdict(list)

    def reverse_complement(seq):
        return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]

    def normalize_pair(r1, s1, r2, s2):
        # Sort by ref name and strand to group symmetric pairs
        s1_r = "-" if s1 == "+" else "+"
        s2_r = "-" if s2 == "+" else "+"
        return min(((r1, s1), (r2, s2)), ((r2, s2_r), (r1, s1_r)))

    for query_name, aln_list in high_identity_alignments.items():
        if len(aln_list) < 2:
            continue

        sorted_alns = sorted(aln_list, key=lambda x: x['query_start'])

        for i in range(len(sorted_alns) - 1):
            for j in range(i+1, len(sorted_alns)):
                aln1 = sorted_alns[i]
                aln2 = sorted_alns[j]

                ref1, strand1 = aln1['ref_id']
                ref2, strand2 = aln2['ref_id']
                identity1 = aln1['identity']
                identity2 = aln2['identity']

                if identity1 < 0.995 and identity2 < 0.995:
                    continue

                if ref1 == ref2:
                    continue

                if ref1[:ref1.rfind("_")] == ref2[:ref2.rfind("_")]:
                    continue

                ref1_len = aln1['length_ref']
                ref2_len = aln2['length_ref']

                # Require ref1 alignment near its end, ref2 near its start
                if not ((ref1_len - aln1['ref_end'] <= 20) and (aln2['ref_start'] <= 20)):
                    continue

                # Extract query bridge
                alignment = aln1['alignment']
                query_seq = alignment.query_sequence
                if query_seq is None:
                    continue

                q_start1 = aln1['query_start']
                q_end1 = aln1['query_end']
                q_start2 = aln2['query_start']
                q_end2 = aln2['query_end']
                if q_end2 < q_end1:
                    continue

                bridge_seq = query_seq[q_start1:q_end2]

                ref1_seq = fasta.fetch(ref1) if strand1 == "+" else reverse_complement(fasta.fetch(ref1))
                ref2_seq = fasta.fetch(ref2) if strand2 == "+" else reverse_complement(fasta.fetch(ref2))
                ref1_seq = ref1_seq[0: aln1['ref_start']]
                ref2_seq = ref2_seq[aln2['ref_end']:]

                connected_seq = ref1_seq + bridge_seq + ref2_seq

                # Normalize for grouping, but store original order for printing
                ref_pair = normalize_pair(ref1, strand1, ref2, strand2)

                gap_size = 0
                if ("tail" in ref1 and strand1 == "+") or ("head" in ref1 and strand1 == "-"):
                    if ("head" in ref2 and strand2 == "+") or ("tail" in ref2 and strand2 == "-"):
                        gap_size = q_start2 - q_end1
                    else:
                        continue
                
                if ("head" in ref1 and strand1 == "+") or ("tail" in ref1 and strand1 == "-"):
                    continue
                
                ref_pair_connections[ref_pair].append({
                    'query_name': query_name,
                    'connected_len': len(connected_seq),
                    'connected_seq': connected_seq,
                    'ref1': (ref1, strand1, aln1['ref_start'], aln1['ref_end']),
                    'ref2': (ref2, strand2, aln2['ref_start'], aln2['ref_end']),
                    'query_span': (q_start1, q_end1, q_start2, q_end2),
                    'gap': gap_size,
                    "idt": (identity1, identity2)
                })

    # filter out conflict pairs
    ref2entries = defaultdict(list)
    ref2ref_pairs = defaultdict(tuple)
    for ref_pair, entries in ref_pair_connections.items():
        contig1_n1, contig1_n2, ori1 = ref_pair[0][0].split("_")
        contig2_n1, contig2_n2, ori2 = ref_pair[1][0].split("_")
        ref_min = min((contig1_n1 + "_" + contig1_n2, contig2_n1 + "_" + contig2_n2), (contig2_n1 + "_" + contig2_n2, contig1_n1 + "_" + contig1_n2))
        if ref_min in ref2entries:
            if len(ref2entries[ref_min]) < len(entries):
                ref2entries[ref_min] = entries
                ref2ref_pairs[ref_min] = ref_pair
            elif len(ref2entries[ref_min]) == len(entries):
                # If same length, keep the one with higher identity
                idt1 = sum(ent['idt'][0] for ent in ref2entries[ref_min]) / len(ref2entries[ref_min])
                idt2 = sum(ent['idt'][0] for ent in entries) / len(entries)
                if idt2 > idt1:
                    ref2entries[ref_min] = entries
                    ref2ref_pairs[ref_min] = ref_pair
        else:
            ref2entries[ref_min] = entries
            ref2ref_pairs[ref_min] = ref_pair

    # Output summary
    with open(output_result, 'w') as out_f:
        for _, ref_pair in ref2ref_pairs.items():
            entries = ref_pair_connections[ref_pair]
            print(f"[Connect] Ref pair: {ref_pair[0]} <--> {ref_pair[1]} ({len(entries)} supporting reads)")

            connected_seq = ""
            max_idt = 0
            contig1_n1, contig1_n2, ori1 = ref_pair[0][0].split("_")
            contig2_n1, contig2_n2, ori2 = ref_pair[1][0].split("_")
            strand1 = ref_pair[0][1]
            strand2 = ref_pair[1][1]

            for ent in entries:
                if ent['idt'][0] > max_idt:
                    connected_seq = ent['connected_seq']
                    max_idt = ent['idt'][0]
                    contig1_n1, contig1_n2, ori1 = ent["ref1"][0].split("_")
                    contig2_n1, contig2_n2, ori2 = ent["ref2"][0].split("_")
                    strand1 = ent["ref1"][1]
                    strand2 = ent["ref2"][1]
            
            if not connected_seq:
                continue

            node1 = ""
            node2 = ""
            if ori1 == "head" and strand1 == "+":
                node1 = contig1_n1
            elif ori1 == "tail" and strand1 == "-":
                node1 = "-" + contig1_n2 if contig1_n2[0] != "-" else contig1_n2[1:]
            elif ori1 == "head" and strand1 == "-":
                node1 = "-" + contig1_n1 if contig1_n1[0] != "-" else contig1_n1[1:]
            elif ori1 == "tail" and strand1 == "+":
                node1 = contig1_n2
            
            if ori2 == "head" and strand2 == "+":
                node2 = contig2_n1
            elif ori2 == "tail" and strand2 == "-":
                node2 = "-" + contig2_n2 if contig2_n2[0] != "-" else contig2_n2[1:]
            elif ori2 == "head" and strand1 == "-":
                node2 = "-" + contig2_n1 if contig2_n1[0] != "-" else contig2_n1[1:]
            elif ori2 == "tail" and strand1 == "+":
                node2 = contig2_n2
            
            print(f"  Nodes: {node1} --> {node2}")
            out_f.write(f"{node1}_{node2}\t{connected_seq}\n")
            print(f'  Nodes: {"-" + node2 if node2[0] != "-" else node2[1:]} --> {"-" + node1 if node1[0] != "-" else node1[1:]}')
            
            sum_gaps = 0
            for ent in entries:
                ref1_info = ent['ref1']
                ref2_info = ent['ref2']
                print(f"  Read: {ent['query_name']}")
                print(f"    Ref1: {ref1_info}, Ref2: {ref2_info}")
                print(f"    Query span: {ent['query_span']}")
                print(f"    Gap: {ent['gap']}")
                print(f"    Connected sequence length: {ent['connected_len']}")
                print(f"    Identities: {ent['idt']}")
                sum_gaps += ent["gap"]
            print (f"[Connect] -- Average gap {sum_gaps/len(entries):.0f} --")
 
def main():
    parser = argparse.ArgumentParser(description="Extract 20kb flanks and align HiFi reads")
    parser.add_argument("fasta", help="FASTA file with contigs")
    parser.add_argument("hifi_reads", help="HiFi reads (FASTQ)")
    parser.add_argument("output_bam", help="Output sorted BAM file")
    parser.add_argument("output_result", help="Output result file")
    parser.add_argument("--flank_size", type=int, default=20000, help="Flank size (default: 20000)")
    parser.add_argument("-t", "--thread", type=int, default=50, help="number of threads to use (default: 50)")
    parser.add_argument("-c", "--compress", help="path to LJA compress")
    parser.add_argument("-i", "--identity", type=float, default=0.9, help="identity threshold for filtering alignments (default: 0.9)")
    args = parser.parse_args()

    flank_fasta = args.fasta[:args.fasta.rfind(".fa")] + ".flanks.fasta"
    if not os.path.isfile(flank_fasta):
        extract_flanks(args.fasta, flank_fasta, args.flank_size)
    if not os.path.isfile(args.output_bam):
        run_minimap2_to_bam(flank_fasta, args.hifi_reads, args.output_bam, args.compress, args.thread)
    
    high_identity_alignments = filter_alignments_with_identity(args.output_bam, 0.9)
    summarize_full_connected_sequences(high_identity_alignments, flank_fasta, args.output_result)

if __name__ == "__main__":
    main()
