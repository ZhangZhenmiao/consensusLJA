#!/usr/bin/env python
import subprocess
import os
from pathlib import Path
import pysam
import numpy as np
from itertools import groupby
import argparse

def calculate_identity(alignment):
    """Calculate percent identity of a read from its CIGAR."""
    cigar = alignment.cigartuples
    if not cigar:
        return 0
    matches = sum(length for op, length in cigar if op == 7)
    mismatches = sum(length for op, length in cigar if op == 8)
    insertions = sum(length for op, length in cigar if op == 1)
    deletions = sum(length for op, length in cigar if op == 2)
    total_bases = matches + mismatches + insertions + deletions
    return matches / total_bases if total_bases > 0 else 0.0

def align_reads_and_detect_chimeric(
    reads: str,
    inprefix: str,
    outprefix: str,
    compress: str,
    dot_file: str,
    threads: int = 80
):
    """
    Full pipeline: compress reads, align with minimap2, sort BAM, detect chimeric reads.
    """
    fasta_file = f"{inprefix}.fasta"
    bam_file = f"{outprefix}.bam"
    sam_file = f"{outprefix}.sam"
    mmi_file = f"{inprefix}.mmi"

    # Step 1: Align if BAM doesn't exist
    if not Path(bam_file).exists():
        # Remove existing FASTA index if any
        fai_file = f"{fasta_file}.fai"
        if Path(fai_file).exists():
            os.remove(fai_file)

        # Index reference FASTA
        subprocess.run(["samtools", "faidx", fasta_file], check=True)

        # Create SAM header
        with open(sam_file, "w") as f_sam:
            subprocess.run(
                f"cut -f1,2 {fai_file} | awk '{{print \"@SQ\\tSN:\"$1\"\\tLN:\"$2}}'",
                shell=True,
                stdout=f_sam,
                check=True
            )

        # Build minimap2 index
        subprocess.run(["minimap2", "-d", mmi_file, "--split-prefix", "refsplit", fasta_file], check=True)

        # Compress reads and align
        compress_proc = subprocess.Popen(
            [compress, "--dimer-compress", "32,32,1", "--reads", reads],
            stdout=subprocess.PIPE
        )

        minimap2_proc = subprocess.Popen(
            ["minimap2", "-t", str(threads), "-ax", "map-hifi", "--eqx", mmi_file, "-"],
            stdin=compress_proc.stdout,
            stdout=subprocess.PIPE
        )
        compress_proc.stdout.close()

        # Filter SAM lines (exclude headers) and append
        with open(sam_file, "a") as f_sam:
            subprocess.run(["grep", "-v", "^@"], stdin=minimap2_proc.stdout, stdout=f_sam, check=True)
        minimap2_proc.wait()

        # Sort SAM to BAM
        subprocess.run(["samtools", "sort", "-@", str(threads), sam_file, "-o", bam_file], check=True)
        os.remove(sam_file)

    # Step 2: Detect chimeric reads
    node2len = {}
    with open(dot_file) as df:
        for line in df:
            if "label" in line and "->" not in line and "+" not in line:
                node = line.split()[0]
                node_len = line[line.find('L')+1: line.find('"', line.find('L')+1)]
                node2len[node] = int(node_len)

    chimeric_edges = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        alignments = {}
        # Collect high-identity alignments
        for read in bam:
            if read.is_unmapped:
                continue
            if calculate_identity(read) > 0.99:
                alignments.setdefault(read.reference_name, []).append({
                    "name": read.query_name,
                    "start": read.reference_start,
                    "end": read.reference_end,
                    "ref_length": bam.get_reference_length(read.reference_name),
                    "aligned": read.query_alignment_length,
                    "idt": calculate_identity(read)
                })

        # Detect chimeric edges
        tolerant_size = 1000
        # print("Contig", "Contig_len", "Node_start", "Node_end", "Read", "Aln_start", "Aln_end", "Aligned_len", "PI", sep="\t", flush=True)
        for r in alignments:
            sorted(alignments[r], key=lambda x: x["start"])
            node1, node2 = r.split('_')
            node1 = node1[:node1.find('.')]
            node2 = node2[:node2.find('.')]
            contig_len = bam.get_reference_length(r)

            # for aln in alignments[r]:
            #     print(r, contig_len, node2len[node1], node2len[node2], aln["name"], aln["start"], aln["end"], aln["aligned"], aln["idt"], sep="\t", flush=True)

            split_coordinate1 = min(node2len[node1], contig_len - node2len[node2])
            split_coordinate2 = max(node2len[node1], contig_len - node2len[node2])
            if split_coordinate1 <= 2*tolerant_size or split_coordinate2 >= contig_len-2*tolerant_size:
                continue
            if node2len[node1] + node2len[node2] - contig_len >= 10000 - tolerant_size:
                continue

            cnt_f = 0
            cnt_r = 0
            cnt_all = 0
            end_counts = {}
            min_internal = tolerant_size
            max_internal = contig_len - tolerant_size
            for aln in alignments[r]:
                if aln["start"] <= max(split_coordinate1 - tolerant_size, 0) and aln["end"] > node2len[node1]:
                    cnt_f += 1
                if aln["start"] < contig_len - node2len[node2] and aln["end"] >= min(split_coordinate2 + tolerant_size, contig_len):
                    cnt_r += 1
                cnt_all += 1
                if min_internal < aln["end"] < max_internal:
                    end_counts[aln["end"]] = end_counts.get(aln["end"], 0) + 1
                
            # Threshold for pileup, can be adjusted
            pileup_threshold = 5
            has_internal_pileup = any(count >= pileup_threshold for count in end_counts.values())

            if (cnt_f <= 0 or cnt_r <= 0) and has_internal_pileup:
                print("[RemoveChimeric]", "contig name", r, "contig len", contig_len, "len 1", node2len[node1], "len 2", node2len[node2], "supporting reads", cnt_f, cnt_r, cnt_all, "is chimeric", flush=True)
                chimeric_edges.extend(r.split('_'))

        # Detect internal chimeric (low coverage regions)
        tolerant_size = 100
        for r in alignments:
            end_counts = {}
            contig_len = bam.get_reference_length(r)
            min_internal = 1000
            max_internal = contig_len - 1000
            coverages = np.zeros(bam.get_reference_length(r), dtype=int)
            for i in alignments[r]:
                if i["end"]-tolerant_size > i["start"]:
                    coverages[i["start"]:i["end"]-tolerant_size] += 1
                if min_internal < i["end"] < max_internal:
                    end_counts[i["end"]] = end_counts.get(i["end"], 0) + 1
            
            pileup_threshold = 3
            has_internal_pileup = any(count >= pileup_threshold for count in end_counts.values())

            potential_pos = []
            if np.average(coverages) >= 5:
                window_size = 20000
                half_window = window_size // 2
                ref_len = bam.get_reference_length(r)
                for i in range(ref_len):
                    if coverages[i] <= 1 and i >= 5001 - 10 and i <= ref_len - 5001 + 10:
                        # Calculate window boundaries
                        left = max(0, i - half_window)
                        right = min(ref_len, i + half_window)
                        window_cov = coverages[left:right]
                        avg_window_cov = np.average(window_cov) if len(window_cov) > 0 else 0
                        # if avg_window_cov >= 5 and has_internal_pileup:
                        if avg_window_cov >= 5:
                            potential_pos.append(i)
            
            filtered_potential_pos = []
            for k, g in groupby(enumerate(potential_pos), lambda x: x[0] - x[1]):
                group = list(map(lambda x: x[1], g))
                if len(group) <= tolerant_size + 10:
                    filtered_potential_pos.extend(group)

            potential_pos = filtered_potential_pos

            # for i in potential_pos:
            #     print("Potential chimeric (low cov) :", i, "in", r)
            
            reads = set()
            for i in alignments[r]:
                for j in potential_pos:
                    if i["start"] <= j-5001+tolerant_size and i["end"] >= j+5001-tolerant_size:
                        reads.add(i["name"])
            
            # print IDs of chimeric reads
            if len(reads) != 0:
                print("[RemoveChimeric]", r, "is chimeric (internal)", flush=True)
                chimeric_edges.extend(r.split('_'))

    return chimeric_edges

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align reads and remove chimeric contigs.")
    parser.add_argument("reads", help="Reads FASTA file")
    parser.add_argument("inprefix", help="Input prefix for reference FASTA")
    parser.add_argument("outprefix", help="Output prefix for BAM and results")
    parser.add_argument("compress", help="Path to compress program")
    parser.add_argument("dot_file", help="DOT file for node lengths")
    parser.add_argument("--threads", type=int, default=50, help="Threads for minimap2 and samtools")
    parser.add_argument("output", help="Output file for chimeric edges")
    args = parser.parse_args()

    if not os.path.isfile(args.output):
        chimeric_edges = align_reads_and_detect_chimeric(
            args.reads,
            args.inprefix,
            args.outprefix,
            args.compress,
            args.dot_file,
            threads=args.threads
        )
        with open(args.output, "w") as fout:
            for e in chimeric_edges:
                fout.write(e + '\n')
