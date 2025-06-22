#!/usr/bin/env python
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

def main(input_bam, contigs_fasta, original_fasta, corrected_fasta, output_fasta):
    print("Loading original read lengths...", flush=True)
    read_lengths = {record.id: len(record.seq) for record in SeqIO.parse(original_fasta, "fasta")}
    print(f"Loaded {len(read_lengths)} reads from original FASTA", flush=True)

    print("Processing BAM alignments...", flush=True)
    # print(f"Read Read_length Aligned_length Mismatches Insertions Deletions PID Aligned_fraction Ref Ref_start Ref_end")
    bam = pysam.AlignmentFile(input_bam, "rb")
    ref = pysam.FastaFile(contigs_fasta)
    valid_alignments = {}
    valid_alignments_cov = {}

    total_reads = 0
    valid_count = 0
    multi_count = 0
    rejected_count = 0

    all_aligned_reads = set()
    for read in bam:
        if read.is_unmapped:
            continue

        read_id = read.query_name
        cigar = read.cigartuples
        all_aligned_reads.add(read_id)

        if not cigar:
            continue

        original_len = read_lengths.get(read_id, 0)
        matches = sum(length for op, length in cigar if op == 7)  # '='
        mismatches = sum(length for op, length in cigar if op == 8)  # 'X'
        insertions = sum(length for op, length in cigar if op == 1)
        deletions = sum(length for op, length in cigar if op == 2)

        aligned_len = matches + mismatches # query_alignment_length is aligned_len + insertions
        total_bases = aligned_len + insertions + deletions
        pid = matches / total_bases if total_bases > 0 else 0.0
        coverage = read.query_alignment_length / original_len if original_len > 0 else 0.0
        assert(read.query_alignment_length == read.query_alignment_end - read.query_alignment_start)

        total_reads += 1
        # print(f"{read_id} {original_len} {read.query_alignment_length} {mismatches} {insertions} {deletions} {pid:.4f} {coverage:.3f} {read.reference_name} {read.reference_start} {read.reference_end}", flush=True)

        if pid > 0.99 and coverage > 0.99:
            if read_id in valid_alignments:
                if coverage > valid_alignments_cov[read_id]:
                    valid_alignments[read_id] = read
                    valid_alignments_cov[read_id] = coverage
                    # print(f"[{read_id}] -> Multi-mapped: repalced (cov {coverage})", flush=True)
                # else:
                #     print(f"[{read_id}] -> Multi-mapped: discard (cov {coverage})", flush=True)
                multi_count += 1
            else:
                valid_alignments[read_id] = read
                valid_alignments_cov[read_id] = coverage
                valid_count += 1
        else:
            rejected_count += 1

    low_quality_aligned_reads = set()
    for r in all_aligned_reads:
        if r not in valid_alignments:
            low_quality_aligned_reads.add(r)
    
    print(f"Low quality aligned reads: {len(low_quality_aligned_reads)}", flush=True)

    print("Writing corrected reads to output FASTA...", flush=True)
    written = 0
    with open(output_fasta, "w") as f_out:
        written_reads = set()
        for record in SeqIO.parse(corrected_fasta, "fasta"):
            read_id = record.id
            # if read_id in multi_mapped:
            #     print(f"[{read_id}] -> Multi-mapped: writing original", flush=True)
            #     SeqIO.write(record, f_out, "fasta")
            if read_id in valid_alignments:
                read = valid_alignments[read_id]
                ref_seq = ref.fetch(read.reference_name, read.reference_start, read.reference_end)
                record.seq = Seq(ref_seq)
                # print(f"[{read_id}] -> Corrected from reference", flush=True)
                SeqIO.write(record, f_out, "fasta")
                written += 1
                written_reads.add(read_id)
            elif read_id not in low_quality_aligned_reads:
                # print(f"[{read_id}] -> No valid alignment: writing original", flush=True)
                SeqIO.write(record, f_out, "fasta")
                written += 1
                written_reads.add(read_id)

        for read_id in valid_alignments:
            if read_id not in written_reads:
                read = valid_alignments[read_id]
                ref_seq = ref.fetch(read.reference_name, read.reference_start, read.reference_end)
                # print(f"[{read_id}] -> Corrected from reference", flush=True)
                f_out.write(f">{read_id}\n{ref_seq}\n")
                written += 1
                written_reads.add(read_id)
        
        for record in SeqIO.parse(original_fasta, "fasta"):
            read_id = record.id
            if read_id not in written_reads:
                SeqIO.write(record, f_out, "fasta")
                written += 1
                written_reads.add(read_id)

    bam.close()
    ref.close()

    print("\nSummary:", flush=True)
    print(f"  Total alignments: {total_reads}", flush=True)
    print(f"  Valid alignments: {valid_count}", flush=True)
    print(f"  Multi-mapped alignments: {multi_count}", flush=True)
    print(f"  Rejected alignments: {rejected_count}", flush=True)
    print(f"  Total reads written to output: {written}", flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct reads based on high-identity alignments using only CIGAR (with =/X).")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file with --eqx CIGAR format")
    parser.add_argument("-c", "--contigs", required=True, help="FASTA file of reference contigs")
    parser.add_argument("-i", "--input_fasta", required=True, help="Original reads in FASTA format")
    parser.add_argument("-l", "--corrected_fasta", required=True, help="Corrected LJA reads in FASTA format")
    parser.add_argument("-o", "--output_fasta", required=True, help="Output corrected reads in FASTA format")
    args = parser.parse_args()

    main(args.bam, args.contigs, args.input_fasta, args.corrected_fasta, args.output_fasta)
