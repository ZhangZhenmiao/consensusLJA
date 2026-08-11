#!/usr/bin/env python
import argparse
from pathlib import Path
import subprocess
import os
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

def correct_reads_pipeline(
    reads_ori: str,
    reads_corr: str,
    high_contig: str,
    outprefix: str,
    compress: str,
    threads: int = 50
):
    """Integrated read correction pipeline: compress, align, and correct reads."""
    ori_fasta = f"{outprefix}.ori.fasta"
    ori_fasta_tmp = f"{ori_fasta}.tmp"
    bam_file = f"{outprefix}.bam"
    bam_file_tmp = f"{bam_file}.tmp"
    corrected_fasta = f"{outprefix}.corrected.fasta"
    corrected_fasta_tmp = f"{corrected_fasta}.tmp"

    # Step 1: Compress reads if needed
    if not Path(ori_fasta).exists():
        with open(ori_fasta_tmp, "w") as output:
            subprocess.run(
                [compress, "--dimer-compress", "32,32,1", "--reads", reads_ori],
                stdout=output,
                check=True
            )
        os.replace(ori_fasta_tmp, ori_fasta)

    # Step 2: Align reads and sort BAM if needed
    if not Path(bam_file).exists():
        minimap2_cmd = [
            "minimap2", "-t", str(threads), "-ax", "map-hifi", "--eqx",
            high_contig, ori_fasta
        ]
        samtools_sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", bam_file_tmp]

        p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_sort_cmd, stdin=p1.stdout)
        p1.stdout.close()
        sort_returncode = p2.wait()
        minimap2_returncode = p1.wait()
        if minimap2_returncode != 0:
            raise subprocess.CalledProcessError(minimap2_returncode, minimap2_cmd)
        if sort_returncode != 0:
            raise subprocess.CalledProcessError(sort_returncode, samtools_sort_cmd)
        subprocess.run(["samtools", "index", bam_file_tmp], check=True)
        os.replace(bam_file_tmp, bam_file)
        os.replace(f"{bam_file_tmp}.bai", f"{bam_file}.bai")

    # Step 3: Run read correction logic
    print("[CorrectHigh] Loading original read lengths...", flush=True)
    read_lengths = {record.id: len(record.seq) for record in SeqIO.parse(ori_fasta, "fasta")}
    print(f"[CorrectHigh] Loaded {len(read_lengths)} reads from original FASTA", flush=True)

    print("[CorrectHigh] Processing BAM alignments...", flush=True)
    bam = pysam.AlignmentFile(bam_file, "rb")
    ref = pysam.FastaFile(high_contig)
    valid_alignments = {}
    valid_alignments_cov = {}
    all_aligned_reads = set()
    total_reads = valid_count = multi_count = rejected_count = 0

    for read in bam:
        if read.is_unmapped:
            continue
        read_id = read.query_name
        cigar = read.cigartuples
        all_aligned_reads.add(read_id)
        if not cigar:
            continue

        original_len = read_lengths.get(read_id, 0)
        matches = sum(length for op, length in cigar if op == 7)
        mismatches = sum(length for op, length in cigar if op == 8)
        insertions = sum(length for op, length in cigar if op == 1)
        deletions = sum(length for op, length in cigar if op == 2)

        aligned_len = matches + mismatches
        total_bases = aligned_len + insertions + deletions
        pid = matches / total_bases if total_bases > 0 else 0.0
        coverage = read.query_alignment_length / original_len if original_len > 0 else 0.0
        assert(read.query_alignment_length == read.query_alignment_end - read.query_alignment_start)

        total_reads += 1
        if pid > 0.99 and coverage > 0.95:
            if read_id in valid_alignments:
                if coverage > valid_alignments_cov[read_id]:
                    valid_alignments[read_id] = read
                    valid_alignments_cov[read_id] = coverage
                multi_count += 1
            else:
                valid_alignments[read_id] = read
                valid_alignments_cov[read_id] = coverage
                valid_count += 1
        else:
            rejected_count += 1

    low_quality_aligned_reads = {r for r in all_aligned_reads if r not in valid_alignments}
    print(f"[CorrectHigh] Low quality aligned reads: {len(low_quality_aligned_reads)}", flush=True)

    # Write corrected reads
    print("[CorrectHigh] Writing corrected reads to output FASTA...", flush=True)
    written = 0
    with open(corrected_fasta_tmp, "w") as f_out:
        written_reads = set()
        for record in SeqIO.parse(reads_corr, "fasta"):
            read_id = record.id
            if read_id in valid_alignments:
                read = valid_alignments[read_id]
                ref_seq = ref.fetch(read.reference_name, read.reference_start, read.reference_end)
                record.seq = Seq(ref_seq)
            SeqIO.write(record, f_out, "fasta")
            written += 1
            written_reads.add(read_id)

    os.replace(corrected_fasta_tmp, corrected_fasta)

    bam.close()
    ref.close()

    # Step 4: Cleanup
    if Path(ori_fasta).exists():
        os.remove(ori_fasta)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Integrated read correction pipeline")
    parser.add_argument("reads_ori", help="Original reads FASTA")
    parser.add_argument("reads_corr", help="Corrected LJA reads FASTA (input to correction step)")
    parser.add_argument("high_contig", help="High-coverage contigs FASTA")
    parser.add_argument("outprefix", help="Output prefix for intermediate and final files")
    parser.add_argument("compress", help="Compression program path")
    parser.add_argument("--threads", type=int, default=50, help="Number of threads for minimap2 and samtools")
    args = parser.parse_args()

    correct_reads_pipeline(
        args.reads_ori,
        args.reads_corr,
        args.high_contig,
        args.outprefix,
        args.compress,
        threads=args.threads
    )
