#!/usr/bin/env python3
import sys
import subprocess
import re

def read_single_fasta(file_path):
    """Reads a single-sequence FASTA file and returns the sequence."""
    with open(file_path, 'r') as f:
        seq = []
        for line in f:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    import random
    sequence = ''.join(seq)
    # Replace N with random nucleotide
    def replace_N(s):
        nucleotides = ['A', 'C', 'G', 'T']
        return ''.join(random.choice(nucleotides) if c == 'N' else c for c in s)
    return replace_N(sequence)

def reverse_complement(seq):
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def parse_cigar(cigar):
    """Parse CIGAR string and calculate identity."""
    cigar_operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    matches = 0
    total_aligned_bases = 0
    for length, operation in cigar_operations:
        length = int(length)
        if operation in ['M', '=']:
            matches += length
            total_aligned_bases += length
        elif operation in ['X', 'I', 'D']:
            total_aligned_bases += length
    # Use matches divided by the length of the shorter input sequence
    global seq1, seq2
    shorter_len = min(len(seq1), len(seq2))
    identity = matches / shorter_len if shorter_len > 0 else 0
    return identity

def calculate_identities_from_cigar(cigar, gap_threshold=10):
    """
    Calculate identity metrics from a CIGAR string.

    Returns:
    - identity = matches / min(query_len, ref_len)
    - identity_nogap = (matches - mismatches) / (aligned_bases - long gaps)
    """
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    long_gaps = 0

    query_len = 0
    ref_len = 0
    aligned_bases = 0

    cigar_operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    for length_str, op in cigar_operations:
        length = int(length_str)

        if op in ['M', '=', 'X']:
            query_len += length
            ref_len += length
            aligned_bases += length
            if op == '=':
                matches += length
            elif op == 'X':
                mismatches += length
            elif op == 'M':
                matches += length  # assume M = match by default (could overestimate if mismatches not marked as X)

        elif op == 'I':
            query_len += length
            insertions += length
            aligned_bases += length
            if length >= gap_threshold:
                long_gaps += length

        elif op == 'D':
            ref_len += length
            deletions += length
            aligned_bases += length
            if length >= gap_threshold:
                long_gaps += length

        # S/H/N/P ignored for alignment stats

    identity = matches / (matches + mismatches + insertions + deletions)
    identity_nogap = matches / (matches + mismatches + insertions + deletions - long_gaps)

    return identity,identity_nogap

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <fasta1> <fasta2>")
    sys.exit(1)

fasta1 = sys.argv[1]
fasta2 = sys.argv[2]

seq1 = read_single_fasta(fasta1)[:1000000]
seq2 = read_single_fasta(fasta2)[:1000000]

# Write sequences to temp files for unialigner
with open("seq1_tmp.fasta", "w") as f1:
    f1.write(">seq1\n" + seq1 + "\n")
with open("seq2_tmp.fasta", "w") as f2:
    f2.write(">seq2\n" + seq2 + "\n")

output_dir = "unialigner_out"
cigar_file_path = f"{output_dir}/cigar.txt"
cmd = [
    "/Poppy/zmzhang/software/unialigner_new/tandem_aligner/build/bin/tandem_aligner",
    "--first", "seq1_tmp.fasta",
    "--second", "seq2_tmp.fasta",
    "-o", output_dir
]
result = subprocess.run(cmd, capture_output=False, text=True)

if result.returncode == 0:
    try:
        with open(cigar_file_path, "r") as cigar_file:
            cigar_string = cigar_file.read().strip()
        identity,identity_ng = calculate_identities_from_cigar(cigar_string)
        print(f"CIGAR: {cigar_string}")
        print(f"Identity: {identity:.4f}")
        print(f"Identity_ng: {identity_ng:.4f}")
    except FileNotFoundError:
        print("CIGAR file not found. Alignment may have failed.")
else:
    print("Unialigner failed:", result.stderr.strip())

# Cleanup temp files if desired
# import os
# os.remove("seq1_tmp.fasta")
# os.remove("seq2_tmp.fasta")
# import shutil
# shutil.rmtree(output_dir)