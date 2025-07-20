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

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <fasta1> <fasta2>")
    sys.exit(1)

fasta1 = sys.argv[1]
fasta2 = sys.argv[2]

max_len = 100_000
seq1 = reverse_complement(read_single_fasta(fasta1))[:max_len]
seq2 = read_single_fasta(fasta2)[:max_len]

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
        identity = parse_cigar(cigar_string)
        print(f"CIGAR: {cigar_string}")
        print(f"Identity: {identity:.4f}")
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