#!/usr/bin/env python
import sys
import edlib
import subprocess
import re
from multiprocessing import Pool
import random

def replace_N(seq: str) -> str:
    """Replace N bases with random A/C/G/T nucleotides."""
    bases = ["A", "C", "G", "T"]
    return "".join(random.choice(bases) if c == "N" else c for c in seq.upper())

def parse_cigar(cigar):
    # Regular expression to extract the operations and lengths from the CIGAR string
    cigar_operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    total_aligned_bases = 0
    exact_matches = 0
    total_aligned_bases_no_gap = 0

    for length, operation in cigar_operations:
        length = int(length)
        
        if operation in ['M', '=', 'X']:  # Matches/mismatches (affect both ref and query)
            if operation == '=' or operation == 'M':
                exact_matches += length  # Exact matches
            total_aligned_bases += length
            total_aligned_bases_no_gap += length
        elif operation == 'D':  # Deletions (affect only the reference)
            total_aligned_bases += length
            if length < 10:
                total_aligned_bases_no_gap += length
        elif operation == 'I':  # Insertions (affect only the query)
            total_aligned_bases += length
            if length < 10:
                total_aligned_bases_no_gap += length

    # Calculate identity as the ratio of exact matches to total aligned bases
    if total_aligned_bases > 0:
        identity = exact_matches / total_aligned_bases
    else:
        identity = 0
    
    if total_aligned_bases_no_gap > 0:
        identity_no_gap = exact_matches / total_aligned_bases_no_gap
    else:
        identity_no_gap = 0

    return identity, identity_no_gap

# Function to read sequences from a FASTA-like file
def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        current_seq_id = None
        current_seq = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if current_seq_id is not None:
                    sequences[current_seq_id] = ''.join(current_seq)
                current_seq_id = line[1:]  # Remove '>' and get the ID
                current_seq = []
            else:
                current_seq.append(line)
        if current_seq_id is not None:
            sequences[current_seq_id] = ''.join(current_seq)
    return sequences

def reverse_complement(dna_seq):
    # Create a dictionary with complementary bases
    dna_seq = dna_seq.upper()
    valid_bases = set('ATCGN')
    if not set(dna_seq.upper()).issubset(valid_bases):
        raise ValueError("Invalid DNA sequence. Only 'A', 'T', 'C', and 'G' are allowed.")
    
    # Create a dictionary with complementary bases
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    
    # Reverse the DNA sequence and get the complementary bases
    reverse_complement_seq = ''.join([complement[base] for base in reversed(dna_seq)])
    
    return reverse_complement_seq

def align_contig_chromosome(args):
    chromosome, contig_id, contig_sequence, chromosome_sequence, strand = args
    if contig_sequence and chromosome_sequence:
        if strand == '-':
            contig_sequence = reverse_complement(contig_sequence)
        print("aligning", chromosome, contig_id, strand, flush=True)
        # Perform global alignment with edlib
        alignment_result = edlib.align(contig_sequence, chromosome_sequence, mode="NW", task="path")
        cigar_string = alignment_result['cigar']
        print("finish", chromosome, contig_id, *parse_cigar(cigar_string), flush=True)
        return chromosome, cigar_string
    else:
        return chromosome, "Sequences not found"

def run_unialigner(args):
    chromosome, contig_id, contig_sequence, chromosome_sequence, strand = args  # Unpack the tuple here
    chromosome = chr2ref[chromosome] if chromosome in chr2ref else "None"
    if contig_sequence and chromosome_sequence:
        if strand == '-':
            contig_sequence = reverse_complement(contig_sequence)
        contig_sequence = replace_N(contig_sequence)
        chromosome_sequence = replace_N(chromosome_sequence)
        print("aligning", chromosome, contig_id, strand, flush=True)
        # Write the contig and chromosome sequences to temporary files
        with open(f"{contig_id}.fasta", 'w') as contig_file, open(f"{chromosome}.fasta", 'w') as chromosome_file:
            contig_file.write(f">{contig_id}\n{contig_sequence}\n")
            chromosome_file.write(f">{chromosome}\n{chromosome_sequence}\n")

        output_dir = f"{chromosome}_{contig_id}_unialigner"
        cigar_file_path = f"{output_dir}/cigar.txt"
        # Run unialigner via subprocess
        cmd = ["/Poppy/zmzhang/software/unialigner_new/tandem_aligner/build/bin/tandem_aligner", "--first", f"{chromosome}.fasta", "--second", f"{contig_id}.fasta", "-o", chromosome + '_' + contig_id+ '_unialigner']
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the command succeeded and retrieve the CIGAR string
        if result.returncode == 0:
            try:
                with open(cigar_file_path, 'r') as cigar_file:
                    cigar_string = cigar_file.read().strip()  # Read the cigar.txt content
            except FileNotFoundError:
                cigar_string = "CIGAR file not found"
        else:
            cigar_string = f"Error: {result.stderr.strip()}"

        # Cleanup: remove the temporary files
        # subprocess.run(["rm", f"{contig_id}.fasta", f"{chromosome}.fasta"])
        # subprocess.run(["rm", "-r", chromosome + '_' + contig_id+ '_unialigner'])
        print("finish", chromosome, contig_id, *parse_cigar(cigar_string), flush=True)
        return chromosome, cigar_string
    else:
        return chromosome, "Sequences not found"

# Step 1: Read the input files from the command line
if len(sys.argv) != 7:
    print("Usage: python script.py <alignment_file> <contig_sequence_file> <chromosome_sequence_file> <processes> <alignment_tool> <output_file>")
    print("Alignment tool options: 'edlib' or 'unialigner'")
    sys.exit(1)

alignment_file = sys.argv[1]
contig_sequence_file = sys.argv[2]
chromosome_sequence_file = sys.argv[3]
processes = int(sys.argv[4])
alignment_tool = sys.argv[5].lower()  # 'edlib' or 'unialigner'
output_file = sys.argv[6]

# Step 2: Parse the alignment file and find the largest contig for each chromosome
global chr2ref
chr2ref = {
    "CM075330.1": "Chr1",
    "CM075331.1": "Chr2",
    "CM075332.1": "Chr3",
    "CM075333.1": "Chr4",
    "CM075334.1": "Chr5",
    "CM075335.1": "Chr6",
    "CM075336.1": "Chr7",
    "CM075337.1": "Chr8",
    "CM075338.1": "Chr9",
    "CM075339.1": "Chr10",
    "CM075340.1": "Chr11",
    "CM075341.1": "Chr12",
    "CM075342.1": "Chr13",
    "CM075343.1": "Chr14",
    "NC_091245.1": "Chr1",
    "NC_091246.1": "Chr2",
    "NC_091247.1": "Chr3",
    "NC_091248.1": "Chr4",
    "NC_091249.1": "Chr5",
    "NC_091250.1": "Chr6",
    "NC_091251.1": "Chr7",
    "NC_091252.1": "Chr8",
    "NC_091253.1": "Chr9",
    "NC_091254.1": "Chr10",
    "NC_091255.1": "Chr11",
    "NC_091256.1": "Chr12",
    "NC_091257.1": "Chr13",
    "NC_091258.1": "Chr14",
    "NC_091259.1": "Chr15",
    "NC_091260.1": "Chr16",
    "NC_091261.1": "Chr17",
    "NC_091262.1": "Chr18",
    "NC_091263.1": "Chr19",
    "NC_091264.1": "Chr20",
    "NC_091265.1": "Chr21",
    "NC_091266.1": "Chr22",
    "NC_091267.1": "Chr23",
    "NC_091268.1": "Chr24",
    "NC_091269.1": "Chr25",
    "NC_091270.1": "Chr26",
    "NC_091271.1": "Y",
    "NC_091727.1": "X"
}
chromosome_contig_dict = {}
with open(alignment_file, 'r') as file:
    for line in file:
        fields = line.split()
        contig_id = fields[0]
        # chromosome = fields[1][:-1]  # Remove the last character (A or B)
        chromosome = fields[1]
        if chromosome[0] == "-":
            chromosome = chromosome[1:]
        chromosome = chr2ref[chromosome] if chromosome in chr2ref else "None"
        contig_length = int(fields[2])

        # Store the contig with the largest length for each chromosome
        if chromosome not in chromosome_contig_dict or contig_length > chromosome_contig_dict[chromosome][1]:
            chromosome_contig_dict[chromosome] = (contig_id, contig_length)

chromosome_strand_dict = {}
with open(alignment_file, 'r') as file:
    for line in file:
        strand = '+'
        fields = line.split()
        contig_id = fields[0]
        # chromosome = fields[1][:-1]  # Remove the last character (A or B)
        chromosome = fields[1]
        if chromosome[0] == "-":
            chromosome = chromosome[1:]
            strand = '-'
        chromosome = chr2ref[chromosome] if chromosome in chr2ref else "None"
        contig_length = int(fields[2])

        # Store the contig with the largest length for each chromosome
        if chromosome not in chromosome_strand_dict or contig_length == chromosome_contig_dict[chromosome][1]:
            chromosome_strand_dict[chromosome] = (contig_id, strand)

for c in chromosome_strand_dict:
    print(c, chromosome_contig_dict[c], chromosome_strand_dict[c])

print("\n=== Largest Contig Length Per Chromosome ===")
print("{:<15} {:<20} {:>10}".format("Chromosome", "Contig ID", "Length"))
for chromosome, (contig_id, length) in sorted(chromosome_contig_dict.items()):
    print("{:<15} {:<20} {:>10}".format(chromosome, contig_id, length))
print("============================================\n", flush=True)

# Step 3: Read contig sequences
contig_sequences = read_fasta(contig_sequence_file)

# Step 4: Read chromosome sequences
chromosome_sequences = read_fasta(chromosome_sequence_file)

# Step 5: Perform global alignments in parallel using ThreadPoolExecutor
args_list = []
ref2chr = {v: k for k, v in chr2ref.items()}
for chromosome, (contig_id, _) in chromosome_contig_dict.items():
    if chromosome not in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr15", "Chr18", "Chr20", "Chr21", "Chr22", "Chr23", "Chr24", "Chr26", "X", "Y"]:
        continue
    
    r_name = ref2chr[chromosome]
    args_list.append((r_name, contig_id, contig_sequences.get(contig_id, ""), chromosome_sequences.get(r_name, ""), chromosome_strand_dict[chromosome][1]))

# Step 6: Run the alignments in parallel using multiprocessing Pool
with Pool(processes=processes) as pool:
    if alignment_tool == 'edlib':
        results = pool.map(align_contig_chromosome, args_list)
    else:
        results = pool.map(run_unialigner, args_list)

# Step 7: Write the results to a file
with open(output_file, 'w') as f:
    for chromosome, cigar in results:
        f.write(f"Chromosome: {chromosome}, CIGAR: {cigar}\n")

print(chromosome_contig_dict)
print(contig_sequences.keys())
print(chromosome_sequences.keys())