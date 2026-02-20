#!/usr/bin/env python
import sys
import edlib
import subprocess
import re
from multiprocessing import Pool
import random
import os

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
        contig_sequence = replace_N(contig_sequence)
        chromosome_sequence = replace_N(chromosome_sequence)
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
    if contig_sequence and chromosome_sequence:
        output_dir = f"{chromosome}_{contig_id}_unialigner"
        cigar_file_path = f"{output_dir}/cigar.txt"
        if not os.path.exists(output_dir):
            if strand == '-':
                contig_sequence = reverse_complement(contig_sequence)
            contig_sequence = replace_N(contig_sequence)
            chromosome_sequence = replace_N(chromosome_sequence)
            print("aligning", chromosome, contig_id, strand, flush=True)
            # Write the contig and chromosome sequences to temporary files
            with open(f"{contig_id}.fasta", 'w') as contig_file, open(f"{chromosome}.fasta", 'w') as chromosome_file:
                contig_file.write(f">{contig_id}\n{contig_sequence}\n")
                chromosome_file.write(f">{chromosome}\n{chromosome_sequence}\n")
            # Run unialigner via subprocess
            cmd = ["/scratch/zvz5647/software/unialigner/tandem_aligner/src/projects/tandem_aligner/tandem_aligner", "--first", f"{chromosome}.fasta", "--second", f"{contig_id}.fasta", "-o", chromosome + '_' + contig_id+ '_unialigner']
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
        else:
            try:
                with open(cigar_file_path, 'r') as cigar_file:
                    cigar_string = cigar_file.read().strip()  # Read the cigar.txt content
            except FileNotFoundError:
                cigar_string = "CIGAR file not found"

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

id_map = {
    "chr1_mat_hsa1": "1M", 
    "chr2_mat_hsa3": "2M", 
    "chr3_mat_hsa4": "3M", 
    "chr4_mat_hsa5": "4M", 
    "chr5_mat_hsa6": "5M", 
    "chr6_mat_hsa7": "6M", 
    "chr7_mat_hsa8": "7M", 
    "chr8_mat_hsa10": "8M", 
    "chr9_mat_hsa11": "9M", 
    "chr10_mat_hsa12": "10M", 
    "chr11_mat_hsa9": "11M", 
    "chr12_mat_hsa2a": "12M", 
    "chr13_mat_hsa2b": "13M", 
    "chr14_mat_hsa13": "14M", 
    "chr14_mat_hsa13_random_utig4-822": "14M", 
    "chr14_mat_hsa13_random_utig4-823": "14M", 
    "chr14_mat_hsa13_random_utig4-824": "14M", 
    "chr14_mat_hsa13_random_utig4-825": "14M", 
    "chr14_mat_hsa13_random_utig4-826": "14M", 
    "chr14_mat_hsa13_random_utig4-827": "14M", 
    "chr14_mat_hsa13_random_utig4-2241": "14M", 
    "chr14_mat_hsa13_random_utig4-2242": "14M", 
    "chr15_mat_hsa14": "15M", 
    "chr16_mat_hsa15": "16M", 
    "chr17_mat_hsa18": "17M", 
    "chr18_mat_hsa16": "18M", 
    "chr19_mat_hsa17": "19M", 
    "chr20_mat_hsa19": "20M", 
    "chr21_mat_hsa20": "21M", 
    "chr22_mat_hsa21": "22M", 
    "chr23_mat_hsa22": "23M", 
    "chrX_mat_hsaX": "X", 
    "chr1_pat_hsa1": "1P", 
    "chr2_pat_hsa3": "2P", 
    "chr3_pat_hsa4": "3P", 
    "chr4_pat_hsa5": "4P", 
    "chr5_pat_hsa6": "5P", 
    "chr6_pat_hsa7": "6P", 
    "chr7_pat_hsa8": "7P", 
    "chr8_pat_hsa10": "8P", 
    "chr9_pat_hsa11": "9P", 
    "chr10_pat_hsa12": "10P", 
    "chr11_pat_hsa9": "11P", 
    "chr12_pat_hsa2a": "12P", 
    "chr13_pat_hsa2b": "13P", 
    "chr14_pat_hsa13": "14P", 
    "chr15_pat_hsa14": "15P", 
    "chr16_pat_hsa15": "16P", 
    "chr17_pat_hsa18": "17P", 
    "chr18_pat_hsa16": "18P", 
    "chr19_pat_hsa17": "19P", 
    "chr20_pat_hsa19": "20P", 
    "chr21_pat_hsa20": "21P", 
    "chr22_pat_hsa21": "22P", 
    "chr23_pat_hsa22": "23P", 
    "chr1522_pat_hsa1421_random_utig4-95": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-96": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-97": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-99": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-100": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-145": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-147": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-325": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-327": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-839": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-996": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-997": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-1011": "1522P", 
    "chr1522_pat_hsa1421_random_utig4-2055": "1522P", 
    "chrY_pat_hsaY": "Y",
    "chr1_mat": "1M", 
    "chr2_mat": "2M", 
    "chr3_mat": "3M", 
    "chr4_mat": "4M", 
    "chr5_mat": "5M", 
    "chr6_mat": "6M", 
    "chr7_mat": "7M", 
    "chr8_mat": "8M", 
    "chr9_mat": "9M", 
    "chr10_mat": "10M", 
    "chr11_mat": "11M", 
    "chr12_mat": "12M", 
    "chr13_mat": "13M", 
    "chr14_mat": "14M", 
    "chrX_mat": "X",
    "chr1_pat": "1P",
    "chr2_pat": "2P",
    "chr3_pat": "3P",
    "chr4_pat": "4P",
    "chr5_pat": "5P",
    "chr6_pat": "6P",
    "chr7_pat": "7P",
    "chr8_pat": "8P",
    "chr9_pat": "9P",
    "chr10_pat": "10P",
    "chr11_pat": "11P",
    "chr12_pat": "12P",
    "chr13_pat": "13P",
    "chr14_pat": "14P",
    "chrY_pat": "Y",
    "chrM": "mtDNA",
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
    "NC_091727.1": "X",
    "CP139523.2": "1M",
    "CP139519.2": "2M",
    "CP139518.2": "3M",
    "CP139517.2": "4M",
    "CP139516.2": "5M",
    "CP139515.2": "6M",
    "CP139514.2": "7M",
    "CP139513.2": "8M",
    "CP139512.2": "9M",
    "CP139533.2": "10M",
    "CP139532.2": "11M",
    "CP139531.2": "12M",
    "CP139530.2": "13M",
    "CP139529.2": "14M",
    "CP139528.2": "15M",
    "CP139527.2": "16M",
    "CP139526.2": "17M",
    "CP139525.2": "18M",
    "CP139524.2": "19M",
    "CP139522.2": "20M",
    "CP139521.2": "21M",
    "CP139520.2": "22M",
    "CP139511.2": "X",
    "CP139510.1": "mtDNA",
    "CP139546.2": "1P",
    "CP139542.2": "2P",
    "CP139541.2": "3P",
    "CP139540.2": "4P",
    "CP139539.2": "5P",
    "CP139538.2": "6P",
    "CP139537.2": "7P",
    "CP139536.2": "8P",
    "CP139535.2": "9P",
    "CP139556.2": "10P",
    "CP139555.2": "11P",
    "CP139554.2": "12P",
    "CP139553.2": "13P",
    "CP139552.2": "14P",
    "CP139551.2": "15P",
    "CP139550.2": "16P",
    "CP139549.2": "17P",
    "CP139548.2": "18P",
    "CP139547.2": "19P",
    "CP139545.2": "20P",
    "CP139544.2": "21P",
    "CP139543.2": "22P",
    "CP139534.2": "Y"
}

# Step 2: Parse the alignment file and find the largest contig for each chromosome
chromosome_contig_dict = {}
with open(alignment_file, 'r') as file:
    for line in file:
        fields = line.split()
        contig_id = fields[0]
        # chromosome = fields[1][:-1]  # Remove the last character (A or B)
        # chromosome = fields[1][:fields[1].find("_")] if fields[1].find("_") != -1 else fields[1]
        chromosome = fields[1]
        if chromosome[0] == "-":
            chromosome = chromosome[1:]
        chromosome = id_map.get(chromosome, chromosome)
        if "M" in chromosome or "P" in chromosome:
            chromosome = chromosome[:-1]
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
        chromosome = id_map.get(chromosome, chromosome)
        if "M" in chromosome or "P" in chromosome:
            chromosome = chromosome[:-1]
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
chr_names = chromosome_sequences.keys()
chr_names2chr_id = {name: id_map.get(name.split()[0], name) for name in chr_names}
chr_id2chr_name = {v: k for k, v in chr_names2chr_id.items()}
for chromosome, (contig_id, _) in chromosome_contig_dict.items():
    # for giraffe MGA
    # if chromosome not in ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr9", "chr10", "chr12", "chr13", "chrX", "chrY"]:
    #     continue
    # # for bonobo MGA
    # if chromosome not in ["chr3", "chr6", "chr8", "chr9", "chr11", "chr13", "chr16", "chr18", "chr21", "chr22", "chr23"]:
    #     continue
    # for HG002
    if chromosome not in ["2","4", "6", "7", "8", "9", "10", "11", "12", "17", "18", "20", "X"]:
        continue
    if chromosome not in ['X', 'Y']:
        m_name = ""
        p_name = ""
        for name in chr_names:
            m_name = chr_id2chr_name.get(chromosome + "M", "")
            p_name = chr_id2chr_name.get(chromosome + "P", "")
        print(chromosome, contig_id, m_name, p_name, chromosome_strand_dict[chromosome])
        args_list.append((chromosome + "M", contig_id, contig_sequences.get(contig_id, ""), chromosome_sequences.get(m_name, ""), chromosome_strand_dict[chromosome][1]))
        args_list.append((chromosome + "P", contig_id, contig_sequences.get(contig_id, ""), chromosome_sequences.get(p_name, ""), chromosome_strand_dict[chromosome][1]))
    else:
        name = chr_id2chr_name.get(chromosome, "")
        print(chromosome, contig_id, name, chromosome_strand_dict[chromosome])
        args_list.append((chromosome, contig_id, contig_sequences.get(contig_id, ""), chromosome_sequences.get(name, ""), chromosome_strand_dict[chromosome][1]))

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