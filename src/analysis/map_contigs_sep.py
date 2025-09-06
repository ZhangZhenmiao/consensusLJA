#!/usr/bin/env python3
import os
import sys
from collections import defaultdict
from Bio import SeqIO
import subprocess
from Bio.SeqRecord import SeqRecord

if len(sys.argv) != 5:
    print(f"Usage: {sys.argv[0]} <table.txt> <contigs.fasta> <reference.fasta> <outdir>")
    sys.exit(1)

table_file, contigs_fasta, reference_fasta, outdir = sys.argv[1:]

os.makedirs(outdir, exist_ok=True)

# Step 1: Read table, collect contigs per chromosome
chrom_to_contigs = defaultdict(list)
with open(table_file) as f:
    for line in f:
        if line.strip() == "":
            continue
        cols = line.strip().split("\t")
        contig, chrom = cols[0], cols[1]
        if not chrom.startswith("-"):  # filter
            chrom_to_contigs[chrom].append(contig)

print(f"Found {len(chrom_to_contigs)} chromosomes with contigs")

# Step 2: Load contigs fasta
contig_dict = SeqIO.to_dict(SeqIO.parse(contigs_fasta, "fasta"))

# Step 3: Load reference fasta
ref_dict = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))

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
    "CM075343.1": "Chr14"   
}

new_ref_dict = {}
for old_name, record in ref_dict.items():
    if old_name in id_map:
        new_name = id_map[old_name]
        record.id = new_name              # change FASTA header
        record.name = new_name
        record.description = new_name
        new_ref_dict[new_name] = record
    else:
        # keep old if no mapping
        new_ref_dict[old_name] = record
ref_dict = new_ref_dict

# Step 4: Write per-chromosome contig & reference FASTA files
chrom_fastas = {}
for chrom, contigs in chrom_to_contigs.items():
    contig_out = os.path.join(outdir, f"{chrom}.contigs.fasta")
    ref_out = os.path.join(outdir, f"{chrom}.ref.fasta")

    # Example: when writing per-chromosome FASTA
    values = contig_dict.values()
    if chrom.endswith("M"):
        contigs.extend(chrom_to_contigs[chrom[:-1] + "P"])
    elif chrom.endswith("P"):
        contigs.extend(chrom_to_contigs[chrom[:-1] + "M"])
    contigs = list(set(contigs))  # unique
    with open(contig_out, "w") as out:
        for c in contigs:   # c is the contig ID from table
            found = False
            for rec in values:
                # header like xx_yy -> split by "_"
                if "_" not in rec.id:
                    continue
                xx, yy = rec.id.split("_", 1)
                if yy == c:    # match contig id from table
                    found = True
                    # reverse complement
                    rev = rec.reverse_complement(id=rec.id, description=rec.description)
                    out.write(f">{c}\n")
                    out.write(str(rev.seq) + "\n")
                    break
                elif xx == c:
                    found = True
                    out.write(f">{c}\n")
                    out.write(str(rec.seq) + "\n")
                    break
            if not found:
                print(f"Error: contig {c} not found in contigs fasta")
                exit(1)


    if chrom not in ref_dict:
        print(f"Error: chromosome {chrom} not found in reference fasta")
        exit(1)
    SeqIO.write(ref_dict[chrom], open(ref_out, "w"), "fasta")

    chrom_fastas[chrom] = (contig_out, ref_out)

# Step 5: Run minimap2 in parallel
procs = []
for chrom, (contig_fa, ref_fa) in chrom_fastas.items():
    bam_out = os.path.join(outdir, f"{chrom}.bam")
    print(f"Mapping {chrom} ...")

    cmd1 = ["minimap2", "-ax", "asm20", "-t", "20", ref_fa, contig_fa]
    cmd2 = ["samtools", "view", "-bS", "-"]
    cmd3 = ["samtools", "sort", "-@", "20", "-o", bam_out]

    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(cmd3, stdin=p2.stdout)
    p1.stdout.close()
    p2.stdout.close()
    procs.append((chrom, p3, bam_out))

# Wait for all to finish & index BAM
for chrom, p, bam_out in procs:
    p.communicate()
    subprocess.run(["samtools", "index", bam_out])
    print(f"{chrom} done -> {bam_out}")
