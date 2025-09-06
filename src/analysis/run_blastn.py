#!/usr/bin/env python3
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import argparse
from collections import defaultdict

# -----------------------------
# Helper function to merge intervals
# -----------------------------
def merge_intervals(intervals):
    """Merge overlapping intervals and return total length"""
    if not intervals:
        return 0
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return sum(end - start + 1 for start, end in merged)

# -----------------------------
# Parse command line arguments
# -----------------------------
parser = argparse.ArgumentParser(description="BLAST contigs ≤1Mb against local NT and summarize results")
parser.add_argument("fasta_files", nargs=1, help="Three FASTA files to process")
parser.add_argument("prefixes", nargs=1, help="Prefix for each FASTA file, in same order")
parser.add_argument("-o", "--outdir", default="blast_output", help="Output folder")
parser.add_argument("--db", default="/Poppy/zmzhang/database/nt", help="Local NT database path")
parser.add_argument("--threads", type=int, default=8, help="Number of threads for BLAST")
parser.add_argument("--evalue", type=float, default=1e-10, help="E-value cutoff for BLAST")
args = parser.parse_args()

fasta_files = args.fasta_files
prefixes = args.prefixes
outdir = args.outdir
blast_db = args.db
threads = args.threads
evalue_cutoff = args.evalue

# Create output folder
os.makedirs(outdir, exist_ok=True)

tmp_fasta = os.path.join(outdir, "all_contigs.fa")
blast_out = os.path.join(outdir, "blast_results.tsv")
summary_csv = os.path.join(outdir, "blast_summary.csv")

# -----------------------------
# Step 1: Collect all contigs ≤1Mb into one FASTA with user-defined prefixes
# -----------------------------
with open(tmp_fasta, "w") as out_f:
    for fasta_file, prefix in zip(fasta_files, prefixes):
        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.seq) <= 1_000_000:
                new_id = f"{prefix}_{record.id}"
                out_f.write(f">{new_id}\n{record.seq}\n")

print(f"Collected contigs ≤1Mb into {tmp_fasta}")

# -----------------------------
# Step 2: Run local BLASTN once
# -----------------------------
cmd = [
    "blastn",
    "-task", "megablast",
    "-query", tmp_fasta,
    "-db", blast_db,
    "-out", blast_out,
    "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore",
    "-num_threads", str(threads),
]

print("Running BLASTN...")
subprocess.run(cmd, check=True)
print(f"BLASTN finished. Results saved to {blast_out}")

# -----------------------------
# Step 3: Parse BLAST results, combine hits by species with proper coverage & weighted PI
# -----------------------------
contig_hits = defaultdict(list)

with open(blast_out) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
        qseqid, sseqid, taxid, sciname, ref_len, qlen, qstart, qend, sstart, send, pident, evalue = parts
        qlen = int(qlen)
        ref_len = int(ref_len)
        qstart = int(qstart)
        qend = int(qend)
        pident = float(pident)
        evalue = float(evalue)

        if evalue > evalue_cutoff:
            continue

        contig_hits[qseqid].append({
            "species": sciname,
            "interval": (qstart, qend),
            "ref_len": ref_len,
            "PI": pident
        })

# Combine intervals and compute weighted PI per species
best_hits = {}
for qseqid, hits in contig_hits.items():
    species_intervals = defaultdict(list)
    species_best_hit = {}

    for hit in hits:
        species = hit["species"]
        species_intervals[species].append(hit["interval"])
        # Keep best hit for ref_len
        if species not in species_best_hit or (hit["interval"][1]-hit["interval"][0]) > (species_best_hit[species]["interval"][1]-species_best_hit[species]["interval"][0]):
            species_best_hit[species] = hit

    # Compute combined coverage and weighted PI
    species_cov = {}
    species_weighted_PI = {}
    for species, intervals in species_intervals.items():
        combined_len = merge_intervals(intervals)
        species_cov[species] = combined_len / qlen

        total_weight = 0
        weighted_sum = 0
        for hit in [h for h in hits if h["species"] == species]:
            start, end = hit["interval"]
            length = end - start + 1
            weighted_sum += hit["PI"] * length
            total_weight += length
        species_weighted_PI[species] = weighted_sum / total_weight if total_weight > 0 else 0

    # Select species with largest combined coverage
    best_species = max(species_cov.items(), key=lambda x: x[1])[0]
    best_hit = species_best_hit[best_species]
    best_hits[qseqid] = {
        "Contig file": qseqid.split("_")[0],
        "Contig name": "_".join(qseqid.split("_")[1:]),
        "Contig length": qlen,
        "Scientific name": best_species,
        "Query coverage": round(species_cov[best_species], 3),
        "Ref length": best_hit["ref_len"],
        "PI": round(species_weighted_PI[best_species], 2)
    }

# -----------------------------
# Step 4: Save summary table
# -----------------------------
df = pd.DataFrame(best_hits.values())
df.to_csv(summary_csv, index=False)
print(f"Summary saved to {summary_csv}")