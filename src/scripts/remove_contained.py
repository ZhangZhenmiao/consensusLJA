import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Remove contigs fully contained in others (allowing mismatches)")
    parser.add_argument("input_fasta", help="Input FASTA file with contigs")
    parser.add_argument("output_fasta", help="Output FASTA file with filtered contigs")
    parser.add_argument("--threads", type=int, default=50, help="Number of threads for minimap2 (default: 8)")
    parser.add_argument("--min_coverage", type=float, default=0.9, help="Minimum coverage to consider contained (default: 0.9)")
    args = parser.parse_args()

    paf_file = Path(args.output_fasta).with_suffix(".paf")

    if not paf_file.exists():

        # Step 1: Run minimap2
        cmd = [
            "minimap2", "-x", "asm20", "-t", str(args.threads),
            args.input_fasta, args.input_fasta
        ]
        with open(paf_file, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)

    # Step 3: Write filtered FASTA
    name2len = {}
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        name2len[record.id] = len(record)

    # Step 2: Parse PAF and detect contained contigs
    contained = {}
    with open(paf_file) as paf:
        for line in paf:
            fields = line.strip().split("\t")
            query, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            target = fields[5]

            if query == target:
                continue

            aligned_len = qend - qstart
            coverage = aligned_len / qlen

            if coverage >= args.min_coverage:
                if query not in contained or coverage > contained[query][1]:
                    if name2len[query] < name2len[target]:
                        contained[query] = (target, coverage)

    # Step 3: Write filtered FASTA
    with open(args.output_fasta, "w") as out_f:
        for record in SeqIO.parse(args.input_fasta, "fasta"):
            if record.id not in contained:
                SeqIO.write(record, out_f, "fasta")

    # Step 4: Report
    print(f"\nContained contigs (≥ {args.min_coverage:.2f} coverage):")
    for query, (target, coverage) in contained.items():
        print(f"{query} len {name2len[query]} is contained in {target} len {name2len[target]} (coverage = {coverage:.3f})")

    print(f"\nRemoved {len(contained)} contained contigs.")
    print(f"Filtered FASTA: {args.output_fasta}")
    print(f"PAF file: {paf_file}")

if __name__ == "__main__":
    main()
