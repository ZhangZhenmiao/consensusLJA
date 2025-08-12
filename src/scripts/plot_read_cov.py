import pysam
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict

def calculate_identity(alignment):
    cigar = alignment.cigartuples
    if not cigar:
        return 0

    matches = sum(length for op, length in cigar if op == 7)  # '='
    mismatches = sum(length for op, length in cigar if op == 8)  # 'X'
    insertions = sum(length for op, length in cigar if op == 1)
    long_gaps_i = sum(length for op, length in cigar if op == 1 and length >= 10)
    deletions = sum(length for op, length in cigar if op == 2)
    long_gaps_d = sum(length for op, length in cigar if op == 2 and length >= 10)

    aligned_len = matches + mismatches # query_alignment_length is aligned_len + insertions
    total_bases = aligned_len + insertions + deletions
    pid = matches / total_bases if total_bases > 0 else 0.0
    total_bases = total_bases - long_gaps_i - long_gaps_d  # Exclude long gaps from total bases
    pid_nogap = matches / total_bases if total_bases > 0 else 0.0
    return pid


def plot_coverage_along_contig(bam_path, contig_name, identity_threshold=0.9, output_png=None):
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    if contig_name not in bamfile.references:
        raise ValueError(f"Contig '{contig_name}' not found in BAM file.")

    contig_len = bamfile.get_reference_length(contig_name)
    coverage = [0] * contig_len

    for aln in bamfile.fetch(contig_name):
        if aln.is_unmapped:
            continue

        identity = calculate_identity(aln)
        if identity < identity_threshold:
            continue

        for ref_pos in aln.get_reference_positions():
            coverage[ref_pos] += 1

    bamfile.close()

    # Plot coverage along contig
    plt.figure(figsize=(12, 4))
    plt.plot(range(1, contig_len + 1), coverage, color="steelblue", linewidth=0.7)
    plt.title(f"Coverage along {contig_name} (identity ≥ {identity_threshold})")
    plt.xlabel("Position on contig")
    plt.ylabel("Coverage depth")
    plt.grid(True, linestyle="--", alpha=0.4)

    if output_png:
        plt.savefig(output_png, dpi=300)
    else:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot per-base coverage along a given contig with identity filtering.")
    parser.add_argument("bam", help="Input BAM file (sorted and indexed)")
    parser.add_argument("contig", help="Contig name to analyze")
    parser.add_argument("-i", "--identity", type=float, default=0.9, help="Identity threshold (default: 0.9)")
    parser.add_argument("-o", "--output", help="Output PNG file (if not given, show interactively)")
    args = parser.parse_args()

    plot_coverage_along_contig(args.bam, args.contig, args.identity, args.output)