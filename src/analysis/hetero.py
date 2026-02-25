import sys

def parse_paftools_call(call_file, min_mapq=60):
    """
    Parses paftools.js call output to calculate heterozygosity.
    Specifically designed for haploid-to-haploid (e.g., Maternal vs Paternal T2T) comparisons.
    """
    callable_bases = 0
    num_snps = 0
    num_indels = 0

    try:
        with open(call_file, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue

                fields = line.strip().split('\t')
                rec_type = fields[0]

                # R (Region) records define the 'callable' denominator
                if rec_type == "R":
                    # fields[2] = start, fields[3] = end (0-based half-open)
                    start = int(fields[2])
                    end = int(fields[3])
                    callable_bases += (end - start)

                # V (Variant) records
                elif rec_type == "V":
                    # fields[4] = depth (number of query sequences covering this ref base)
                    # fields[5] = MAPQ (mapping quality)
                    # fields[6] = ref allele, fields[7] = alt allele
                    depth = int(fields[4])
                    mapq = int(fields[5])
                    ref_allele = fields[6]
                    alt_allele = fields[7]

                    # 1. Depth == 1 ensures a unique 1-to-1 orthologous mapping
                    # 2. MAPQ filter ensures high-confidence alignments (especially in T2T repeats)
                    if depth == 1 and mapq >= min_mapq:
                        if len(ref_allele) == 1 and len(alt_allele) == 1:
                            num_snps += 1
                        else:
                            num_indels += 1

    except FileNotFoundError:
        print(f"Error: File {call_file} not found.")
        sys.exit(1)

    return callable_bases, num_snps, num_indels


def main():
    if len(sys.argv) < 2:
        print("Usage: python heterozygosity.py <paftools_call_output> [min_mapq]")
        sys.exit(1)

    call_file = sys.argv[1]
    # Default MAPQ to 60 for T2T alignments unless specified otherwise
    min_mapq = int(sys.argv[2]) if len(sys.argv) > 2 else 60

    print(f"--- Analyzing: {call_file} ---")
    print(f"Filtering for MAPQ >= {min_mapq} and Depth == 1")

    total_bases, snps, indels = parse_paftools_call(call_file, min_mapq)

    if total_bases == 0:
        print("Error: No callable bases found. Check if your file contains 'R' records.")
        return

    snp_rate = snps / total_bases
    total_het_rate = (snps + indels) / total_bases

    print(f"\nResults:")
    print(f"{'Callable Bases:':<25} {total_bases:,}")
    print(f"{'Heterozygous SNPs:':<25} {snps:,}")
    print(f"{'Heterozygous Indels:':<25} {indels:,}")
    print("-" * 40)
    print(f"{'SNP Heterozygosity:':<25} {snp_rate:.8f} (per base)")
    print(f"{'Total Heterozygosity:':<25} {total_het_rate:.8f} (per base)")
    print(f"\nInterpretation: 1 SNP every {int(1/snp_rate) if snp_rate > 0 else 0:,} bases.")


if __name__ == "__main__":
    main()