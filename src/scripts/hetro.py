import sys

def count_aligned_bases(paf_file):
    total_bases = 0
    with open(paf_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                fields = line.strip().split('\t')
                aligned_len = int(fields[10])  # Col 11 is aligned block length
                total_bases += aligned_len
    return total_bases

def count_variants(vcf_file):
    num_variants = 0
    with open(vcf_file) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                num_variants += 1
    return num_variants

def calculate_heterozygosity(num_variants, aligned_bases):
    if aligned_bases == 0:
        return 0.0
    return num_variants / aligned_bases

def main():
    if len(sys.argv) != 3:
        print("Usage: python calculate_heterozygosity.py <paf_file> <vcf_file>")
        sys.exit(1)

    paf_file = sys.argv[1]
    vcf_file = sys.argv[2]

    print("Reading PAF file:", paf_file)
    total_aligned_bases = count_aligned_bases(paf_file)
    print("Total aligned bases:", total_aligned_bases)

    print("Reading VCF file:", vcf_file)
    total_variants = count_variants(vcf_file)
    print("Total variants (SNPs + indels):", total_variants)

    heterozygosity = calculate_heterozygosity(total_variants, total_aligned_bases)
    print(f"Heterozygosity rate: {heterozygosity:.6f} variants per base")

if __name__ == "__main__":
    main()
