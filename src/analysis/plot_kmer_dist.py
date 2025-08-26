from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt

def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def print_top_kmers(kmer_freq, top_n=1000):
    sorted_kmers = sorted(kmer_freq.items(), key=lambda x: x[1], reverse=True)[:top_n]
    print(f"Top {top_n} most frequent 501-mers:\n")
    for kmer, freq in sorted_kmers:
        print(f"{kmer}\t{freq}")

def get_kmer_counts(fasta_file, k=501, flank_size=5000):
    kmer_freq = defaultdict(int)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < 2 * flank_size:
            continue

        prefix = seq[:flank_size]
        suffix = seq[-flank_size:]
        regions = [prefix, suffix]

        for region in regions:
            for i in range(len(region) - k + 1):
                kmer = region[i:i + k]
                rc_kmer = reverse_complement(kmer)
                canonical = min(kmer, rc_kmer)
                kmer_freq[canonical] += 1
    
    print_top_kmers(kmer_freq)

def plot_kmer_distribution(kmer_freq, output_file="501mer_freq_distribution.pdf", top_n=1000):
    # Get top N most frequent k-mers
    top_frequencies = sorted(kmer_freq.values(), reverse=True)[:top_n]

    plt.figure(figsize=(10, 6))
    plt.hist(top_frequencies, bins=50, color='steelblue', edgecolor='black')
    plt.yscale('log')  # Optional: log scale helps if highly skewed
    plt.xlabel('k-mer Frequency (Top {})'.format(top_n))
    plt.ylabel('Count')
    plt.title('Top {} Most Frequent 501-mers'.format(top_n))
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


# === Run the analysis ===
fasta_path = "/Poppy/zmzhang/cLJA_Project/Bonobo/genome/mPanPan1.compressed.fasta"  # change this to your actual file
kmer_counts = get_kmer_counts(fasta_path)
plot_kmer_distribution(kmer_counts)
