#!/usr/bin/env python
import random
import wotplot
from Bio import SeqIO
from matplotlib import pyplot
import sys

def replace_N(seq: str) -> str:
    """Replace N bases with random A/C/G/T nucleotides."""
    bases = ["A", "C", "G", "T"]
    return "".join(random.choice(bases) if c == "N" else c for c in seq)

# Load and preprocess sequences
sequences_ori = str(list(SeqIO.parse(sys.argv[1], "fasta"))[0].seq).upper()
sequences_rm = str(list(SeqIO.parse(sys.argv[2], "fasta"))[0].seq).upper()

# Replace Ns with random nucleotides
sequences_ori = replace_N(sequences_ori)
sequences_rm = replace_N(sequences_rm)

print(len(sequences_ori), len(sequences_rm))

# Create dot plot
fig, ax = pyplot.subplots(1, 1)
em = wotplot.DotPlotMatrix(sequences_ori, sequences_rm, 200, verbose=True)
wotplot.viz_spy(em, markersize=0.01, ax=ax, title="")

ax.set_xlabel(f"{sys.argv[1][:sys.argv[1].rfind('.')]} ({len(sequences_ori)/1e6:.2f} Mb)")
ax.set_ylabel(f"{sys.argv[2][:sys.argv[2].rfind('.')]} ({len(sequences_rm)/1e6:.2f} Mb)")

fig.set_size_inches(16, 16)
fig.tight_layout()
pyplot.savefig(sys.argv[3], bbox_inches="tight", dpi=1200)
