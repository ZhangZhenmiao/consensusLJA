#!/usr/bin/env python
import wotplot
from Bio import SeqIO
from cigar import Cigar
import edlib
import sys

# Load the sequences from a FASTA file
sequences_ori = str(list(SeqIO.parse(sys.argv[1], "fasta"))[0].seq).upper()
sequences_rm = str(list(SeqIO.parse(sys.argv[2], "fasta"))[0].seq).upper()
from matplotlib import pyplot
print(len(sequences_ori), len(sequences_rm))
ems = []
fig, ax = pyplot.subplots(1, 1)

em = wotplot.DotPlotMatrix(sequences_ori, sequences_rm, 15, verbose=True)
wotplot.viz_spy(
    em, markersize=0.01, ax=ax, title=f""
)
ax.set_xlabel(f"{sys.argv[1][:sys.argv[1].rfind('.')]} ({len(sequences_ori)/1e3:.2f} Kbp)")
ax.set_ylabel(f"{sys.argv[2][:sys.argv[2].rfind('.')]} ({len(sequences_rm)/1e3:.2f} Kbp)")
        
fig.set_size_inches(16, 16)
fig.tight_layout()
pyplot.savefig(sys.argv[3], format="pdf")