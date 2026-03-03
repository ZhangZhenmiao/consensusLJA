
# Mosaic Genome Assembler (MGA)

## Overview

**MGA** is a consensus assemblier developed on [LJA](https://github.com/AntonBankevich/LJA). MGA generates near-complete consensus assemblies using HiFi reads alone.

## Installation

To install **MGA**, follow these steps:

```bash
git clone https://github.com/ZhangZhenmiao/consensusLJA.git
```

MGA depends on the following packages:

- C++ (C++20 support)
- minimap2 (tested v2.21)
- samtools (tested v1.11)
- pysam (tested v0.22.1)
- biopython (tested v1.81)
- numpy (tested v2.2.5)
- CMake (tested v3.22.1)
- GNU Make
- zlib

We provide a conda command to install all dependencies (this creates a new environment named mga; installation typically takes 1–5 minutes):

```bash
conda env create -f requirements.yml
conda activate mga
```


Build MGA:

```bash
cd consensusLJA && cmake . && make -j 8
# MGA executable will be located at bin after building
```

## Usage

Run MGA with the following command:

```bash
usage: MGA --reads=<path_to_reads> --output=<output_folder> [options] ... 
options:
  -r, --reads      path to reads (string)
  -o, --output     the output directory (string)
  -t, --threads    number of threads (int [=50])
  -?, --help       print this message
```

Reads can be compressed or uncompressed, and can be provided in FASTQ or FASTA format. Read names must not contain whitespace. A valid fastq read name example:

```bash
@m84124_230731_175605_s1/251662659/ccs
```
An invalid fastq read name example (contains whitespace and will cause issues for LJA):

```bash
@m84124_230731_175605_s1/251662659/ccs  ML:B:C,255,254,250,247,255,255,120,231,131,218,255,102,89,255,254,249,254,217,255,255,254,244,253,255,251,254,247,230,235,255,207,255,255,238,253   MM:Z:C+m?,30,49,6,29,3,29,39,49,48,51,9,34,13,34,24,33,8,10,8,1,10,39,1,28,27,56,60,8,6,14,1,3,1,0,9;
```

You can use `src/analysis/trim_header` to trim read names if they contain whitespace. For example,

```bash
gcc -o src/analysis/trim_header src/analysis/trim_header.c -lz
# trim_header only supports gzipped fastq files as input
src/analysis/trim_header reads.fastq.gz reads.cleaned.fastq
```

The results will be in <output_folder>/5_polishing/assembly.fasta.

## MGA Parameters (internal)

Like all assembly tools, MGA includes numerous parameters, which can be broadly divided into two categories: assembly parameters and alignment parameters. MGA uses the default **LJA** parameters for constructing assembly graphs, along with several additional parameters:


---

### Parameter Summary

| Category | Parameter | Default Value |
| :--- | :--- | :--- |
| **Graph Cleaning** | `highCovMultiplier` | 10 |
| | `chimera` | 2 |
| **Iterative Simplification** | `maxPathSize` | 8 |
| | `minPI` | 60% |
| | `PI_Tip` | 60% |
| | `TipRepairingPI` | 80% |
| | `SharedLength` | 5 Mb |
| | `SufficientlyShort` | 200 kb |
| | `ShortEdgeLength` | 20 kb |
| **Scaffolding** | `Scaffolding` | (L=5k, k=501); (L=5k, k=301); (L=20k, k=501) |
| | `ScaffoldingPI` | (PI_span_low=90.0%, PI_span_high=99.5%) |
| | `PI_Cognate` | 70% |
| | `Span_Cognate` | 85% |

---

Below we describe how these parameters are used in various modules of MGA.


#### 1. Graph Cleaning
* **Coverage parameter (highCovMultiplier)**: The graph cleaning module automatically estimates the average coverage ($Cov$) from the histogram of read coverage in the graph $LJA_k(Reads)$. The threshold is defined as $highCov = highCovMultiplier \times Cov$, which is used to correct reads in high-coverage regions.
* **Removing chimeric reads (chimera)**: An edge is classified as chimeric if it is supported by fewer than chimera (default = 2) reads.

#### 2. Iterative Graph Simplification
* **Detouring (maxPathSize and minPI)**: A detour is "short" if both of its paths contain at most maxPathSize (default: 8) edges. A detour (that is not a simple bubble) is valid if the percent identity (PI) between its paths is at least minPI (default: 60%).
* **Repairing broken tips (TipRepairingPI)**: A prefix of an edge is its first 1 Mb (or the entire edge if shorter). An out-tip $(v,w)$ is repaired if the percent identity between its prefix and the prefix of its repair edge is at least TipRepairingPI (default: 80%). MGA does not repair tips $\ge 10$ Mb, since such long tips often originate from a different chromosome than their repair edges. The parameters for repairing in-tips are defined similarly.
* **Repairing broken tips plus (SharedLength, SufficientlyShort, and PI_Tip)**: SharedLength (default: 5 Mb) and SufficientlyShort (default: 200 kb) determine if an in-tip is shared by multiple incoming edges. PI_Tip (default: 60%) is the minimum PI to collapse two out-tips, in the operation RepairingTips+.
* **Contracting short edges (ShortEdgeLength)**: An edge is classified as short if its length does not exceed ShortEdgeLength (default: 20 kb).
#### 3. Scaffolding

#### Connecting parameters (ScaffoldingParameters and ScaffoldingPI)
The Connect(G,Reads*) operation forms the set OpenEndsL by generating starting and ending segments of length L for all contigs in the graph multiDBConsensus_k(Reads*). It then identifies weakly-overlapping contigs in OpenEndsL by constructing the graph DB_k(OpenEnds) with small k-mer sizes for several (L, k) combinations. By default, MGA uses the following parameter sets in order: (L=5,000, k=501); (L=5,000, k=301); (L=20,000, k=501). These two combinations are designed for detecting overlaps at most 5,000 bp. MGA detects larger overlaps up to 20,000 bp using the setting (L=20,000, k=501); we select k=501 instead of k=301 for this combination to keep the graph relatively simple. These parameter sets are collectively referred to as ScaffoldingParameters. Each setting is repeated iteratively until no further connections can be made.

The Connect(G,Reads*) module also checks whether any two strings in OpenEndsL (for L=20,000) can be connected by spanning reads. It aligns all reads to OpenEndsL using minimap2 and analyzes all alignments with percent identity ≥ PI_span_low (default 90%) and a span of at least 3 kb. Strings S and T in OpenEndsL are spanned by a read R if 
* R aligns to both a prefix of S (starting coordinate ≤ 20) and a suffix of T (distance from alignment’s end to the end of string T ≤ 20);
* At least one of these alignments has a percent identity ≥ PI_span_high (default 99.5%).

If spanning reads are found for strings S and T in OpenEndsL, MGA merges them into a single string (edge). The pair of parameters (PI_span_low, PI_span_high)  is collectively referred to as ScaffoldingPI.  

#### Deduplication parameters (PI_Cognate and Span_Cognate)
For all contigs spelled by edges in the graph multiDBConsensusk(Reads*), the deduplication module performs all-vs-all alignments  using minimap2 with option “-x asm20” and “-p 0.1”. 

A contig A is classified a cognate contig of contig B if:
* A is aligned to B with percent identity ≥ PI_Cognate (default: 70%);
* The aligned fraction on A is larger than Span_Cognate (default 85%); and
* A is shorter than B.

The edge of A is then removed in multiDBConsensusk(Reads*), and the non-branching paths in the resulting graph are condensed.