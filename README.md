
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

Reads can be compressed or uncompressed, and can be provided in FASTQ or FASTA format. Read names must not contain whitespace.

A valid fastq read name example:

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
src/analysis/trim_header reads.fastq.gz -o reads.cleaned.fastq
```

The results will be in <output_folder>/5_polishing/assembly.fasta.
