
# consensusLJA

## Overview

**cLJA** is a consensus assemblier developed on [LJA](https://github.com/AntonBankevich/LJA). cLJA generates near-complete consensus assemblies using HiFi reads alone.

## Installation

To install **consensusLJA**, follow these steps:

```bash
git clone https://github.com/ZhangZhenmiao/consensusLJA.git
```

Ensure the following dependencies are installed:

- C++ (C++20 support)
- minimap2 (tested v2.21+)
- samtools (tested v1.11+)
- pysam (tested v0.22.1+)
- CMake (tested v3.16+)
- GNU Make
- zlib

If dependencies are not installed, try below command to install using mamba (faster than conda), or install by yourself:

```bash
mamba env create -f requirements.yml
```


Build cLJA:

```bash
cd cLJA && cmake . && make
```

## Usage

Run cLJA with the following command:

```bash
usage: cLJA --reads=<path_to_reads> --output=<output_folder> [options] ... 
options:
  -r, --reads      path to reads (string)
  -o, --output     the output directory (string)
  -t, --threads    number of threads (int [=50])
  -?, --help       print this message
```

The results will be in <output_folder>/5_polishing/assembly.fasta.