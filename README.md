# consensusLJA

## Overview
**consensusLJA** is an assembler designed to generate high-quality consensus genome assemblies using PacBio high-fidelity (HiFi) reads alone. This tool aims to simplify the assembly process for non-human genomes by producing nearly complete haplotype-mixed consensus assemblies, which are cost-effective and require fewer sequencing technologies compared to diploid assemblies.

## Features
- Generates consensus genome assemblies using only HiFi reads
- Provides nearly complete assemblies with a mosaic of two haplotypes
- Reduces the complexity and cost associated with diploid assembly generation
- Tailored for non-human genomes, though applicable to any genome sequencing projects

## Installation

To install **consensusLJA (cLJA)**, follow the steps below to install the required dependencies, clone the repository, and build the tool.

### Steps:

1. **Install Dependencies**:  
   Install the required dependencies to ensure proper functionality:
   
   - [**minimap2**](https://github.com/lh3/minimap2): v2.21 (For creating BAM files)
   - [**samtools**](https://github.com/samtools/samtools): v1.11 (For processing BAM files)
   - [**pysam**](https://pysam.readthedocs.io/en/latest/installation.html): v0.22.1 (For processing BAM files)

   These dependencies can be installed using Conda, which will create a new environment named `cLJA` (if you already installed these dependencies, skip this step):

   ```bash
   conda env create -f requirements.yml
   ```
2. **Install LJA assembler**:  
   cLJA uses the assembly graph generated by LJA. So LJA assembler (the [experimental](https://github.com/AntonBankevich/LJA/tree/experimental) branch) is required.
   ```bash
   git clone https://github.com/AntonBankevich/LJA.git -b experimental
   cd LJA
   cmake . && make
   ```
3. **Activate the Environment**:  
   Activate the `clja` Conda environment:

   ```bash
   conda activate clja
   ```

4. **Clone the Repository**:  
   Clone the **consensusLJA** repository:

   ```bash
   git clone https://github.com/yourusername/consensusLJA.git
   ```

5. **Build the Tool**:  
   Navigate to the repository and build the tool using `make`:

   ```bash
   cd consensusLJA && make
   ```

6. **Add cLJA to Your PATH**:  
   Add the `cLJA` executable to your environment's `PATH`:

   ```bash
   export PATH="`pwd`":$PATH
   ```

   > _Note_: You can add this line to your `~/.bashrc` for convenience, so that `cLJA` is available in your `PATH` every time you start a new shell session:

   ```bash
   echo -e '\nexport PATH="'$(pwd)'":$PATH' >> ~/.bashrc
   ```
## Usage

With HiFi reads, first assemble using **LJA**. See more details in the [LJA manual](https://github.com/AntonBankevich/LJA/blob/experimental/docs/lja_manual.md):

```bash
lja [options] -o <output_dir> --reads <reads_file> [--reads <reads_file2> ...]
```

The output directory will contain the following structure:

```
<output_dir>/
│
├── 00_CoverageBasedCorrection/
├── 01_TopologyBasedCorrection/
│   ├── final_dbg.dot
│   ├── final_dbg.fasta
├── 02_MDBG/
│   ├── mdbg_edge_seqs.fasta
├── 03_Polishing/
├── assembly.fasta
├── lja.log
├── mdbg.gfa
└── version.txt
```

After assembling with LJA, you can use **cLJA** to process the relevant files (LJA will be integrated to cLJA in later versions).


### cLJA Command Syntax:

```bash
cLJA --dot=<string> --fasta=<string> --multidbg=<string> --output=<string>
```

### Parameters:

- **`-d, --dot <string>`**:  
   This specifies the `graph.dot` file, typically found under `01_TopologyBasedCorrection/final_dbg.dot`. It represents the condensed de Bruijn graph of LJA.  
   *Example*: `--dot=final_dbg.dot, -d final_dbg.dot`

- **`-f, --fasta <string>`**:  
   This expects the `graph.fasta` file, found under `01_TopologyBasedCorrection/final_dbg.fasta`. It contains the sequences of the edges in the graph of LJA.  
   *Example*: `--fasta=final_dbg.fasta, -f final_dbg.fasta`

- **`-m, --multidbg <string>`**:  
   Specifies the `multidbg` file, which contains the edge sequences in the multiplex de Bruijn graph of LJA, found under `02_MDBG/mdbg_edge_seqs.fasta`.  
   *Example*: `--multidbg=mdbg_edge_seqs.fasta, -m mdbg_edge_seqs.fasta`

- **`-o, --output <string>`**:  
   The output directory where the results will be stored. This directory must be new and should not contain any prior data to avoid conflicts or overwrites.  
   *Example*: `--output=/path/to/output_directory, -o /path/to/output_directory`

- **`-?, --help`**:  
   Prints a help message that includes information about all the available options. Use this if you want a quick summary of how to use the tool.  
   *Example*: `--help`

### Example Usage:

In the `example` folder under this repository, there is a small test dataset. The output directory from LJA has already been generated in `LJA_output`. You can run **cLJA** using the command line below:

```bash
cLJA -d LJA_output/01_TopologyBasedCorrection/final_dbg.dot -f LJA_output/01_TopologyBasedCorrection/final_dbg.fasta -m LJA_output/02_MDBG/mdbg_edge_seqs.fasta -o cLJA_output
```

In this example, **cLJA** will process the `final_dbg.dot`, `final_dbg.fasta`, and `mdbg_edge_seqs.fasta` files, and output the results into the `cLJA_output` directory. It will generate the final consensus assembly in the file `cLJA_output/graph.final.fasta`.

If you encounter any issues or errors, please feel free to open an issue on the repository.