#!/bin/bash

set -euo pipefail

# Input arguments
PAT=$1
MAT=$2

if [ ! -f pat_mat.paf ]; then
    echo "Generating PAF file for paternal and maternal assemblies..."
    minimap2 -x asm5 --cs -t 50 $PAT $MAT > pat_mat.paf
else
    echo "PAF file already exists. Skipping generation."
fi

sort -k6,6 -k8,8n pat_mat.paf > pat_mat.srt.paf


paftools.js call pat_mat.srt.paf > pat_mat.vcf

python /Poppy/zmzhang/Consensus_Assembly/cLJA/src/analysis/hetero.py pat_mat.vcf