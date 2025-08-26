#!/bin/bash
ref=$1
asm=$2
out=$3

if [ ! -d $out ]; then
    mkdir $out
fi

ref=$(realpath $ref)
asm=$(realpath $asm)
cd $out

if [ ! -f align.ref.bam ]; then
    minimap2 -ax asm20 $ref $asm -t 100 | samtools sort -o align.ref.bam
fi
/Poppy/zmzhang/Consensus_Assembly/cLJA/src/analysis/get_reference_custom.py align.ref.bam $asm -o align.ref.bam.stats
/Poppy/zmzhang/Consensus_Assembly/cLJA/src/analysis/global_align.py align.ref.bam.stats $asm $ref 36 unialigner results > unialigner.log 2>&1