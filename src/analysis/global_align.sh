#!/bin/bash
ref=$1
asm=$2
out=$3
threads=$4
set -e

if [ ! -d $out ]; then
    mkdir $out
fi

ref=$(realpath $ref)
asm=$(realpath $asm)
cd $out

if [ ! -f align.ref.paf ]; then
    minimap2 -x asm20 $ref $asm -t $threads > align.ref.paf
fi
if [ ! -f align.ref.paf.stats ]; then
    /Poppy/zmzhang/Consensus_Assembly/cLJA/src/analysis/global_align_get_ref.py align.ref.paf $asm -o align.ref.paf.stats
fi
# /Poppy/zmzhang/Consensus_Assembly/cLJA/src/analysis/global_align.sheep.py align.ref.paf.stats $asm $ref 10 unialigner results > unialigner.log 2>&1
/Poppy/zmzhang/Consensus_Assembly/cLJA/src/analysis/global_align.py align.ref.paf.stats $asm $ref 10 unialigner results > unialigner.log 2>&1