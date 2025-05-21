#!/bin/bash
ref=$1
asm=$2
out=$3

if [ ! -d $out ]; then
    mkdir $out
fi
jumboDBG -k 5001 --reads $ref --coverage -o $out/ref.jumboDBG -t 100

GraphAligner -g $out/ref.jumboDBG/graph.gfa -f $ref -a $out/ref.gaf -x dbg -t 100
GraphAligner -g $out/ref.jumboDBG/graph.gfa -f $asm -a $out/asm.gaf -x dbg -t 100

/Poppy/zmzhang/Consensus_Assembly/cLJA/src/scripts/phase_switches.py $out/ref.gaf $out/asm.gaf $out/ref.jumboDBG/graph.dot $out/asm.results