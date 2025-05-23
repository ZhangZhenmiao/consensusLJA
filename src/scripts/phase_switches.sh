#!/bin/bash
# ref=$1
asm=$1
ref_dbg=$2
ref_gaf=$3
out=$4

if [ ! -d $out ]; then
    mkdir $out
fi
# jumboDBG -k 5001 --reads $ref --coverage -o $out/ref.jumboDBG -t 100

# GraphAligner -g $ref_dbg/graph.gfa -f $ref -a $out/ref.gaf -x dbg -t 100
GraphAligner -g $ref_dbg/graph.gfa -f $asm -a $out/asm.gaf -x dbg -t 100

/Poppy/zmzhang/Consensus_Assembly/cLJA/src/scripts/phase_switches.py $ref_gaf $out/asm.gaf $ref_dbg/graph.dot $out/asm.results