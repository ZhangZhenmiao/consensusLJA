#!/bin/bash
reads=$1
inprefix=$2
outprefix=$3

compress=$4
analyze_chimeric=$5

if [ ! -f $outprefix.bam ]; then
    minimap2 -d $inprefix.mmi --split-prefix refsplit $inprefix.fasta
    $compress --dimer-compress 32,32,1 --reads $reads | minimap2 -t 100 -ax map-hifi $inprefix.mmi - | samtools sort -@ 100 -o $outprefix.bam
fi

$analyze_chimeric $outprefix.bam $inprefix.dot $outprefix.chimeric.txt