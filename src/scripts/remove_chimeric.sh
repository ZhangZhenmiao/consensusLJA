#!/bin/bash
reads=$1
inprefix=$2
outprefix=$3

compress=$4
analyze_chimeric=$5

if [ ! -f $outprefix.bam ]; then
    $compress --dimer-compress 32,32,1 --reads $reads | minimap2 -ax map-hifi $inprefix.fasta - -t 100 | samtools sort -@ 100 -o $outprefix.bam
fi

$analyze_chimeric $outprefix.bam $inprefix.dot $outprefix.chimeric.txt