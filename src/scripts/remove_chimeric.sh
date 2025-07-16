#!/bin/bash
reads=$1
inprefix=$2
outprefix=$3

compress=$4
analyze_chimeric=$5

if [ ! -f $outprefix.bam ]; then
    rm -f $inprefix.fasta.fai
    samtools faidx $inprefix.fasta
    cut -f1,2 $inprefix.fasta.fai | awk '{print "@SQ\tSN:"$1"\tLN:"$2}' >  $outprefix.sam

    minimap2 -d $inprefix.mmi --split-prefix refsplit $inprefix.fasta
    $compress --dimer-compress 32,32,1 --reads $reads | minimap2 -t 80 -ax map-hifi --eqx $inprefix.mmi - | grep -v '^@' >> $outprefix.sam
    samtools sort -@ 80 $outprefix.sam -o $outprefix.bam; rm $outprefix.sam
fi

$analyze_chimeric $outprefix.bam $inprefix.dot $outprefix.chimeric.txt