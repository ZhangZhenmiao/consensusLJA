#!/bin/bash
reads_ori=$1
reads_corr=$2
high_contig=$3
outprefix=$4
compress=$5
correct_reads_script=$6

set -e

if [ ! -f "$outprefix.ori.fasta" ]; then
    $compress --dimer-compress 32,32,1 --reads $reads_ori > $outprefix.ori.fasta
fi

if [ ! -f "$outprefix.bam" ]; then
    minimap2 -t 50 -ax map-hifi --eqx $high_contig $outprefix.ori.fasta | samtools sort -@ 50 -o $outprefix.bam
    samtools index $outprefix.bam
fi

$correct_reads_script -b $outprefix.bam -c $high_contig -i $outprefix.ori.fasta -l $reads_corr -o $outprefix.corrected.fasta

rm $outprefix.ori.fasta