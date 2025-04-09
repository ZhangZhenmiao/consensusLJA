reads=$1
asm=$2
outprefix=$3

# compress --dimer-compress 32,32,1 --reads $1 > $outprefix.compressed.fasta
/Poppy/zmzhang/software/anaconda3/envs/clja/bin/minimap2 -ax map-hifi $asm $outprefix.compressed.fasta -t 100 | samtools sort -@ 100 -o $outprefix.bam

python /Poppy/zmzhang/Consensus_Assembly/src/analyze_chimeric.py $outprefix.bam  $reads $outprefix.filtered.fastq