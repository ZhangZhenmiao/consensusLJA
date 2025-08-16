ref=$1
asm=$2

minimap2 -ax asm20 $ref $asm -t 100 | grep -v '^@' > aln.ref.sam

samtools faidx $ref
cut -f1,2 $ref.fai | awk '{print "@SQ\tSN:"$1"\tLN:"$2}' > aln.ref.header.sam

cat aln.ref.header.sam aln.ref.sam | samtools sort -@ 100 -o aln.ref.sorted.bam