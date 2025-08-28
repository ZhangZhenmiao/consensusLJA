ref=$1
asm=$2
prefix=$3

# for contig to reference alignment
minimap2 -ax asm20 $ref $asm -t 100 | grep -v '^@' > $prefix.sam

# for contig to contig alignment (investigate cognate contigs)
# minimap2 -ax asm20 --eqx -Y -p 0.1 $ref $asm -t 100 | grep -v '^@' > $prefix.sam

samtools faidx $ref
cut -f1,2 $ref.fai | awk '{print "@SQ\tSN:"$1"\tLN:"$2}' > $prefix.header.sam

cat $prefix.header.sam $prefix.sam | samtools sort -@ 100 -o $prefix.sorted.bam