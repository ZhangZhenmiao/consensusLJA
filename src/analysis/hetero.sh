#!/bin/bash

set -euo pipefail

# Input arguments
REF=$1          # Diploid assembly FASTA
READS=$2        # HiFi reads (FASTQ.gz)
PREFIX=$3       # Output prefix

# Tools required: minimap2, samtools, bcftools

if [ ! -f ${PREFIX}.bam.bai ]; then
    # Step 1: Index reference
    minimap2 -d ${PREFIX}.mmi $REF

    # Step 2: Align reads
    minimap2 -ax map-hifi ${PREFIX}.mmi $READS -t 50 | \
        samtools sort -@ 50 -o ${PREFIX}.bam
    samtools index ${PREFIX}.bam
fi

# Step 3: Variant calling
/Poppy/zmzhang/software/bcftools-1.21/bcftools mpileup -Ou -f $REF ${PREFIX}.bam | \
    /Poppy/zmzhang/software/bcftools-1.21/bcftools call -mv -Oz -o ${PREFIX}.vcf.gz
/Poppy/zmzhang/software/bcftools-1.21/bcftools index ${PREFIX}.vcf.gz

# Step 4: Extract heterozygous SNVs
/Poppy/zmzhang/software/bcftools-1.21/bcftools view -v snps -g het ${PREFIX}.vcf.gz > ${PREFIX}_het_snps.vcf

# Step 5: Count heterozygous SNVs
HET_SNPS=$(grep -vc "^#" ${PREFIX}_het_snps.vcf)
echo "Heterozygous SNVs: $HET_SNPS"

# Step 6: Calculate effective genome size (non-N bases)
GENOME_SIZE=$(grep -v "^>" $REF | tr -d 'Nn' | wc -c)
echo "Effective genome size: $GENOME_SIZE bp"

# Step 7: Compute heterozygosity
HETEROZYGOSITY=$(awk -v h=$HET_SNPS -v g=$GENOME_SIZE 'BEGIN {printf "%.6e\n", h/g}')
echo "Heterozygosity (het SNVs per bp): $HETEROZYGOSITY"
