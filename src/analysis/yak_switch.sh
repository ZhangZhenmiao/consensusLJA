# Parental T2T references (can be .fa or .fa.gz)
PAT_FA=/scratch/zvz5647/hg002/genome/GCA_018852605.3_hg002v1.1.pat_genomic.fna
MAT_FA=/scratch/zvz5647/hg002/genome/GCA_018852615.3_hg002v1.1.mat_genomic.fna

# Child phased assembly
CHILD_ASM=/scratch/zvz5647/hg002/hifiasm/hifiasm.p_ctg.dedup.fa

# 1) Build parental k-mer DBs (k defaults to 31; -b37 is a good general setting)
if [ ! -f "pat.yak" ]; then
    echo "Building paternal k-mer DB..."
    yak count -K1.5g -t32 -o pat.yak "$PAT_FA"
fi
if [ ! -f "mat.yak" ]; then
    echo "Building maternal k-mer DB..."
    yak count -K1.5g -t32 -o mat.yak "$MAT_FA"
fi

# 2) Evaluate child assembly for switches
yak trioeval -t32 pat.yak mat.yak "$CHILD_ASM" > trioeval.hifiasm.txt
