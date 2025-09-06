# Parental T2T references (can be .fa or .fa.gz)
PAT_FA=/Poppy/zmzhang/cLJA_Project/Bonobo/genome/mPanPan1.mat.cur.20231122.fasta
MAT_FA=/Poppy/zmzhang/cLJA_Project/Bonobo/genome/mPanPan1.pat.cur.20231122.fasta

# Child phased assembly
CHILD_ASM=/Poppy/zmzhang/cLJA_Project/Bonobo/clja_investigate_scaffolding/5_polishing/assembly.filtered.fasta

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
yak trioeval -t32 pat.yak mat.yak "$CHILD_ASM" > trioeval.clja.txt
