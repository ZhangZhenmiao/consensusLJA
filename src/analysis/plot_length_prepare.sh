ref=$1
asm=$2
prefix=$3
threads=$4

set -e

# deduplicate
if [ ! -f $prefix.dedup.fa ]; then
    python /scratch/zvz5647/software/consensusLJA/src/analysis/deduplication.py $asm $prefix.dedup.fa > $prefix.dedup.log 2>&1 &
else
    echo "$prefix.dedup.fa already exists, skipping deduplication."
fi

# get referernce stats
if [ ! -f $prefix.paf ]; then
    minimap2 -x asm20 $ref $asm -t $threads > $prefix.ref.paf
else
    echo "$prefix.paf already exists, skipping minimap2 alignment."
fi

/Poppy/zmzhang/Consensus_Assembly/cLJA/src/analysis/plot_length_get_ref.py $prefix.ref.paf -o $prefix.ref_stats

wait