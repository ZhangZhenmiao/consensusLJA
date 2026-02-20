#!/bin/sh
set -e

# --- Configuration & Paths ---
GRAPHALIGNER="/Poppy/zmzhang/software/anaconda3/envs/graphaligner/bin/GraphAligner"
GFA="chr21_synteny/bulge_removed.gfa.gfa"
INDEX_DIR="chr21_synteny"
PREFIX="${INDEX_DIR}/graph"
INPUT_READS="consensus_hifiasm_chr21.fa"
OUTPUT_GAF="consensus_hifiasm_chr21.gaf"
CLIPPED_OUT="consensus_hifiasm_chr21.ref.fa"
CORRECTED_OUT="consensus_hifiasm_chr21.corrected.fa"
THREADS=24

# Define alignment-specific variables
MEM_WINDOW="--seeds-mxm-windowsize 5000"

# Ensure the output directory exists
mkdir -p $INDEX_DIR

# --- STEP 1: Indexing Phase ---
# We use an empty file to force GraphAligner to just build the index
echo "Starting Indexing Phase..."
touch empty.fasta

# Build the diploid heuristic cache if not haploid
DIPLOID_FLAGS="--diploid-heuristic 21 31 --diploid-heuristic-cache ${INDEX_DIR}/diploid.index"

if [ ! -f ${PREFIX}.index ]; then
  echo "Building index for the graph..."
  $GRAPHALIGNER -t $THREADS -g $GFA -f empty.fasta -a empty.gaf \
    $DIPLOID_FLAGS \
    --seeds-mxm-cache-prefix $PREFIX \
    --bandwidth 15 \
    --seeds-mxm-length 30 \
    --mem-index-no-wavelet-tree \
    --seeds-mem-count 10000 && touch ${PREFIX}.index
  rm -f empty.gaf empty.fasta
  echo "Indexing complete."
fi

# --- STEP 2: Alignment Phase ---
echo "Starting Alignment Phase..."

$GRAPHALIGNER -t $THREADS -g $GFA -f $INPUT_READS -a $OUTPUT_GAF \
  $DIPLOID_FLAGS \
  --seeds-mxm-cache-prefix $PREFIX \
  $MEM_WINDOW \
  --seeds-mxm-length 30 \
  --seeds-mem-count 10000 \
  --bandwidth 15 \
  --multimap-score-fraction 1 \
  --precise-clipping 0.85 \
  --min-alignment-score 5000 \
  --discard-cigar \
  --clip-ambiguous-ends 100 \
  --overlap-incompatible-cutoff 0.15 \
  --max-trace-count 5 \
  --mem-index-no-wavelet-tree \
  --corrected-clipped-out $CLIPPED_OUT \
  --corrected-out $CORRECTED_OUT
# --multimap-score-fraction 0.99 \

echo "Workflow complete. Alignment saved to: $OUTPUT_GAF"