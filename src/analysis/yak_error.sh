ref1=$1
ref2=$2
asm=$3
pref=$4

if [ ! -f mat.yak ]; then
    echo "Building yak k-mer database..."
    yak count -K1.5g -t32 -o mat.yak -k31 $ref1 &
fi

if [ ! -f pat.yak ]; then
    echo "Building yak k-mer database..."
    yak count -K1.5g -t32 -o pat.yak -k31 $ref2 &
fi

wait

yak qv -p -K3.2g -l100k mat.yak $asm > $pref.mat.yak.qv &
yak qv -p -K3.2g -l100k pat.yak $asm > $pref.pat.yak.qv &
wait

# the QV will be the larger one of the two beacause the asm is consensus of the two haplotypes