#!/usr/bin/env bash
<<bws2cv.sh
This script is used to avg over bed for a directory of bigwigs (*.bigwig or *.bigWig). Need bigWigAverageOverBed.
    Usage:
        bws2cv.sh BW_DIR BED OUT_PATH NUM_PROCESS
bws2cv.sh

bw_dir=$1
bed=$2
out_path=$3
nproc=$4

if [[ $bw_dir == "" ]] # no args
then
    echo "Usage: bws2cv.sh BW_DIR BED OUT_PATH"
    exit 1
fi

# generate random tmp dir
tmp_dir=./tmp.`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`
if [[ ! -d $tmp_dir ]]
then
    # dir not exist
    mkdir $tmp_dir
else
    (>&2 echo "Warning: $tmp_dir exist. Overwriting.")
    rm -r $tmp_dir
    mkdir $tmp_dir
fi

# nrow in bed
a=($(wc -l $bed))
nrow=${a[0]}

# iterate through bigwigs
ls $bw_dir/*.bigwig $bw_dir/*.bigWig | xargs -n 1 -P $nproc onebw.sh $bed $tmp_dir $nrow

# Increase maximum number of open file descriptors
ulimit -Sn 5000
ulimit -Hn 5000
echo "binID" > $tmp_dir/binID
cut -f4 $bed | sort -k1,1 >> $tmp_dir/binID
paste $tmp_dir/binID $tmp_dir/*.tab > $out_path
rm -r $tmp_dir

exit 0
