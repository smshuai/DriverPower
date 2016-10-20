#!/usr/bin/env bash
<<narrowpks2tbl.sh
This script is used to generate peak density over bed for a directory of narrowpeaks. Need BEDTOOLS.
    Usage:
        
narrowpks2tbl.sh

# args
np_dir=$1
bed=$2
out_path=$3

# narrowpeaks
nps="$np_dir/*.narrowPeak"

# generate random tmp dir
tmp_dir=$np_dir/tmp.`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`
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

# binID
echo 'binID' > $tmp_dir/binID
sort -k4,4 $bed | cut -f4 >> $tmp_dir/binID
# length


# iterate through nps
for np in $nps
do
    id=${np##*/} # keep only file name
    array=(${id//./ }) # convert to bash array, split by "."
    id=${array[0]}
    (>&2 echo "Narrowpeak: procoessing $id")
    bedtools intersect -a $bed -b $np -c > $tmp_dir/${id}.bed
    # sort by binID and keep peak count only 
    echo $id > $tmp_dir/${id}.tab
    sort -k4,4 $tmp_dir/${id}.bed | cut -f5 >> $tmp_dir/${id}.tab
    
done
