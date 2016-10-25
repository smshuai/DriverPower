#!/usr/bin/env bash
bed_list=$1
dirIN=$2
bed=$3
tot_cg=$4
dirOUT=$5

if [[ -f "$bed_list" ]]
then
    echo "Found BED list $bed_list"
else
    echo "BED list $bed_list not found"
    exit 1
fi

cat $bed_list | while read line
do
    bed_dir="$dirIN/$line/"
    if [[ -d "$bed_dir" ]]
    then
        echo "beds2cv.sh $bed_dir $bed $tot_cg $dirOUT/${line}.tsv"
        beds2cv.sh $bed_dir $bed $tot_cg $dirOUT/${line}.tsv
    else
        echo "WARNING: $bed_dir not found"
    fi
done

exit 0
