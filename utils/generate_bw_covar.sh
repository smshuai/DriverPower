#!/usr/bin/env bash
bw_list=$1
dirIN=$2
bed=$3
dirOUT=$4

if [[ -f "$bw_list" ]]
then
    echo "BW list found!"
else
    echo "BW list is not found"
    exit 1
fi

cat $bw_list | while read line
do
    bw_dir="$dirIN/$line/"
    if [[ -d "$bw_dir" ]]
    then
        echo "bws2cv.sh $bw_dir $bed $dirOUT/${line}.tsv"
        bws2cv.sh "$bw_dir" "$bed" "$dirOUT/${line}.tsv"
    fi
done
