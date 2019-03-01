#!/usr/bin/env bash

bed=$1
tmp_dir=$2
nrow=$3  # nrow in bed
bw=$4

if [[ -f $bw ]]
then
	id=${bw##*/}
	array=(${id//./ })
	id=${array[0]}
    (>&2 echo "Processing $id")
	# Run bigWigAverageOverBed
	bigWigAverageOverBed $bw $bed $tmp_dir/${id}.tab
	# nrow in tab
	b=($(wc -l $tmp_dir/${id}.tab))
	nrow_tab=${b[0]}
    # echo "$nrow and $nrow_tab"
	if [ $nrow == $nrow_tab ]
	then
	    # first column is unique binID, sixth column is mean score
        # sort by binID
        sort -k1,1 $tmp_dir/${id}.tab > $tmp_dir/${id}.tab.sort
        # keep only score
        echo $id > $tmp_dir/${id}.tab # add header first
        cut -f6 $tmp_dir/${id}.tab.sort >> $tmp_dir/${id}.tab # append
        rm $tmp_dir/${id}.tab.sort
	else
		(>&2 echo "Error: table from bigWigAverageOverBed has different number of rows as input BED")
        exit 1
    fi
fi