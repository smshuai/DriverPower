#!/usr/bin/env bash
<<bws2tbl.sh
This script is used to avg over bed for a directory of bigwigs (*.bigwig or *.bigWig). Need bigWigAverageOverBed.
    Usage:
        bws2tbl.sh BW_DIR BED OUT_PATH
bws2tbl.sh

bw_dir=$1
bed=$2
bws="$bw_dir/*.bigwig $bw_dir/*.bigWig"
out_path=$3
# generate random tmp dir
tmp_dir=$bw_dir/tmp.`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`
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
for bw in $bws
do
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
            cut -f6 $tmp_dir/${id}.tab.sort >> $tmp_dir/${id}.tab # append score
            rm $tmp_dir/${id}.tab.sort
		else
		    (>&2 echo "Error: table from bigWigAverageOverBed has different number of rows as input BED")
            exit 1
        fi
	fi
done

echo "binID" > $tmp_dir/binID
cut -f4 $bed | sort -k1,1 >> $tmp_dir/binID
paste $tmp_dir/binID $tmp_dir/*.tab > $out_path
rm -r $tmp_dir

exit 0
