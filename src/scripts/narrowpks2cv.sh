#!/usr/bin/env bash
<<narrowpks2tbl.sh
This script is used to generate peak density over bed for a directory of narrowpeaks. Need BEDTOOLS.
    Usage:
        
narrowpks2tbl.sh

# args
np_dir=$1
bed=$2
tot_cg=$3
out_path=$4

# narrowpeaks
nps="$np_dir/*.narrowPeak $np_dir/*.narrowPeak.gz"
files_path=""
files_head=""
# validate files
for np in $nps
do
    if [[ -f $np ]]
    then
        files_path="$files_path $np"
        id=${np##*/} # keep only file name
        array=(${id//./ }) # convert to bash array, split by "."
        id=${array[0]} # file name, used for header
        files_head="$files_head $id"
    fi
done

array=($files_head)
num_np=${#array[@]}
(>&2 echo "INFO | Found $num_np narrowPeaks:$files_head")

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


(>&2 echo "STEP | Annotating BED with narrowpeaks found")
# Add header to result
echo "chrom start end binID $files_head" | # add chr, start, end
awk 'BEGIN {OFS="\t"} {$1=$1; print $0}' > $tmp_dir/np_count.tsv # convert space to tab

bedtools annotate -counts -i $bed -files $files_path | sort -k4,4 >> $tmp_dir/np_count.tsv

(>&2 echo "STEP | Calculating peak density")
# check binID
res=`diff <(cut -f4 $tmp_dir/np_count.tsv) <(cut -f1 $tot_cg)`
if [ "$res" != "" ]
then
    (>&2 echo "ERROR | binID in TOT_CG is not the same as BED")
    exit 1
fi

paste <(cut -f2 $tot_cg) $tmp_dir/np_count.tsv | # add tot_cg to res.tsv. 
awk 'BEGIN{OFS="\t"}
{   # cols 1-5 are tot_cg, chrom, start, end, binID
    # cols 6:NF are peak counts
    if (NR != 1)
    { # not header
        for (i = 6; i <= NF; i++)
            $i = $1 > 0 ? $i/$1 : 0; # condition ? true-case : false-case 
        print
    }
    else
        print
}' > $tmp_dir/np_density.tsv

# remove chrom start end totcg
cut -f5- $tmp_dir/np_density.tsv > $tmp_dir/res.tsv

mv $tmp_dir/res.tsv $out_path
rm -r $tmp_dir

(>&2 echo "STEP | Narrowpeaks finished")
exit 0
