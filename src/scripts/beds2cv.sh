#!/usr/bin/env bash
<<beds2cv.sh
This script is used to generate bed intersected fraction over bed for a directory of beds. Need BEDTOOLS.
    Usage:

beds2cv.sh

# args
bed_dir=$1
bed=$2
tot_cg=$3
out_path=$4

# beds - narrowPeak is also a kind of BED
beds="$bed_dir/*.bed $bed_dir/*.narrowPeak"
files_path=""
files_head=""

# validate files
for file in $beds
do
    if [[ -f $file ]]
    then
        files_path="$files_path $file"
        id=${file##*/} # keep only file name
        array=(${id//./ }) # convert to bash array, split by "."
        id=${array[0]} # file name, used for header
        files_head="$files_head $id"
    fi
done
array=($files_head)
num_file=${#array[@]}
(>&2 echo "INFO | Found $num_file beds:$files_head")

# generate random tmp dir
tmp_dir=$bed_dir/tmp.`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`
if [[ ! -d $tmp_dir ]]
then
    # dir not exist
    mkdir $tmp_dir
else
    (>&2 echo "WARNING | $tmp_dir exist. Overwriting.")
    rm -r $tmp_dir
    mkdir $tmp_dir
fi

(>&2 echo "STEP | Annotating BED with beds found")
echo "chrom start end binID $files_head" | # add chr, start, end, binID
awk 'BEGIN {OFS="\t"} {$1=$1; print $0}' > $tmp_dir/bed_cg.tsv # convert space to tab
# pct of coverage
bedtools annotate -i $bed -files $files_path | sort -k4,4 >> $tmp_dir/bed_cg.tsv

(>&2 echo "STEP | Calculating adjusted percent coverage")
# check binID
res=`diff <(cut -f4 $tmp_dir/bed_cg.tsv) <(cut -f1 $tot_cg)`
if [ "$res" != "" ]
then
    (>&2 echo "ERROR | binID in TOT_CG is not the same as BED")
    exit 1
fi

paste <(cut -f2 $tot_cg) $tmp_dir/bed_cg.tsv | # add tot_cg to res.tsv.
awk 'BEGIN{OFS="\t"}
{   # cols 1-5 are tot_cg, chrom, start, end, binID
    # cols 6:NF are pct cg
    if (NR != 1)
    { # not header
        for (i = 6; i <= NF; i++)
            $i = $1 > 0 ? $i*($4-$3)/$1 : 0; # condition ? true-case : false-case
        print
    }
    else
        print
}' > $tmp_dir/bed_cg_adj.tsv

# remove chrom start end totcg
cut -f5- $tmp_dir/bed_cg_adj.tsv > $tmp_dir/res.tsv

mv $tmp_dir/res.tsv $out_path
rm -r $tmp_dir
(>&2 echo "STEP | Beds percent coverage finished")
exit 0

