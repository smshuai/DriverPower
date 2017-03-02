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

# check args
if [[ $bed_dir == "" ]] # no args
then
    echo "Usage: beds2cv.sh BED_DIR BED TOT_CG OUT_PATH"
    exit 1
fi
        
# beds - narrowPeak is also a kind of BED
beds="$bed_dir/*.bed $bed_dir/*.narrowPeak $bed_dir/*.bed.gz $bed_dir/*.narrowPeak.gz"
files_path=""
files_head=""

# check bedtools
if hash bedtools 2>/dev/null
then
    version=`bedtools --version`
    (>&2 echo "INFO | Found $version")
else
    (>&2 echo "ERROR | BEDTOOLS not found in path. Aborting.")
    exit 1 
fi
# igenerate random tmp dir
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

# validate bed (expect chrom, start, end, binID)
ncol=`awk '{print NF; exit 0}' $bed`
if [[ $ncol > 4 ]]
then
    (>&2 echo "WARNING | NCOL in $bed is > 4 (chrom, start, end and binID). Truncating.")
    cut -f1-4 $bed > $tmp_dir/bed.tmp
    bed=$tmp_dir/bed.tmp
elif [[ $ncol < 4 ]]
then
    (>&2 echo "ERROR | NCOL in $bed is < 4 (chrom, start, end and binID). Exit.")
    rm -r $tmp_dir
    exit 1
fi

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


(>&2 echo "STEP | Annotating BED with beds found")
echo "chrom start end binID $files_head" | # add chr, start, end, binID
awk 'BEGIN {OFS="\t"} {$1=$1; print $0}' > $tmp_dir/bed_pct.tsv # convert space to tab
# pct of coverage
bedtools annotate -i $bed -files $files_path | sort -k4,4 >> $tmp_dir/bed_pct.tsv
# pct to num of bases
cat $tmp_dir/bed_pct.tsv | awk 'BEGIN{OFS="\t"}
{
    if (NR != 1)
    { # not header. Cols 5-NF are CV pct.
        for ( i = 5; i <= NF; i++)
            $i = $i * ($3 - $2) # pct*(end-start)
        print
    }
    else # header
        print 
}' | cut -f4- > $tmp_dir/bed_cg.tsv # remove chrom, start, end


# test collapse conditions
uniq_binID=`cut -f1 $tmp_dir/bed_cg.tsv | sort -u | wc -l`
all_binID=`cut -f1 $tmp_dir/bed_cg.tsv | sort | wc -l`
if [[ $uniq_binID == $all_binID ]]
then
    (>&2 echo "INFO | binIDs are unique (N=$all_binID). Skip bin collapse.")
else
    (>&2 echo "INFO | Number of unique binID ($uniq_binID) != all binID ($all_binID). Collapsing now.")
    collapse=1
fi

# collpase if necessary
if [[ $collapse == 1 ]]
then
   pivot_cv.py $tmp_dir/bed_cg.tsv $tmp_dir/bed_cg.tsv
fi

(>&2 echo "STEP | Calculating adjusted percent coverage")
# check binID
res=`diff <(cut -f1 $tmp_dir/bed_cg.tsv | sort) <(cut -f1 $tot_cg | sort)`
if [ "$res" != "" ]
then
    (>&2 echo "ERROR | binID in TOT_CG is not the same as BED")
    exit 1
fi

(head -n 1 $tmp_dir/bed_cg.tsv && tail -n +2 $tmp_dir/bed_cg.tsv | sort -k1,1) > $tmp_dir/bed_cg.sorted.tsv
(head -n 1 $tot_cg && tail -n +2 $tot_cg | sort -k1,1) > $tmp_dir/tot_cg.sorted.tsv

join -t $'\t' $tmp_dir/tot_cg.sorted.tsv $tmp_dir/bed_cg.sorted.tsv | # add tot_cg to res.tsv.
awk 'BEGIN{OFS="\t"}
{   # cols 1-2 are binID, tot_cg
    # cols 3:NF are pct cg
    if (NR != 1)
    { # not header
        for (i = 3; i <= NF; i++)
            $i = $2 > 0 ? $i/$2 : 0; # condition ? true-case : false-case
        print
    }
    else
        print
}' > $tmp_dir/bed_cg_adj.tsv

# remove totcg
cut -f1,3- $tmp_dir/bed_cg_adj.tsv > $tmp_dir/res.tsv

mv $tmp_dir/res.tsv $out_path
rm -r $tmp_dir
(>&2 echo "STEP | Beds percent coverage finished")
exit 0

