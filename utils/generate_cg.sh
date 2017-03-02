#!/usr/bin/env bash
<<generate_cg.sh
This script is used to generate coverage over bed for DriverPower.
    Prerequisites:
        bedtools - Test on 2.24.0
        R        - Need bioconductor.
    Usage:
        generate_cg.sh FASTA BED CALLABLE
    Args:
        FASTA     - Genome FASTA file, use data/hs37d5.fa 
        BED       - Input BED, use 4/5 columns
        CALLABLE  - Input callable regions, use data/blacklist_region/callable50p.bed
        OUT_CG    - Output coverage table. 33 columns. 1st col is binID, 2-33 cols are context based coverage.
        OUT_TOTCG - Output total coverage table. 2 columns. col1 is binID, col2 is the sum of OUT_CG.
generate_cg.sh

FASTA=$1
BED=$2
CALLABLE=$3
OUT_CG=$4
OUT_TOTCG=$5
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
(>&2 echo "INFO | save tmp files to `realpath $tmp_dir`")

# filter by CALLABLE
(>&2 echo "STEP | Intersecting BED with CALLABLE BED")
cat $BED | awk '{ sub(/^chr/, "", $0); print }' |
bedtools intersect -a stdin -b $CALLABLE | # intersect CALLABLE
awk '
BEGIN {
    OFS="\t";
    a[1]=249250621;  a[2]=243199373;  a[3]=198022430;
    a[4]=191154276;  a[5]=180915260;  a[6]=171115067;
    a[7]=159138663;  a[8]=146364022;  a[9]=141213431;
    a[10]=135534747; a[11]=135006516; a[12]=133851895;
    a[13]=115169878; a[14]=107349540; a[15]=102531392;
    a[16]=90354753;  a[17]=81195210;  a[18]=78077248;
    a[19]=59128983;  a[20]=63025520;  a[21]=48129895;
    a[22]=51304566;  a["X"]=155270560; a["Y"]=59373566;
    a["M"]=16571
} # chrom size from  ~/hg19/chrom_size_hg19.tsv
{ # extend bed
    if ($2>0) $2 = $2 - 1;
    if ($3<a[$1]) $3 = $3 + 1;
    print
}' > $tmp_dir/bed.callable.bed
# print extended number of covered bases
awk 'BEGIN {SUM=0} {SUM += $3-$2} END{print "INFO | Extended number of covered bases = " SUM > "/dev/stderr"}' $tmp_dir/bed.callable.bed

# extract seq from the genome
(>&2 echo "STEP | Extracting sequence for BED from genome FASTA")
bedtools getfasta -fi $FASTA -bed $tmp_dir/bed.callable.bed -fo $tmp_dir/tmp.fa

# run fa2cg.R
(>&2 echo "STEP | Making coverage table from extracted sequences")
fa2cg.R $tmp_dir/bed.callable.bed $BED $tmp_dir/tmp.fa $OUT_CG $OUT_TOTCG


rm -r $tmp_dir

(>&2 echo "STEP | Finished")
exit 0
