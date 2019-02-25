#!/usr/bin/env bash
<<generate_cg.sh
This script is used to generate coverage over bed for DriverPower.
    Prerequisites:
        bedtools - Test on 2.24.0
        R        - Need bioconductor.
    Usage:
        generate_cg.sh FASTA BED CALLABLE
    Args:
        FASTA     - Genome FASTA file 
        BED       - Input BED, use bed6
        CALLABLE  - Input callable regions
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
bedtools intersect -a $BED -b $CALLABLE |
awk '
BEGIN {
    OFS="\t";
    a["chr1"]=249250621;  a["chr2"]=243199373;  a["chr3"]=198022430;
    a["chr4"]=191154276;  a["chr5"]=180915260;  a["chr6"]=171115067;
    a["chr7"]=159138663;  a["chr8"]=146364022;  a["chr9"]=141213431;
    a["chr10"]=135534747; a["chr11"]=135006516; a["chr12"]=133851895;
    a["chr13"]=115169878; a["chr14"]=107349540; a["chr15"]=102531392;
    a["chr16"]=90354753;  a["chr17"]=81195210;  a["chr18"]=78077248;
    a["chr19"]=59128983;  a["chr20"]=63025520;  a["chr21"]=48129895;
    a["chr22"]=51304566;  a["chrX"]=155270560; a["chrY"]=59373566;
    a["chrM"]=16571
} # chrom size from  chrom_size_hg19.tsv
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
